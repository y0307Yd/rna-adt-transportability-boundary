#!/usr/bin/env python
"""Submission-grade measured spatial-protein validation figure."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


MARKERS = ["CD8A", "HLA-DR", "PD-1"]
COLORS = {"CD8A": "#2878B5", "HLA-DR": "#D96B27", "PD-1": "#3A8D6D"}
TIME_COLORS = {"pre": "#7B8794", "post": "#B24C63"}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "axes.linewidth": 0.7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "legend.frameon": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.13, 1.04, label, transform=ax.transAxes, fontsize=8,
            fontweight="bold", va="bottom", ha="left")


def main() -> None:
    base = Path("results/ai_extension")
    output = base / "figures"
    output.mkdir(parents=True, exist_ok=True)
    imc = pd.read_csv(base / "real_spatial_protein/patient_effects.csv")
    imc = imc.loc[imc["endpoint"].eq("tumour")].copy()
    imc_summary = pd.read_csv(base / "real_spatial_protein/primary_summary.csv").set_index("marker")
    paired = pd.read_csv(base / "gse289326_spatial_protein/section_metrics.csv")

    # Exact 183 x 112 mm double-column canvas, expressed in inches.
    fig = plt.figure(figsize=(7.2047, 4.4094), constrained_layout=False)
    grid = fig.add_gridspec(2, 3, left=0.075, right=0.985, bottom=0.13, top=0.94,
                            width_ratios=[1.25, 1, 1], height_ratios=[1, 1],
                            wspace=0.42, hspace=0.52)
    ax_a = fig.add_subplot(grid[:, 0])
    ax_b = fig.add_subplot(grid[0, 1:])
    ax_c = fig.add_subplot(grid[1, 1:])

    # Panel a: all patient-level IMC effects plus locked bootstrap summaries.
    for x, marker in enumerate(MARKERS):
        values = imc.loc[imc["marker"].eq(marker), "effect"].to_numpy(float)
        order = np.argsort(np.argsort(values, kind="stable"), kind="stable")
        jitter = np.linspace(-0.12, 0.12, len(values))[order]
        ax_a.scatter(np.full(len(values), x) + jitter, values, s=15,
                     facecolor="white", edgecolor=COLORS[marker], linewidth=0.8,
                     alpha=0.95, zorder=2)
        row = imc_summary.loc[marker]
        ax_a.errorbar(x, row["median_effect"],
                      yerr=[[row["median_effect"] - row["ci_low"]],
                            [row["ci_high"] - row["median_effect"]]],
                      fmt="o", ms=4.5, color=COLORS[marker], ecolor=COLORS[marker],
                      elinewidth=1.5, capsize=3, capthick=1.1, zorder=4)
        grade = row["grade"].replace("SPATIAL-PROTEIN ", "")
        ax_a.text(x, 0.272, grade, color=COLORS[marker], fontsize=5.5,
                  ha="center", va="top", fontweight="bold")
    ax_a.axhline(0, color="#6F7780", linewidth=0.7, linestyle=(0, (2, 2)), zorder=1)
    ax_a.set_xticks(range(3), MARKERS)
    ax_a.set_ylim(-0.19, 0.29)
    ax_a.set_ylabel("High-minus-low tumour-neighbour fraction")
    ax_a.set_title("Independent imaging mass cytometry", loc="left", fontweight="bold")
    panel_label(ax_a, "a")

    # Panel b: paired patients across treatment contexts, faceted by marker.
    x_positions = {("CD8A", "pre"): 0, ("CD8A", "post"): 0.65,
                   ("HLA-DR", "pre"): 1.55, ("HLA-DR", "post"): 2.20,
                   ("PD-1", "pre"): 3.10, ("PD-1", "post"): 3.75}
    for marker in MARKERS:
        wide = paired.loc[paired["marker"].eq(marker)].pivot(
            index="patient", columns="timepoint", values="rho_neighbour_epcam")
        xp, xq = x_positions[(marker, "pre")], x_positions[(marker, "post")]
        for _, row in wide.iterrows():
            ax_b.plot([xp, xq], [row["pre"], row["post"]], color="#C6CBD0",
                      linewidth=0.7, zorder=1)
        for timepoint, x in (("pre", xp), ("post", xq)):
            values = wide[timepoint].to_numpy(float)
            ax_b.scatter(np.full(6, x), values, s=14, color=TIME_COLORS[timepoint],
                         edgecolor="white", linewidth=0.35, zorder=3)
            ax_b.plot([x - 0.12, x + 0.12], [np.median(values)] * 2,
                      color="#20262D", linewidth=1.5, zorder=4)
    ax_b.axhline(0, color="#6F7780", linewidth=0.7, linestyle=(0, (2, 2)))
    ax_b.set_xticks([0.325, 1.875, 3.425], MARKERS)
    ax_b.set_ylabel("rho with neighbouring EPCAM")
    ax_b.set_title("Paired spatial antibody capture", loc="left", fontweight="bold")
    ax_b.scatter([], [], s=18, color=TIME_COLORS["pre"], label="Pretreatment")
    ax_b.scatter([], [], s=18, color=TIME_COLORS["post"], label="Post-treatment")
    ax_b.legend(loc="upper left", ncol=2, columnspacing=1.1, handletextpad=0.4)
    panel_label(ax_b, "b")

    # Panel c: the prespecified communication boundary using measured proteins.
    pd1 = paired.loc[paired["marker"].eq("PD-1")].pivot(
        index="patient", columns="timepoint", values="rho_neighbour_pdl1")
    for _, row in pd1.iterrows():
        ax_c.plot([0, 1], [row["pre"], row["post"]], color="#C6CBD0",
                  linewidth=0.8, zorder=1)
    for timepoint, x in (("pre", 0), ("post", 1)):
        values = pd1[timepoint].to_numpy(float)
        ax_c.scatter(np.full(6, x), values, s=18, color=TIME_COLORS[timepoint],
                     edgecolor="white", linewidth=0.4, zorder=3)
        ax_c.plot([x - 0.14, x + 0.14], [np.median(values)] * 2,
                  color="#20262D", linewidth=1.6, zorder=4)
    ax_c.axhline(0, color="#6F7780", linewidth=0.7, linestyle=(0, (2, 2)))
    ax_c.set_xlim(-0.45, 1.45)
    ax_c.set_xticks([0, 1], ["Pretreatment", "Post-treatment"])
    ax_c.set_ylabel("PD-1 rho with neighbouring PD-L1")
    ax_c.set_title("Measured-protein communication stress test", loc="left", fontweight="bold")
    ax_c.text(1.43, 0.03, "No positive communication claim", ha="right", va="bottom",
              fontsize=6, color="#59636E")
    panel_label(ax_c, "c")

    for ax in (ax_a, ax_b, ax_c):
        ax.tick_params(length=2.5, width=0.6)
        ax.grid(axis="y", color="#E8EAED", linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)

    stem = output / "Figure4_real_spatial_protein_validation"
    fig.savefig(stem.with_suffix(".svg"))
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".tiff"), dpi=600)
    fig.savefig(stem.with_suffix(".png"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
