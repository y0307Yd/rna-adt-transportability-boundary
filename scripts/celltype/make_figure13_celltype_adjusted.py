"""Figure 13: cell-type composition confounding of pooled Spearman.

Scatter: pooled Spearman (x) vs cell-type-adjusted Spearman (y) per target,
diagonal = no confounding. Points below diagonal = pooled inflated by cell composition.
Color = composition-dominant (ct_baseline > pooled).
Two panels: GSE334503 (discovery) and GSE335494 (validation).
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\celltype_annotation_20260905")
FIG = Path(r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\figures")

d1 = pd.read_csv(BASE / "GSE334503_celltype_adjusted_transport.csv")
d2 = pd.read_csv(BASE / "GSE335494_celltype_adjusted_transport.csv")

fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.6), constrained_layout=True)

for ax, d, title in [(axes[0], d1, "GSE334503 (discovery, 30 donors)"),
                     (axes[1], d2, "GSE335494 (validation, 10 donors)")]:
    comp_dom = d.ct_baseline > d.pooled
    lim = [min(d.pooled.min(), d.ct_adjusted.min()) - 0.08,
           max(d.pooled.max(), d.ct_adjusted.max()) + 0.08]
    ax.plot(lim, lim, "--", color="0.5", lw=1, zorder=1)
    ax.scatter(d.loc[~comp_dom, "pooled"], d.loc[~comp_dom, "ct_adjusted"],
               s=42, c="#2c7fb8", alpha=0.85, edgecolors="white", linewidths=0.5,
               label=f"not comp-dominant ({int((~comp_dom).sum())})", zorder=2)
    ax.scatter(d.loc[comp_dom, "pooled"], d.loc[comp_dom, "ct_adjusted"],
               s=42, c="#d95f02", alpha=0.85, edgecolors="white", linewidths=0.5,
               label=f"composition-dominant ({int(comp_dom.sum())})", zorder=2)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("Pooled Spearman (RNA→ADT)", fontsize=10)
    ax.set_ylabel("Cell-type-adjusted Spearman", fontsize=10)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_aspect("equal", adjustable="box")
    ax.axhline(0, color="0.85", lw=0.7, zorder=0)
    ax.axvline(0, color="0.85", lw=0.7, zorder=0)
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    # annotate medians
    mp = d.pooled.median(); ma = d.ct_adjusted.median()
    ax.text(0.03, 0.95,
            f"pooled median = {mp:.3f}\nadjusted median = {ma:.3f}\nΔ = {mp-ma:+.3f}",
            transform=ax.transAxes, va="top", fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", alpha=0.85))

fig.suptitle("Cell-type composition inflates pooled RNA→ADT transportability",
             fontsize=12.5, fontweight="bold", y=1.04)
fig.savefig(FIG / "Figure13_celltype_adjusted_transport.png", dpi=200,
            bbox_inches="tight", facecolor="white", pad_inches=0.25)
fig.savefig(FIG / "Figure13_celltype_adjusted_transport.svg",
            bbox_inches="tight", facecolor="white", pad_inches=0.25)
print("saved Figure13_celltype_adjusted_transport.png/.svg")

# also emit source data
sd = pd.concat([
    d1.assign(cohort="GSE334503"),
    d2.assign(cohort="GSE335494"),
], ignore_index=True)
sd.to_csv(FIG / "Figure13_source_data.csv", index=False, encoding="utf-8")
print("saved Figure13_source_data.csv")
