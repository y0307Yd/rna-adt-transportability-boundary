#!/usr/bin/env python
"""Locked patient-level IMC neighbourhood analysis for Zenodo 4911135."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


EXPECTED_COUNTS_MD5 = "d64fed7d354d47944bb601d04161e9c3"
MARKERS = {"CD8A": "CD8a", "HLA-DR": "HLA-DR", "PD-1": "PD1"}
IMMUNE_TYPES = {"T_NK", "B_cell", "myeloid", "pDC", "aDC", "plasma_cell"}
SEED = 20260820


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh_adjust(values: pd.Series) -> pd.Series:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return pd.Series(out, index=values.index)


def exact_sign_flip_p(effects: np.ndarray) -> float:
    effects = effects[np.isfinite(effects)]
    observed = abs(float(np.median(effects)))
    null = []
    for signs in itertools.product((-1.0, 1.0), repeat=len(effects)):
        null.append(abs(float(np.median(effects * np.asarray(signs)))))
    null = np.asarray(null)
    return float(np.mean(null >= observed - 1e-15))


def bootstrap_ci(effects: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    n = len(effects)
    sampled = rng.choice(effects, size=(10_000, n), replace=True)
    medians = np.median(sampled, axis=1)
    return tuple(np.quantile(medians, [0.025, 0.975]).tolist())


def read_selected_markers(path: Path) -> pd.DataFrame:
    wanted = set(MARKERS.values())
    rows: dict[str, np.ndarray] = {}
    columns: list[str] | None = None
    with path.open("r", encoding="utf-8") as handle:
        columns = [x.strip('"') for x in handle.readline().rstrip("\n\r").split("\t")]
        for line in handle:
            label, values = line.rstrip("\n\r").split("\t", 1)
            label = label.strip('"')
            if label in wanted:
                rows[label] = np.fromstring(values, sep="\t", dtype=np.float32)
    missing = wanted.difference(rows)
    if missing:
        raise ValueError(f"Missing locked IMC markers: {sorted(missing)}")
    if any(len(v) != len(columns) for v in rows.values()):
        raise ValueError("Count matrix row length does not match cell header")
    return pd.DataFrame(rows, index=columns)


def build_neighbour_outcomes(meta: pd.DataFrame, edges: pd.DataFrame) -> pd.DataFrame:
    eligible = meta.loc[meta["include"].astype(str).str.lower().eq("true")].copy()
    eligible["image"] = eligible["Protein_panel_ImageNumber"].astype(np.int32)
    eligible["object"] = eligible["CellNumber"].astype(np.int32)
    target = eligible[["image", "object", "celltype"]].rename(
        columns={"image": "target_image", "object": "target_object", "celltype": "target_type"}
    )
    source = eligible.loc[eligible["celltype"].eq("T_NK"), ["image", "object", "cellID"]].rename(
        columns={"image": "source_image", "object": "source_object"}
    )
    edges = edges.rename(
        columns={
            "First Image Number": "source_image",
            "First Object Number": "source_object",
            "Second Image Number": "target_image",
            "Second Object Number": "target_object",
        }
    )[["source_image", "source_object", "target_image", "target_object"]]
    for column in edges.columns:
        edges[column] = edges[column].astype(np.int32)
    linked = edges.merge(source, on=["source_image", "source_object"], how="inner", validate="many_to_one")
    linked = linked.merge(target, on=["target_image", "target_object"], how="inner", validate="many_to_one")
    linked["tumour"] = linked["target_type"].eq("tumor").astype(np.float32)
    linked["myeloid"] = linked["target_type"].eq("myeloid").astype(np.float32)
    linked["immune"] = linked["target_type"].isin(IMMUNE_TYPES).astype(np.float32)
    return linked.groupby("cellID", sort=False)[["tumour", "myeloid", "immune"]].mean().reset_index()


def roi_effects(cells: pd.DataFrame, exclude_tls: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    if exclude_tls:
        cells = cells.loc[cells["TLS"].eq("no")].copy()
    effects = []
    continuous = []
    for (sample, roi), group in cells.groupby(["sample", "ROI"], sort=False):
        if len(group) < 40:
            continue
        base = {"sample": sample, "ROI": roi, "IE": group["IE"].iloc[0], "TLS": group["TLS"].iloc[0]}
        for marker in MARKERS:
            values = group[marker].to_numpy(float)
            q1, q3 = np.quantile(values, [0.25, 0.75])
            low = values <= q1
            high = values >= q3
            if low.sum() < 10 or high.sum() < 10 or q1 == q3:
                continue
            for endpoint in ("tumour", "myeloid", "immune"):
                effect = float(group.loc[high, endpoint].mean() - group.loc[low, endpoint].mean())
                effects.append({**base, "marker": marker, "endpoint": endpoint, "effect": effect,
                                "n_low": int(low.sum()), "n_high": int(high.sum())})
            rho = stats.spearmanr(values, group["tumour"].to_numpy(float)).statistic
            continuous.append({**base, "marker": marker, "rho_tumour": float(rho), "n_cells": len(group)})
    return pd.DataFrame(effects), pd.DataFrame(continuous)


def patient_aggregate(roi: pd.DataFrame) -> pd.DataFrame:
    return (roi.groupby(["sample", "IE", "marker", "endpoint"], as_index=False)
            .agg(effect=("effect", "median"), n_rois=("ROI", "nunique")))


def summarize_primary(patient: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    rows = []
    primary = patient.loc[patient["endpoint"].eq("tumour")]
    for marker, group in primary.groupby("marker", sort=False):
        effects = group["effect"].to_numpy(float)
        median = float(np.median(effects))
        direction = np.sign(median)
        agreement = float(np.mean(np.sign(effects) == direction)) if direction else 0.0
        ci_low, ci_high = bootstrap_ci(effects, rng)
        try:
            wilcoxon_p = float(stats.wilcoxon(effects, alternative="two-sided", method="exact").pvalue)
        except ValueError:
            wilcoxon_p = 1.0
        rows.append({"marker": marker, "n_patients": len(effects), "median_effect": median,
                     "ci_low": ci_low, "ci_high": ci_high, "directional_agreement": agreement,
                     "sign_flip_p": exact_sign_flip_p(effects), "wilcoxon_p_descriptive": wilcoxon_p})
    summary = pd.DataFrame(rows)
    summary["sign_flip_q_bh"] = bh_adjust(summary["sign_flip_p"])
    summary["grade"] = "SPATIAL-PROTEIN NO-GO"
    exploratory = ((summary["n_patients"] >= 8) & (summary["median_effect"].ne(0)) &
                   (summary["directional_agreement"] >= 0.60))
    go = (exploratory & (summary["median_effect"].abs() >= 0.02) &
          (summary["sign_flip_q_bh"] <= 0.10))
    summary.loc[exploratory, "grade"] = "SPATIAL-PROTEIN EXPLORATORY"
    summary.loc[go, "grade"] = "SPATIAL-PROTEIN GO"
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("data/spatial_protein/zenodo_4911135"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/ai_extension/real_spatial_protein"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    counts_path = args.data_dir / "Protein_panel_singlecell_counts.txt"
    observed_md5 = md5(counts_path)
    if observed_md5 != EXPECTED_COUNTS_MD5:
        raise ValueError(f"Count matrix MD5 mismatch: {observed_md5}")

    meta = pd.read_csv(args.data_dir / "Protein_panel_singlecell_metadata.csv")
    roi_info = pd.read_csv(args.data_dir / "Protein_panel_ROI_info.csv")
    edges = pd.read_csv(args.data_dir / "Protein_panel_direct_neighbour_cells.csv")
    marker_values = read_selected_markers(counts_path).rename(columns={v: k for k, v in MARKERS.items()})
    if not meta["cellID"].is_unique or not marker_values.index.is_unique:
        raise ValueError("Cell identifiers are not unique")

    outcomes = build_neighbour_outcomes(meta, edges)
    cells = (meta.loc[meta["include"].astype(str).str.lower().eq("true") & meta["celltype"].eq("T_NK")]
             .merge(marker_values, left_on="cellID", right_index=True, how="inner", validate="one_to_one")
             .merge(outcomes, on="cellID", how="inner", validate="one_to_one"))
    cells = cells.merge(roi_info[["sample", "ROI", "TLS"]], on=["sample", "ROI"], how="left",
                        suffixes=("", "_roi"), validate="many_to_one")
    if "TLS_roi" in cells:
        cells["TLS"] = cells["TLS_roi"]

    roi_all, continuous = roi_effects(cells, exclude_tls=False)
    roi_no_tls, _ = roi_effects(cells, exclude_tls=True)
    patient_all = patient_aggregate(roi_all)
    patient_no_tls = patient_aggregate(roi_no_tls)
    primary = summarize_primary(patient_all)
    continuous_patient = (continuous.groupby(["sample", "IE", "marker"], as_index=False)
                          .agg(rho_tumour=("rho_tumour", "median"), n_rois=("ROI", "nunique")))

    roi_all.to_csv(args.output_dir / "roi_effects.csv", index=False)
    patient_all.to_csv(args.output_dir / "patient_effects.csv", index=False)
    primary.to_csv(args.output_dir / "primary_summary.csv", index=False)
    roi_no_tls.to_csv(args.output_dir / "roi_effects_no_tls.csv", index=False)
    patient_no_tls.to_csv(args.output_dir / "patient_effects_no_tls.csv", index=False)
    continuous.to_csv(args.output_dir / "continuous_roi_correlations.csv", index=False)
    continuous_patient.to_csv(args.output_dir / "continuous_patient_correlations.csv", index=False)

    manifest = {
        "seed": SEED,
        "counts_md5": observed_md5,
        "n_metadata_cells": int(len(meta)),
        "n_included_t_nk_with_neighbours": int(len(cells)),
        "n_patients": int(cells["sample"].nunique()),
        "n_rois": int(cells[["sample", "ROI"]].drop_duplicates().shape[0]),
        "eligible_markers": list(MARKERS),
        "structurally_unavailable_locked_markers": ["CD103"],
        "primary_grades": dict(zip(primary["marker"], primary["grade"])),
    }
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(primary.to_string(index=False))
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
