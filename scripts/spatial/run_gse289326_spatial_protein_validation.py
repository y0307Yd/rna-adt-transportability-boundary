#!/usr/bin/env python
"""Locked paired spatial antibody-capture validation in GSE289326."""

from __future__ import annotations

import argparse
import itertools
import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import csc_matrix


SEED = 20260820
N_PERM = 1000
MARKERS = {"CD8A": "CD8A", "HLA-DR": "HLA-DRA", "PD-1": "PDCD1"}
EXPECTED_DIRECTIONS = {"CD8A": 1, "HLA-DR": -1, "PD-1": 1}
REQUIRED = set(MARKERS.values()) | {"CD3E", "EPCAM", "CD274"}
SAMPLES = {
    "GSM8789203": ("P1", "pre"), "GSM8789204": ("P1", "post"),
    "GSM8789205": ("P2", "pre"), "GSM8789206": ("P2", "post"),
    "GSM8789207": ("P3", "pre"), "GSM8789208": ("P3", "post"),
    "GSM8789209": ("P4", "pre"), "GSM8789210": ("P4", "post"),
    "GSM8789212": ("P5", "pre"), "GSM8789213": ("P5", "post"),
    "GSM8789215": ("P6", "pre"), "GSM8789216": ("P6", "post"),
}
HEX_OFFSETS = ((0, -2), (0, 2), (-1, -1), (-1, 1), (1, -1), (1, 1))


def decode(values: np.ndarray) -> np.ndarray:
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in values])


def bh_adjust(values: pd.Series) -> pd.Series:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return pd.Series(out, index=values.index)


def exact_sign_flip_p(values: np.ndarray) -> float:
    observed = abs(float(np.mean(values)))
    null = [abs(float(np.mean(values * np.asarray(signs))))
            for signs in itertools.product((-1.0, 1.0), repeat=len(values))]
    return float(np.mean(np.asarray(null) >= observed - 1e-15))


def find_file(data_dir: Path, gsm: str, suffix: str) -> Path:
    hits = list(data_dir.glob(f"{gsm}_*{suffix}"))
    if len(hits) != 1:
        raise FileNotFoundError(f"Expected one {gsm} *{suffix}, found {len(hits)}")
    return hits[0]


def load_section(data_dir: Path, gsm: str) -> pd.DataFrame:
    h5_path = find_file(data_dir, gsm, "_filtered_feature_bc_matrix.h5")
    position_path = find_file(data_dir, gsm, "_tissue_positions.csv.gz")
    with h5py.File(h5_path, "r") as handle:
        matrix = handle["matrix"]
        names = decode(matrix["features/name"][:])
        secondary = decode(matrix["features/secondary_name"][:])
        types = decode(matrix["features/feature_type"][:])
        isotype = decode(matrix["features/isotype_control"][:])
        barcodes = decode(matrix["barcodes"][:])
        sparse = csc_matrix((matrix["data"][:], matrix["indices"][:], matrix["indptr"][:]),
                            shape=tuple(matrix["shape"][:]))
        antibody_idx = np.flatnonzero((types == "Antibody Capture") & (isotype != "TRUE"))
        antibody_names = names[antibody_idx]
        antibody_secondary = secondary[antibody_idx]
        duplicated = pd.Series(antibody_names).duplicated(keep=False).to_numpy()
        antibody_names = np.asarray([
            f"{name}__{second}" if is_duplicate else name
            for name, second, is_duplicate in zip(antibody_names, antibody_secondary, duplicated)
        ])
        missing = REQUIRED.difference(antibody_names)
        if missing:
            raise ValueError(f"Missing required antibodies in {gsm}: {sorted(missing)}")
        raw = sparse[antibody_idx, :].toarray().astype(np.float64).T

    positions = pd.read_csv(position_path)
    positions = positions.loc[positions["in_tissue"].eq(1)].copy()
    lookup = pd.Series(np.arange(len(barcodes)), index=barcodes)
    if not positions["barcode"].isin(lookup.index).all():
        raise ValueError(f"Position barcode mismatch in {gsm}")
    raw = raw[lookup.loc[positions["barcode"]].to_numpy(), :]

    low = np.quantile(raw, 0.05, axis=0)
    high = np.quantile(raw, 0.95, axis=0)
    winsor = np.clip(raw, low, high)
    log_values = np.log1p(winsor)
    clr = log_values - log_values.mean(axis=1, keepdims=True)
    proteins = pd.DataFrame(clr, columns=antibody_names, index=positions.index)
    return pd.concat([positions.reset_index(drop=True), proteins.reset_index(drop=True)], axis=1)


def neighbour_mean(values: np.ndarray, coordinates: list[tuple[int, int]]) -> np.ndarray:
    coord_to_idx = {coord: i for i, coord in enumerate(coordinates)}
    output = np.full(len(coordinates), np.nan)
    for i, (row, col) in enumerate(coordinates):
        indices = [coord_to_idx[(row + dr, col + dc)] for dr, dc in HEX_OFFSETS
                   if (row + dr, col + dc) in coord_to_idx]
        if indices:
            output[i] = float(np.mean(values[indices]))
    return output


def section_metrics(section: pd.DataFrame, gsm: str, rng: np.random.Generator) -> list[dict]:
    if len(section) < 100:
        raise ValueError(f"{gsm} has fewer than 100 in-tissue spots")
    coordinates = list(zip(section["array_row"].astype(int), section["array_col"].astype(int)))
    neighbour_epcam = neighbour_mean(section["EPCAM"].to_numpy(float), coordinates)
    neighbour_pdl1 = neighbour_mean(section["CD274"].to_numpy(float), coordinates)
    threshold = float(section["CD3E"].quantile(0.75))
    anchors = section["CD3E"].to_numpy(float) >= threshold
    anchors &= np.isfinite(neighbour_epcam)
    if anchors.sum() < 25:
        raise ValueError(f"{gsm} has fewer than 25 eligible T-cell-rich anchors")

    rows = []
    patient, timepoint = SAMPLES[gsm]
    for label, feature in MARKERS.items():
        marker = section[feature].to_numpy(float)[anchors]
        tumour = neighbour_epcam[anchors]
        rho = float(stats.spearmanr(marker, tumour).statistic)
        null = np.empty(N_PERM)
        for i in range(N_PERM):
            null[i] = stats.spearmanr(rng.permutation(marker), tumour).statistic
        empirical_p = float((1 + np.sum(np.abs(null) >= abs(rho))) / (N_PERM + 1))
        pdl1_rho = float(stats.spearmanr(marker, neighbour_pdl1[anchors]).statistic)
        rows.append({"gsm": gsm, "patient": patient, "timepoint": timepoint,
                     "marker": label, "rho_neighbour_epcam": rho,
                     "permutation_p": empirical_p, "rho_neighbour_pdl1": pdl1_rho,
                     "median_anchor_marker": float(np.median(marker)),
                     "n_tissue_spots": len(section), "n_anchors": int(anchors.sum())})
    return rows


def grade_markers(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for marker, group in metrics.groupby("marker", sort=False):
        direction = EXPECTED_DIRECTIONS[marker]
        item = {"marker": marker, "expected_direction": direction}
        for timepoint in ("pre", "post"):
            values = group.loc[group["timepoint"].eq(timepoint), "rho_neighbour_epcam"].to_numpy(float)
            item[f"{timepoint}_median_rho"] = float(np.median(values))
            item[f"{timepoint}_agreement"] = float(np.mean(np.sign(values) == direction))
            item[f"{timepoint}_sign_flip_p"] = exact_sign_flip_p(values)
        same_median = all(np.sign(item[f"{t}_median_rho"]) == direction for t in ("pre", "post"))
        agreement = all(item[f"{t}_agreement"] >= 4 / 6 for t in ("pre", "post"))
        magnitude = max(abs(item["pre_median_rho"]), abs(item["post_median_rho"])) >= 0.10
        if same_median and agreement and magnitude:
            item["grade"] = "PAIRED-SPATIAL-PROTEIN GO"
        elif same_median:
            item["grade"] = "PAIRED-SPATIAL-PROTEIN EXPLORATORY"
        else:
            item["grade"] = "PAIRED-SPATIAL-PROTEIN NO-GO"
        rows.append(item)
    result = pd.DataFrame(rows)
    for timepoint in ("pre", "post"):
        result[f"{timepoint}_sign_flip_q_bh"] = bh_adjust(result[f"{timepoint}_sign_flip_p"])
    return result


def paired_changes(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for marker, group in metrics.groupby("marker", sort=False):
        wide = group.pivot(index="patient", columns="timepoint",
                           values=["rho_neighbour_epcam", "median_anchor_marker"])
        for endpoint in ("rho_neighbour_epcam", "median_anchor_marker"):
            changes = (wide[(endpoint, "post")] - wide[(endpoint, "pre")]).to_numpy(float)
            rows.append({"marker": marker, "endpoint": endpoint,
                         "median_post_minus_pre": float(np.median(changes)),
                         "mean_post_minus_pre": float(np.mean(changes)),
                         "paired_sign_flip_p": exact_sign_flip_p(changes)})
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("data/spatial_protein/GSE289326"))
    parser.add_argument("--output-dir", type=Path,
                        default=Path("results/ai_extension/gse289326_spatial_protein"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    rows = []
    for gsm in SAMPLES:
        rows.extend(section_metrics(load_section(args.data_dir, gsm), gsm, rng))
    metrics = pd.DataFrame(rows)
    summary = grade_markers(metrics)
    changes = paired_changes(metrics)
    metrics.to_csv(args.output_dir / "section_metrics.csv", index=False)
    summary.to_csv(args.output_dir / "replication_summary.csv", index=False)
    changes.to_csv(args.output_dir / "paired_treatment_changes.csv", index=False)
    manifest = {
        "seed": SEED,
        "n_permutations": N_PERM,
        "n_patients": 6,
        "n_unique_sections": 12,
        "technical_replicates_excluded": ["GSM8789211", "GSM8789214"],
        "eligible_markers": list(MARKERS),
        "structurally_unavailable_locked_markers": ["CD103"],
        "grades": dict(zip(summary["marker"], summary["grade"])),
        "claim_boundary": "measured spatial-protein association; not communication or causality",
    }
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(summary.to_string(index=False))
    print(changes.to_string(index=False))


if __name__ == "__main__":
    main()
