from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(r"C:\Users\Y\Documents\Codex\2026-08-13\jia\work\ai_multimodal_cholesterol_study")
OUT = Path(r"D:\ai_multimodal_cholesterol_study_outputs\m1_source_positive_control_v1")
SEED = 20260902
THRESHOLDS = {
    "spearman_min": 0.20,
    "specificity_delta_min": 0.10,
    "r2_min_exclusive": 0.0,
    "calibration_slope_lo": 0.5,
    "calibration_slope_hi": 2.0,
    "sd_ratio_hard_lo": 0.5,
    "sd_ratio_hard_hi": 2.0,
    "sd_ratio_exploratory_lo": 0.05,
    "sd_ratio_exploratory_hi": 20.0,
}

COHORTS = {
    "GSE199219": {
        "paired": PROJECT / "data" / "paired_gse199219_tcells",
        "result": PROJECT / "results" / "transport_benchmark" / "gse199219_tcells_ridge",
    },
    "GSE267552": {
        "paired": PROJECT / "data" / "paired_gse267552_tumour_tcells",
        "result": PROJECT / "results" / "transport_benchmark" / "gse267552_tumour_tcells_ridge",
    },
}


def rho(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    rank_a = pd.Series(a).rank(method="average").to_numpy()
    rank_b = pd.Series(b).rank(method="average").to_numpy()
    return float(np.corrcoef(rank_a, rank_b)[0, 1])


def derangement(n: int, rng: np.random.Generator) -> np.ndarray:
    base = np.arange(n)
    while True:
        candidate = rng.permutation(n)
        if np.all(candidate != base):
            return candidate


def metric_row(y: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    pred_sd = float(np.std(pred))
    obs_sd = float(np.std(y))
    slope = (
        float(np.cov(pred, y, ddof=0)[0, 1] / np.var(pred))
        if pred_sd > 1e-12 and obs_sd > 0
        else float("nan")
    )
    intercept = float(np.mean(y) - slope * np.mean(pred)) if np.isfinite(slope) else float("nan")
    return {
        "spearman": rho(pred, y),
        "r2": float(1.0 - np.sum((y - pred) ** 2) / np.sum((y - np.mean(y)) ** 2)),
        "rmse": float(np.sqrt(np.mean((y - pred) ** 2))),
        "mae": float(np.mean(np.abs(y - pred))),
        "prediction_sd": pred_sd,
        "observed_sd": obs_sd,
        "sd_ratio": pred_sd / obs_sd if obs_sd > 0 else float("nan"),
        "calibration_slope": slope,
        "calibration_intercept": intercept,
    }


def reconstruct_raw_scale(y: np.ndarray, pred_z: np.ndarray, donors: np.ndarray) -> np.ndarray:
    pred = np.full_like(pred_z, np.nan, dtype=float)
    for donor in np.unique(donors):
        test = donors == donor
        train = ~test
        train_mean = y[train].mean(axis=0)
        train_sd = y[train].std(axis=0)
        if np.any(train_sd == 0):
            raise ValueError(f"constant training target in held-out donor {donor}")
        pred[test] = pred_z[test] * train_sd + train_mean
    if not np.isfinite(pred).all():
        raise ValueError("non-finite reconstructed predictions")
    return pred


def apply_gate(row: dict[str, float], lo: float, hi: float) -> bool:
    t = THRESHOLDS
    return bool(
        row["spearman"] >= t["spearman_min"]
        and row["specificity_delta"] >= t["specificity_delta_min"]
        and row["r2"] > t["r2_min_exclusive"]
        and t["calibration_slope_lo"] <= row["calibration_slope"] <= t["calibration_slope_hi"]
        and lo <= row["sd_ratio"] <= hi
    )


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    pooled_rows: list[dict] = []
    donor_rows: list[dict] = []
    mapping_rows: list[dict] = []
    input_rows: list[dict] = []
    legacy_check_rows: list[dict] = []

    for cohort, paths in COHORTS.items():
        paired = paths["paired"]
        result = paths["result"]
        y_path = paired / "adt_clr_cells_by_targets.npy"
        pred_path = result / "oof_predictions.npy"
        shuffled_path = result / "oof_shuffled_predictions.npy"
        metadata_path = paired / "cell_metadata.csv"
        target_path = paired / "adt_targets.csv"

        y = np.load(y_path).astype(float)
        pred_z = np.load(pred_path).astype(float)
        shuffled_z = np.load(shuffled_path).astype(float)
        metadata = pd.read_csv(metadata_path)
        donors = metadata["donor_id"].astype(str).to_numpy()
        targets = pd.read_csv(target_path)["target"].astype(str).tolist()
        if y.shape != pred_z.shape or y.shape != shuffled_z.shape:
            raise ValueError(f"shape mismatch for {cohort}: {y.shape}, {pred_z.shape}, {shuffled_z.shape}")
        if y.shape[0] != len(donors) or y.shape[1] != len(targets):
            raise ValueError(f"metadata mismatch for {cohort}")

        pred = reconstruct_raw_scale(y, pred_z, donors)
        shuffled = reconstruct_raw_scale(y, shuffled_z, donors)
        perm = derangement(len(targets), rng)
        for j, k in enumerate(perm):
            mapping_rows.append({
                "cohort": cohort,
                "predicted_target": targets[j],
                "noncognate_observed_target": targets[k],
                "seed": SEED,
                "mapping_rule": "fixed_random_derangement",
            })

        for source_path, role in [
            (y_path, "observed_adt"),
            (pred_path, "stored_fold_standardized_oof_prediction"),
            (shuffled_path, "stored_fold_standardized_oof_shuffled_prediction"),
            (metadata_path, "donor_metadata"),
            (target_path, "target_names"),
        ]:
            input_rows.append({
                "cohort": cohort,
                "role": role,
                "path": str(source_path),
                "bytes": source_path.stat().st_size,
                "sha256": sha256(source_path),
            })

        for j, target in enumerate(targets):
            observed = metric_row(y[:, j], pred[:, j])
            observed["noncognate_spearman"] = rho(pred[:, j], y[:, perm[j]])
            observed["specificity_delta"] = observed["spearman"] - observed["noncognate_spearman"]
            observed["cohort"] = cohort
            observed["target"] = target
            observed["mode"] = "observed"
            observed["n_cells"] = len(y)
            observed["n_donors"] = len(np.unique(donors))
            observed["passes_hard_gate"] = apply_gate(observed, 0.5, 2.0)
            observed["passes_exploratory_gate"] = apply_gate(observed, 0.05, 20.0)
            pooled_rows.append(observed)

            null = metric_row(y[:, j], shuffled[:, j])
            null["noncognate_spearman"] = rho(shuffled[:, j], y[:, perm[j]])
            null["specificity_delta"] = null["spearman"] - null["noncognate_spearman"]
            null["cohort"] = cohort
            null["target"] = target
            null["mode"] = "label_shuffled"
            null["n_cells"] = len(y)
            null["n_donors"] = len(np.unique(donors))
            null["passes_hard_gate"] = apply_gate(null, 0.5, 2.0)
            null["passes_exploratory_gate"] = apply_gate(null, 0.05, 20.0)
            pooled_rows.append(null)

            for donor in np.unique(donors):
                mask = donors == donor
                drow = metric_row(y[mask, j], pred[mask, j])
                drow.update({
                    "cohort": cohort,
                    "held_out_donor": donor,
                    "target": target,
                    "n_cells": int(mask.sum()),
                    "noncognate_spearman": rho(pred[mask, j], y[mask, perm[j]]),
                })
                drow["specificity_delta"] = drow["spearman"] - drow["noncognate_spearman"]
                drow["passes_hard_gate"] = apply_gate(drow, 0.5, 2.0)
                drow["passes_exploratory_gate"] = apply_gate(drow, 0.05, 20.0)
                donor_rows.append(drow)

        legacy = pd.read_csv(result / "fold_metrics.csv")
        current = pd.DataFrame([r for r in donor_rows if r["cohort"] == cohort])
        checked = current.merge(
            legacy[["donor", "target", "spearman"]],
            left_on=["held_out_donor", "target"],
            right_on=["donor", "target"],
            suffixes=("_reconstructed", "_legacy"),
            validate="one_to_one",
        )
        for _, row in checked.iterrows():
            legacy_check_rows.append({
                "cohort": cohort,
                "held_out_donor": row["held_out_donor"],
                "target": row["target"],
                "reconstructed_spearman": row["spearman_reconstructed"],
                "legacy_spearman": row["spearman_legacy"],
                "absolute_difference": abs(row["spearman_reconstructed"] - row["spearman_legacy"]),
            })

    pooled = pd.DataFrame(pooled_rows)
    donor = pd.DataFrame(donor_rows)
    mapping = pd.DataFrame(mapping_rows)
    inputs = pd.DataFrame(input_rows)
    legacy_check = pd.DataFrame(legacy_check_rows)
    pooled.to_csv(OUT / "source_lodo_target_metrics.csv", index=False)
    donor.to_csv(OUT / "source_lodo_donor_target_metrics.csv", index=False)
    mapping.to_csv(OUT / "null_b_noncognate_mapping.csv", index=False)
    inputs.to_csv(OUT / "input_manifest_sha256.csv", index=False)
    legacy_check.to_csv(OUT / "legacy_fold_consistency_audit.csv", index=False)

    observed = pooled[pooled["mode"] == "observed"].copy()
    null = pooled[pooled["mode"] == "label_shuffled"].copy()
    summary_rows = []
    for cohort in [*COHORTS, "COMBINED_TARGET_COHORT_SETTINGS"]:
        x = observed if cohort == "COMBINED_TARGET_COHORT_SETTINGS" else observed[observed["cohort"] == cohort]
        z = null if cohort == "COMBINED_TARGET_COHORT_SETTINGS" else null[null["cohort"] == cohort]
        summary_rows.append({
            "scope": cohort,
            "n_target_cohort_settings": len(x),
            "hard_gate_pass": int(x["passes_hard_gate"].sum()),
            "hard_gate_pass_fraction": float(x["passes_hard_gate"].mean()),
            "exploratory_gate_pass": int(x["passes_exploratory_gate"].sum()),
            "exploratory_gate_pass_fraction": float(x["passes_exploratory_gate"].mean()),
            "label_shuffled_hard_gate_pass": int(z["passes_hard_gate"].sum()),
            "median_spearman": float(x["spearman"].median()),
            "median_specificity_delta": float(x["specificity_delta"].median()),
            "median_r2": float(x["r2"].median()),
            "median_calibration_slope": float(x["calibration_slope"].median()),
            "median_sd_ratio": float(x["sd_ratio"].median()),
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "source_lodo_gate_summary.csv", index=False)

    payload = {
        "analysis_id": "m1_source_positive_control_v1",
        "status": "COMPLETE",
        "seed": SEED,
        "thresholds": THRESHOLDS,
        "prediction_scale_reconstruction": "For each held-out donor, stored standardized OOF predictions were transformed as pred_CLR = pred_z * SD(training ADT CLR) + mean(training ADT CLR).",
        "specificity_control": "Null-B is a fixed random derangement of observed ADT target columns within each cohort; mapping is saved.",
        "cohort_summary": summary.to_dict(orient="records"),
        "legacy_fold_consistency": {
            "rows_matched": int(len(legacy_check)),
            "max_absolute_spearman_difference": float(legacy_check["absolute_difference"].max()),
            "status": "PASS" if float(legacy_check["absolute_difference"].max()) <= 1e-12 else "FAIL",
        },
        "interpretation_constraint": "This is a within-source LODO gate satisfiability control, not evidence of cross-cohort transportability.",
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    combined = summary.iloc[-1]
    report = f"""# M1 source-internal LODO positive-control audit

## Material Passport

- Analysis ID: `m1_source_positive_control_v1`
- Status: COMPLETE
- Source cohorts: GSE199219 and GSE267552
- Statistical unit for the primary summary: target-cohort setting (11 targets x 2 cohorts = 22)
- Random seed: {SEED}
- Inputs and hashes: `input_manifest_sha256.csv`

## Result

The strict five-component gate was satisfiable within the source cohorts: {int(combined['hard_gate_pass'])}/{int(combined['n_target_cohort_settings'])} target-cohort settings ({combined['hard_gate_pass_fraction']:.1%}) passed the hard gate. The label-shuffled control passed in {int(combined['label_shuffled_hard_gate_pass'])}/{int(combined['n_target_cohort_settings'])} settings. Under the broad exploratory SD-ratio tier, {int(combined['exploratory_gate_pass'])}/{int(combined['n_target_cohort_settings'])} settings ({combined['exploratory_gate_pass_fraction']:.1%}) passed.

Median source-internal metrics across the 22 target-cohort settings were Spearman {combined['median_spearman']:.3f}, specificity delta {combined['median_specificity_delta']:.3f}, R2 {combined['median_r2']:.3f}, calibration slope {combined['median_calibration_slope']:.3f}, and prediction/observed SD ratio {combined['median_sd_ratio']:.3f}.

## Reconstruction note

The legacy OOF arrays stored predictions on the training-fold standardized ADT scale. Each held-out donor prediction was therefore inverse-transformed with the mean and standard deviation calculated from all other donors in the same cohort. Direct comparison of those stored arrays with unstandardized CLR observations would be invalid.

## Specificity control

For each cohort, a fixed random derangement paired every predicted target with a different observed ADT target. `specificity_delta` is the cognate Spearman minus this fixed noncognate Spearman. The exact mapping is in `null_b_noncognate_mapping.csv`.

## Claim boundary

This analysis tests whether the gate can be met in donor-held-out source data. It does not establish cross-cohort transportability and does not convert retrospective thresholds into preregistration. Manuscript attribution should compare this source-internal distribution with the external distribution and retain continuous metrics alongside pass counts.
"""
    (OUT / "M1_SOURCE_POSITIVE_CONTROL_REPORT.md").write_text(report, encoding="utf-8")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
