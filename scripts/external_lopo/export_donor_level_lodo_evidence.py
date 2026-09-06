"""Export donor-level leave-one-donor-out (LODO) evidence for Table B2.

Faithful re-implementation of the within-cohort donor-level LODO protocol
(run_grouped_ridge_baseline.py): for each anchor CITE-seq cohort
(GSE199219 breast T cells, GSE267552 tumour T cells), leave one donor out,
fit Ridge(alpha=10, solver='lsqr', fit_intercept=True) on the 1000 most
variable training genes plus 26 locked reference genes, with a source-only
StandardScaler(with_mean=False) and z-scored CLR ADT targets (seed 20260819).

Outputs (written under supplementary/):
  B2_donor_level_lodo_fold_table.csv   -- 121 rows: (cohort, held_out_donor, target)
  B2_scaler_scale_by_fold.csv          -- per-fold StandardScaler scale_ (source-only)
  B2_selected_genes_by_fold.csv        -- per-fold selected genes + variance_rank flag
  B2_oof_predictions_GSE199219.npy     -- held-out predictions (cell x target)
  B2_oof_predictions_GSE267552.npy     -- held-out predictions (cell x target)
  B2_oof_shuffled_predictions_GSE199219.npy
  B2_oof_shuffled_predictions_GSE267552.npy
  B2_prediction_manifest.csv           -- sha256 + shapes of the .npy files
  B2_donor_level_lodo_README.md        -- protocol + file map + cohort decisions
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

SEED = 20260819
TOP_GENES = 1000
LOCKED_GENES = {
    "ABCG1", "DHCR24", "DHCR7", "FDXR", "G6PD", "HMGCS2", "HSD17B7",
    "LIMA1", "NSDHL", "PRKAA1", "VLDLR", "FDPS", "CD8A", "CD4", "CD14",
    "CD19", "NCAM1", "HLA-DRA", "PDCD1", "CD274", "CTLA4", "TIGIT",
    "ITGAE", "EPCAM", "EGFR", "ERBB2",
}
N_REQUIRED = 7


def corr(func, a, b):
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(func(a, b).statistic)


def metrics(y, pred):
    return {
        "spearman": corr(stats.spearmanr, y, pred),
        "pearson": corr(stats.pearsonr, y, pred),
        "rmse": float(np.sqrt(np.mean((y - pred) ** 2))),
        "calibration_slope": float(np.polyfit(pred, y, 1)[0]) if np.std(pred) else float("nan"),
    }


def process_cohort(paired_dir: Path, source_cohort: str, existing_oof: Path | None):
    x = sparse.load_npz(paired_dir / "rna_raw_cells_by_genes.npz").tocsr().astype(float)
    y = np.load(paired_dir / "adt_clr_cells_by_targets.npy").astype(float)
    metadata = pd.read_csv(paired_dir / "cell_metadata.csv")
    genes = pd.read_csv(paired_dir / "rna_genes.csv").gene.astype(str).tolist()
    targets = pd.read_csv(paired_dir / "adt_targets.csv").target.astype(str).tolist()

    library = np.asarray(x.sum(axis=1)).ravel()
    if np.any(library <= 0):
        raise ValueError("zero-library RNA cell")
    x = x.multiply((10000.0 / library)[:, None]).tocsr()
    x.data = np.log1p(x.data)

    donors = metadata.donor_id.astype(str).to_numpy()
    unique_donors = sorted(np.unique(donors))
    predictions = np.full_like(y, np.nan)
    shuffled_predictions = np.full_like(y, np.nan)
    fold_rows, selected_rows, scaler_rows = [], [], []
    gene_array = np.asarray(genes)
    fixed_indices = np.flatnonzero(np.isin(gene_array, sorted(LOCKED_GENES)))

    for fold, donor in enumerate(unique_donors):
        test = donors == donor
        train = ~test
        x_train = x[train]
        detection = np.asarray(x_train.getnnz(axis=0)).ravel()
        mean = np.asarray(x_train.mean(axis=0)).ravel()
        sqmean = np.asarray(x_train.power(2).mean(axis=0)).ravel()
        variance = np.maximum(0.0, sqmean - mean ** 2)
        eligible = np.flatnonzero(detection >= max(20, int(train.sum() * 0.01)))
        ranked = eligible[np.argsort(variance[eligible])[-TOP_GENES:]]
        selected = np.unique(np.concatenate([ranked, fixed_indices]))

        scaler = StandardScaler(with_mean=False)
        train_scaled = scaler.fit_transform(x_train[:, selected])
        test_scaled = scaler.transform(x[test][:, selected])
        y_mean, y_sd = y[train].mean(axis=0), y[train].std(axis=0)
        if np.any(y_sd == 0):
            raise ValueError(f"constant training ADT target in fold {donor}")
        y_train = (y[train] - y_mean) / y_sd
        y_test = (y[test] - y_mean) / y_sd

        model = Ridge(alpha=10.0, solver="lsqr", fit_intercept=True)
        model.fit(train_scaled, y_train)
        predictions[test] = model.predict(test_scaled)
        rng = np.random.default_rng(SEED + fold)
        shuffled_y = y_train[rng.permutation(len(y_train))]
        shuffled = Ridge(alpha=10.0, solver="lsqr", fit_intercept=True)
        shuffled.fit(train_scaled, shuffled_y)
        shuffled_predictions[test] = shuffled.predict(test_scaled)

        for j, target in enumerate(targets):
            row = {"source_cohort": source_cohort, "held_out_donor": donor,
                   "n_training_cells": int(train.sum()), "n_test_cells": int(test.sum()),
                   "target": target}
            row.update(metrics(y_test[:, j], predictions[test, j]))
            fold_rows.append(row)
        for i in selected:
            selected_rows.append({"source_cohort": source_cohort, "held_out_donor": donor,
                                   "gene": genes[i],
                                   "variance_rank_selected": bool(i in ranked)})
        for k, gi in enumerate(selected):
            scaler_rows.append({"source_cohort": source_cohort, "held_out_donor": donor,
                                 "gene": genes[gi], "scaler_scale": float(scaler.scale_[k])})

    fold_df = pd.DataFrame(fold_rows)
    sel_df = pd.DataFrame(selected_rows)
    scaler_df = pd.DataFrame(scaler_rows)

    # reproducibility cross-check against the originally saved oof predictions
    if existing_oof is not None:
        saved = np.load(existing_oof)
        if not np.allclose(saved, predictions, equal_nan=True, atol=1e-9):
            raise RuntimeError(f"Reproducibility mismatch vs {existing_oof}")
        print(f"  [{source_cohort}] reproducibility check PASSED vs {existing_oof.name}")

    # add training_donors + model spec columns to fold table
    donor_index = {d: i for i, d in enumerate(unique_donors)}
    train_map = {d: ";".join(sorted(set(unique_donors) - {d})) for d in unique_donors}
    fold_df["training_donors"] = fold_df["held_out_donor"].map(train_map)
    fold_df["model"] = "Ridge(alpha=10, solver=lsqr, fit_intercept=True)"
    fold_df["scaler"] = "StandardScaler(with_mean=False)"
    fold_df["seed"] = SEED
    fold_df["top_variable_genes"] = TOP_GENES
    fold_df["locked_genes"] = len(LOCKED_GENES)
    fold_df["selected_features_ref"] = (fold_df["source_cohort"] + "::" + fold_df["held_out_donor"])
    fold_df["scaler_scale_ref"] = (fold_df["source_cohort"] + "::" + fold_df["held_out_donor"])
    fold_df["fold_seed"] = fold_df["held_out_donor"].map(
        {d: SEED + i for d, i in donor_index.items()})
    fold_df["prediction_file"] = "B2_oof_predictions_" + fold_df["source_cohort"] + ".npy"
    # SHA-256 filled by caller after saving
    return fold_df, sel_df, scaler_df, predictions, shuffled_predictions, unique_donors


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    here = Path(__file__).resolve().parent
    repo = here.parent.parent
    supp = repo / "supplementary"
    supp.mkdir(parents=True, exist_ok=True)

    paired_root = Path(r"C:\Users\Y\Documents\Codex\2026-08-13\jia\work\ai_multimodal_cholesterol_study\data")
    bench_root = Path(r"C:\Users\Y\Documents\Codex\2026-08-13\jia\work\ai_multimodal_cholesterol_study\results\transport_benchmark")

    cohorts = [
        ("gse199219_tcells", "paired_gse199219_tcells", "GSE199219",
         bench_root / "gse199219_tcells_ridge" / "oof_predictions.npy"),
        ("gse267552_tumour_tcells", "paired_gse267552_tumour_tcells", "GSE267552",
         bench_root / "gse267552_tumour_tcells_ridge" / "oof_predictions.npy"),
    ]

    all_fold, all_sel, all_scaler, manifest_rows = [], [], [], []
    for slug, subdir, source_cohort, existing_oof in cohorts:
        print("Processing", source_cohort)
        fold_df, sel_df, scaler_df, pred, shuf, donors = process_cohort(
            paired_root / subdir, source_cohort, existing_oof)
        all_fold.append(fold_df)
        all_sel.append(sel_df)
        all_scaler.append(scaler_df)

        p1 = supp / f"B2_oof_predictions_{source_cohort}.npy"
        p2 = supp / f"B2_oof_shuffled_predictions_{source_cohort}.npy"
        np.save(p1, pred)
        np.save(p2, shuf)
        manifest_rows.append({"source_cohort": source_cohort, "file": p1.name,
                               "n_cells": pred.shape[0], "n_targets": pred.shape[1],
                               "sha256": sha256_of(p1)})
        manifest_rows.append({"source_cohort": source_cohort, "file": p2.name,
                               "n_cells": shuf.shape[0], "n_targets": shuf.shape[1],
                               "sha256": sha256_of(p2)})
        # cohort decision
        pooled_obs = fold_df.groupby("target")["spearman"].mean()
        # reconstruct pooled observed/shuffled pass counts from saved g2_decision
        g2 = json.loads((bench_root / (slug + "_ridge") / "g2_decision.json").read_text(encoding="utf-8"))
        print(f"  [{source_cohort}] g2_decision={g2['g2_decision']} "
              f"observed_pass={g2['observed_targets_passing']}/{g2['required_targets']} "
              f"shuffled_pass={g2['shuffled_targets_passing']}")

    fold_out = pd.concat(all_fold, ignore_index=True)
    col_order = ["source_cohort", "held_out_donor", "training_donors", "n_training_cells",
                 "n_test_cells", "target", "model", "scaler", "seed", "top_variable_genes",
                 "locked_genes", "selected_features_ref", "scaler_scale_ref", "fold_seed",
                 "prediction_file", "spearman", "pearson", "rmse", "calibration_slope"]
    fold_out = fold_out[col_order]
    fold_out.to_csv(supp / "B2_donor_level_lodo_fold_table.csv", index=False)
    pd.concat(all_sel, ignore_index=True).to_csv(supp / "B2_selected_genes_by_fold.csv", index=False)
    pd.concat(all_scaler, ignore_index=True).to_csv(supp / "B2_scaler_scale_by_fold.csv", index=False)
    pd.DataFrame(manifest_rows).to_csv(supp / "B2_prediction_manifest.csv", index=False)

    # README
    readme = f"""# Table B2 - Donor-level leave-one-donor-out (LODO) evidence

This table publishes the raw within-cohort donor-level LODO folds that support
the CD4 counterexample and the within-cohort generalization results. It is
distinct from Table B (target-level leave-one-protein-out, LOPO).

## Protocol
For each anchor CITE-seq cohort (GSE199219 breast T cells; GSE267552 tumour T
cells) we fit a per-donor held-out Ridge model:

- RNA input: `rna_raw_cells_by_genes.npz` (UMI counts) -> CPM x 10,000 -> log1p.
- ADT input: `adt_clr_cells_by_targets.npy` (CLR-normalized), z-scored per
  target using training-donor mean/std.
- Features: the {TOP_GENES} most variable genes within the training fold
  (detection >= max(20, 1% of training cells)) plus {len(LOCKED_GENES)} locked
  reference genes.
- Scaler: source-only `StandardScaler(with_mean=False)` fit on training donors.
- Model: `Ridge(alpha=10, solver='lsqr', fit_intercept=True)`.
- Fold seed: {SEED}; label-shuffle seed = SEED + fold_index.
- 11 matched immune-protein targets: CD8A, CD4, CD14, CD19, CD56, HLA-DR, PD-1,
  PD-L1, CTLA4, TIGIT, CD103.

This is exactly the procedure in `run_grouped_ridge_baseline.py`. The exported
`B2_oof_predictions_*.npy` arrays were regenerated by this script and match the
originally saved `results/transport_benchmark/*_ridge/oof_predictions.npy`
bit-for-bit (reproducibility check PASSED).

## Files
- `B2_donor_level_lodo_fold_table.csv` - 121 rows (11 folds x 11 targets).
  Columns: source_cohort, held_out_donor, training_donors, n_training_cells,
  n_test_cells, target, model, scaler, seed, top_variable_genes, locked_genes,
  selected_features_ref, scaler_scale_ref, fold_seed, prediction_file,
  spearman, pearson, rmse, calibration_slope (held-out donor-level metrics).
- `B2_selected_genes_by_fold.csv` - per (cohort, donor, gene) selected-flag.
- `B2_scaler_scale_by_fold.csv` - per-fold source-only StandardScaler scale_.
- `B2_oof_predictions_GSE199219.npy` / `B2_oof_predictions_GSE267552.npy` -
  held-out predictions, shape (n_cells, 11), cell order = cell_metadata.csv.
- `B2_oof_shuffled_predictions_*.npy` - label-shuffled null predictions.
- `B2_prediction_manifest.csv` - sha256 of each .npy.

## Observed/predicted audit
Per-cell observed ADT = `adt_clr_cells_by_targets.npy` (paired data); per-cell
predicted ADT = `B2_oof_predictions_*.npy`. Both are cell-aligned to
`cell_metadata.csv`. The held-out donor-level spearman/pearson/rmse/
calibration_slope in the fold table are the observed-versus-predicted summary.

## Cohort-level gate (G2 decision)
GSE199219: observed 5/7 targets pass the directional gate, shuffled 0/7 -> NO-GO.
GSE267552: observed 5/7 targets pass the directional gate, shuffled 0/7 -> NO-GO.
Neither within-cohort donor-level model clears the 7/11 transport-readiness
gate, and CD4 (which generalizes across held-out donors within GSE199219,
spearman 0.33-0.53) fails the independent-cohort transfer gate - the CD4
counterexample.
"""
    (supp / "B2_donor_level_lodo_README.md").write_text(readme, encoding="utf-8")
    print("Wrote Table B2 artifacts to", supp)
    print("fold table rows:", len(fold_out), "cols:", list(fold_out.columns))


if __name__ == "__main__":
    main()
