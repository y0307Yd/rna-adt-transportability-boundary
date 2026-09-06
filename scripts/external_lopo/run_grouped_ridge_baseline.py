from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler


SEED = 20260819
LOCKED_GENES = {
    "ABCG1", "DHCR24", "DHCR7", "FDXR", "G6PD", "HMGCS2", "HSD17B7",
    "LIMA1", "NSDHL", "PRKAA1", "VLDLR", "FDPS", "CD8A", "CD4", "CD14",
    "CD19", "NCAM1", "HLA-DRA", "PDCD1", "CD274", "CTLA4", "TIGIT",
    "ITGAE", "EPCAM", "EGFR", "ERBB2",
}


def corr(func, a: np.ndarray, b: np.ndarray) -> float:
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(func(a, b).statistic)


def metrics(y: np.ndarray, pred: np.ndarray) -> dict:
    return {
        "spearman": corr(stats.spearmanr, y, pred),
        "pearson": corr(stats.pearsonr, y, pred),
        "rmse": float(np.sqrt(np.mean((y - pred) ** 2))),
        "calibration_slope": float(np.polyfit(pred, y, 1)[0]) if np.std(pred) else float("nan"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paired_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--top-genes", type=int, default=1000)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    x = sparse.load_npz(args.paired_dir / "rna_raw_cells_by_genes.npz").tocsr().astype(float)
    y = np.load(args.paired_dir / "adt_clr_cells_by_targets.npy").astype(float)
    metadata = pd.read_csv(args.paired_dir / "cell_metadata.csv")
    genes = pd.read_csv(args.paired_dir / "rna_genes.csv").gene.astype(str).tolist()
    targets = pd.read_csv(args.paired_dir / "adt_targets.csv").target.astype(str).tolist()

    library = np.asarray(x.sum(axis=1)).ravel()
    if np.any(library <= 0):
        raise ValueError("zero-library RNA cell")
    x = x.multiply((10000.0 / library)[:, None]).tocsr()
    x.data = np.log1p(x.data)

    donors = metadata.donor_id.astype(str).to_numpy()
    unique_donors = sorted(np.unique(donors))
    predictions = np.full_like(y, np.nan)
    shuffled_predictions = np.full_like(y, np.nan)
    fold_rows, selected_rows = [], []
    gene_array = np.asarray(genes)
    fixed_indices = np.flatnonzero(np.isin(gene_array, sorted(LOCKED_GENES)))

    for fold, donor in enumerate(unique_donors):
        test = donors == donor
        train = ~test
        x_train = x[train]
        detection = np.asarray(x_train.getnnz(axis=0)).ravel()
        mean = np.asarray(x_train.mean(axis=0)).ravel()
        sqmean = np.asarray(x_train.power(2).mean(axis=0)).ravel()
        variance = np.maximum(0.0, sqmean - mean**2)
        eligible = np.flatnonzero(detection >= max(20, int(train.sum() * 0.01)))
        ranked = eligible[np.argsort(variance[eligible])[-args.top_genes:]]
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
            row = {"donor": donor, "target": target, "n_test": int(test.sum())}
            row.update(metrics(y_test[:, j], predictions[test, j]))
            fold_rows.append(row)
        selected_rows.extend(
            {"donor": donor, "gene": genes[i], "variance_rank_selected": bool(i in ranked)}
            for i in selected
        )

    pooled_rows = []
    for mode, pred in [("observed", predictions), ("shuffled", shuffled_predictions)]:
        for j, target in enumerate(targets):
            row = {"mode": mode, "target": target, "n": len(y)}
            row.update(metrics(y[:, j], pred[:, j]))
            within = [corr(stats.spearmanr, y[donors == donor, j], pred[donors == donor, j])
                      for donor in unique_donors]
            row["median_within_donor_spearman"] = float(np.nanmedian(within))
            row["passes_directional_gate"] = bool(
                row["spearman"] > 0.10 and row["median_within_donor_spearman"] > 0.05
            )
            pooled_rows.append(row)

    pooled = pd.DataFrame(pooled_rows)
    observed_pass = int(pooled.query("mode == 'observed'").passes_directional_gate.sum())
    shuffled_pass = int(pooled.query("mode == 'shuffled'").passes_directional_gate.sum())
    decision = "GO" if observed_pass >= 7 and shuffled_pass < 7 else "NO-GO"

    np.save(args.output_dir / "oof_predictions.npy", predictions)
    np.save(args.output_dir / "oof_shuffled_predictions.npy", shuffled_predictions)
    pd.DataFrame(fold_rows).to_csv(args.output_dir / "fold_metrics.csv", index=False)
    pooled.to_csv(args.output_dir / "pooled_metrics.csv", index=False)
    pd.DataFrame(selected_rows).to_csv(args.output_dir / "selected_genes_by_fold.csv", index=False)
    payload = {
        "seed": SEED, "folds": unique_donors,
        "model": "Ridge(alpha=10, solver=lsqr)",
        "top_variable_genes_per_training_fold": args.top_genes,
        "observed_targets_passing": observed_pass,
        "shuffled_targets_passing": shuffled_pass,
        "required_targets": 7, "g2_decision": decision,
    }
    (args.output_dir / "g2_decision.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
