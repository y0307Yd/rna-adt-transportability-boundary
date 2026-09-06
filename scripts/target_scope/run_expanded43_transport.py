from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent / "shared43_datasets_v1"
OUT = Path(__file__).resolve().parent / "p3_transport_expanded43_v1"
OUT.mkdir(parents=True, exist_ok=True)
TARGETS = pd.read_csv(Path(__file__).resolve().parent / "P3_SHARED_TARGET_RNA_MAPPING_20260902.csv")


def load(cohort):
    d = ROOT / cohort
    base = Path(r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\p1_external_lopo_validation") / cohort / "v1_reconstruction" / "transport_dataset"
    genes = pd.read_csv(base / "rna_genes.csv").gene.astype(str).tolist()
    gene_map = {g.upper(): i for i, g in enumerate(genes)}
    idx = []
    for value in TARGETS.rna_gene:
        candidates = str(value).split(";")
        found = [gene_map[g.upper()] for g in candidates if g.upper() in gene_map]
        idx.append(found[0] if found else None)
    if any(x is None for x in idx):
        raise ValueError("Missing mapped RNA genes")
    R = sparse.load_npz(base / "rna_raw_cells_by_genes.npz")[:, idx].toarray().astype(np.float32)
    Y = np.load(d / "adt_clr_cells_by_targets.npy").astype(np.float32)
    M = pd.read_csv(base / "cell_metadata.csv")
    return R, Y, M


def corr(a, b):
    return float(spearmanr(a, b).statistic) if np.std(a) > 0 and np.std(b) > 0 else np.nan


def fit_eval(Rtr, Ytr, Rte, Yte, Mte, label, shuffle=False):
    rows = []
    donors = Mte.donor_id.astype(str).to_numpy()
    for j, target in enumerate(TARGETS.gse334503_name):
        y = Ytr[:, j].copy()
        if shuffle:
            y = y[np.random.default_rng(20260902 + j).permutation(len(y))]
        scaler = StandardScaler().fit(Rtr[:, [j]])
        model = Ridge(alpha=1.0).fit(scaler.transform(Rtr[:, [j]]), y)
        pred = model.predict(scaler.transform(Rte[:, [j]]))
        rows.append({"direction": label, "mode": "shuffled" if shuffle else "observed", "target": target,
                     "spearman": corr(pred, Yte[:, j]), "mae": float(np.mean(np.abs(pred - Yte[:, j]))),
                     "n_cells": len(pred), "n_donors": int(pd.Series(donors).nunique()),
                     "median_donor_spearman": float(np.nanmedian([corr(pred[donors == d], Yte[donors == d, j]) for d in np.unique(donors)]))})
    return rows


if __name__ == "__main__":
    np.random.seed(20260902)
    R1, Y1, M1 = load("GSE334503")
    R2, Y2, M2 = load("GSE335494")
    rows = []
    rows += fit_eval(R1, Y1, R1, Y1, M1, "GSE334503_self", False)
    rows += fit_eval(R1, Y1, R1, Y1, M1, "GSE334503_self", True)
    rows += fit_eval(R1, Y1, R2, Y2, M2, "GSE334503_to_GSE335494", False)
    rows += fit_eval(R1, Y1, R2, Y2, M2, "GSE334503_to_GSE335494", True)
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "P3_EXPANDED43_TRANSPORT_METRICS.csv", index=False)
    manifest = {"version": "p3_transport_expanded43_v1", "discovery": "GSE334503", "external_validation": "GSE335494", "ridge_alpha": 1.0,
                "feature_rule": "mapped target cognate RNA only; source-fitted StandardScaler", "seed": 20260902,
                "targets": TARGETS.gse334503_name.tolist(), "status": "complete"}
    (OUT / "P3_EXPANDED43_TRANSPORT_MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(table.groupby(["direction", "mode"]).size())
