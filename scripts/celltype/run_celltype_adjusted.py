"""Cell-type adjusted transportability + cell-type baseline + worst-donor analysis.

Addresses reviewer concerns:
  MC5  pooled Spearman may be driven by cell composition -> within-cell-type + adjusted.
  R30  cell-type/pseudobulk baseline -> does model learn RNA->protein or just cell identity?
  R18  worst-donor stability -> median / 10th pct / min donor Spearman.

Design (per cohort, cognate RNA only, matches locked rule alpha=1 source-only):
  - pooled   : Spearman(log1p(cognateRNA), CLR(ADT)) over all cells
  - within   : Spearman within each cell_type (n>=100), report median & min across cell_types
  - adjusted : Spearman after removing cell_type mean (cell-type-adjusted residual)
  - baseline : cell-type one-hot mean baseline Spearman (predict ADT = cell-type mean ADT)
  - donor    : median / p10 / min of per-donor Spearman
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr

BASE = Path(r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap")
OUT = BASE / "celltype_annotation_20260905"
ANN = OUT  # annotated metadata lives here
MAPPING = pd.read_csv(BASE / "p3_target_scope_20260902/P3_SHARED_TARGET_RNA_MAPPING_20260902.csv")


def load_cohort(cohort):
    td = BASE / f"p1_external_lopo_validation/{cohort}/v1_reconstruction/transport_dataset"
    genes = pd.read_csv(td / "rna_genes.csv").gene.astype(str).tolist()
    gmap = {g.upper(): i for i, g in enumerate(genes)}

    # resolve cognate RNA index per target (first candidate gene)
    rna_idx = []
    targets = MAPPING["gse334503_name" if cohort == "GSE334503" else "gse335494_name"].tolist()
    for i, gene_str in enumerate(MAPPING["rna_gene"]):
        cands = [g.strip().upper() for g in str(gene_str).split(";")]
        found = None
        for c in cands:
            if c in gmap:
                found = gmap[c]
                break
        rna_idx.append(found)

    R = sparse.load_npz(td / "rna_raw_cells_by_genes.npz")[:, rna_idx].toarray().astype(np.float32)
    Y = np.load(BASE / f"p3_target_scope_20260902/shared43_datasets_v1/{cohort}/adt_clr_cells_by_targets.npy").astype(np.float32)
    meta = pd.read_csv(ANN / f"{cohort}_cell_metadata_annotated.csv")
    return R, Y, meta, targets


def spearman(a, b):
    if np.std(a) == 0 or np.std(b) == 0 or len(a) < 3:
        return np.nan
    return float(spearmanr(a, b).statistic)


def analyze(cohort):
    R, Y, meta, targets = load_cohort(cohort)
    X = np.log1p(R)  # cognate RNA log1p
    donors = meta["donor_id"].astype(str).to_numpy()
    ctypes = meta["cell_type"].to_numpy()
    uniq_donors = np.unique(donors)
    uniq_ct = [c for c in np.unique(ctypes) if c != "Ambiguous"]

    # precompute cell-type one-hot means for baseline
    ct_codes = pd.Categorical(ctypes, categories=uniq_ct)
    rows = []
    for j, t in enumerate(targets):
        x = X[:, j]
        y = Y[:, j]
        # exclude Ambiguous (low-confidence) cells for cell-type-aware analyses
        valid = np.isfinite(y) & np.isfinite(x) & (ctypes != "Ambiguous")
        xv, yv, dv, cv = x[valid], y[valid], donors[valid], ct_codes[valid]

        pooled = spearman(xv, yv)

        # within-cell-type
        within = []
        for c in uniq_ct:
            m = cv == c
            if m.sum() >= 100:
                s = spearman(xv[m], yv[m])
                if np.isfinite(s):
                    within.append(s)
        within_median = float(np.median(within)) if within else np.nan
        within_min = float(np.min(within)) if within else np.nan

        # cell-type adjusted residual
        ct_mean_y = {c: yv[cv == c].mean() for c in uniq_ct}
        y_adj = yv - np.array([ct_mean_y[c] for c in cv])
        adjusted = spearman(xv, y_adj)

        # cell-type baseline: predict ADT = cell-type mean ADT
        pred_baseline = np.array([ct_mean_y[c] for c in cv])
        baseline = spearman(pred_baseline, yv)

        # donor-level
        dspear = []
        for d in uniq_donors:
            m = dv == d
            if m.sum() >= 50:
                s = spearman(xv[m], yv[m])
                if np.isfinite(s):
                    dspear.append(s)
        dspear = np.array(dspear)
        donor_median = float(np.median(dspear)) if len(dspear) else np.nan
        donor_p10 = float(np.percentile(dspear, 10)) if len(dspear) else np.nan
        donor_min = float(np.min(dspear)) if len(dspear) else np.nan
        donor_n = int(len(dspear))

        rows.append({
            "cohort": cohort, "target": t,
            "pooled": pooled,
            "within_ct_median": within_median, "within_ct_min": within_min,
            "ct_adjusted": adjusted,
            "ct_baseline": baseline,
            "donor_median": donor_median, "donor_p10": donor_p10, "donor_min": donor_min,
            "donor_n": donor_n,
            "n_cells": int(len(xv)),
        })

    out = pd.DataFrame(rows)
    out.to_csv(OUT / f"{cohort}_celltype_adjusted_transport.csv", index=False)
    print(f"=== {cohort} ({len(targets)} targets) ===")
    print(out[["target", "pooled", "within_ct_median", "ct_adjusted", "ct_baseline", "donor_median", "donor_p10", "donor_min"]].round(3).to_string(index=False))
    print()
    return out


if __name__ == "__main__":
    for c in ["GSE334503", "GSE335494"]:
        analyze(c)
