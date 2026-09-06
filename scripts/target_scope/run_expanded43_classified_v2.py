"""Apply canonical transport classification (TRANSPORTABILITY_RULE_FORMALIZATION §2)
to the P3 43-target main result (GSE334503 -> GSE335494 cross-cohort transport).

P3 v1 only reported pooled Spearman + MAE + median_donor_spearman and lacked:
  - Null-B noncognate specificity_delta
  - R2, calibration slope, prediction/observed SD ratio
  - worst-donor (donor_min / donor_p10) on the VALIDATION cohort

This script computes all of the above and emits a classified main-results table,
making the 43-target conclusion reproducible under the canonical rule.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score

BASE = Path(r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap")
P3 = BASE / "p3_target_scope_20260902"
OUT = P3 / "p3_transport_expanded43_v2_classified"
OUT.mkdir(parents=True, exist_ok=True)
MAPPING = pd.read_csv(P3 / "P3_SHARED_TARGET_RNA_MAPPING_20260902.csv")
SRC43 = P3 / "shared43_datasets_v1"
SRCBASE = BASE / "p1_external_lopo_validation"


def load(cohort):
    td = SRCBASE / cohort / "v1_reconstruction" / "transport_dataset"
    genes = pd.read_csv(td / "rna_genes.csv").gene.astype(str).tolist()
    lut = {g.upper(): i for i, g in enumerate(genes)}
    idx = []
    for s in MAPPING.rna_gene.astype(str):
        h = [lut[x.upper()] for x in s.split(";") if x.upper() in lut]
        idx.append(h[0] if h else None)
    if any(x is None for x in idx):
        raise ValueError(f"{cohort} missing RNA gene")
    R = sparse.load_npz(td / "rna_raw_cells_by_genes.npz")[:, idx].toarray().astype(np.float32)
    Y = np.load(SRC43 / cohort / "adt_clr_cells_by_targets.npy").astype(np.float32)
    M = pd.read_csv(td / "cell_metadata.csv")
    return R, Y, M


def spearman(a, b):
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float(spearmanr(a, b).statistic)


def classify(rho, delta, r2, slope, sd_ratio):
    rank = (rho >= 0.20) and (delta >= 0.10)
    numeric = rank and (r2 > 0) and (0.5 <= slope <= 2.0) and (0.05 <= sd_ratio <= 20.0)
    if numeric:
        return "numerically_transportable"
    if rank:
        return "rank_invariant_but_uncalibrated"
    if rho >= 0.10:
        return "weak_context_dependent"
    return "non_specific_or_non_transportable"


def main():
    R1, Y1, M1 = load("GSE334503")   # discovery
    R2, Y2, M2 = load("GSE335494")   # validation
    targets = MAPPING.gse334503_name.tolist()

    rng = np.random.default_rng(20260902)
    perm = rng.permutation(43)  # fixed noncognate pairing (Null-B)

    donors = M2.donor_id.astype(str).to_numpy()
    uniq_donors = np.unique(donors)

    rows = []
    for j, target in enumerate(targets):
        sc = StandardScaler().fit(R1[:, [j]])
        model = Ridge(alpha=1.0).fit(sc.transform(R1[:, [j]]), Y1[:, j])
        pred = model.predict(sc.transform(R2[:, [j]]))
        y_obs = Y2[:, j]
        y_nc = Y2[:, perm[j]]

        rho_obs = spearman(pred, y_obs)
        rho_nc = spearman(pred, y_nc)
        delta = rho_obs - rho_nc if np.isfinite(rho_obs) and np.isfinite(rho_nc) else np.nan

        r2 = float(r2_score(y_obs, pred))
        ps, ys = float(np.std(pred)), float(np.std(y_obs))
        sd_ratio = ps / ys if ys > 0 else np.nan
        if ps > 1e-8 and ys > 0:
            slope = float(np.cov(pred, y_obs, ddof=0)[0, 1] / np.var(pred))
        else:
            slope = np.nan

        # worst-donor on validation cohort
        dspear = []
        for d in uniq_donors:
            m = donors == d
            if m.sum() >= 50:
                s = spearman(pred[m], y_obs[m])
                if np.isfinite(s):
                    dspear.append(s)
        dspear = np.array(dspear)
        donor_min = float(np.min(dspear)) if len(dspear) else np.nan
        donor_p10 = float(np.percentile(dspear, 10)) if len(dspear) else np.nan
        donor_median = float(np.median(dspear)) if len(dspear) else np.nan

        cls = classify(rho_obs, delta, r2, slope, sd_ratio)

        # worst-donor downgrade (§4.3): numeric but collapses in worst donor
        if cls == "numerically_transportable" and (donor_min < 0 or donor_p10 < 0.05):
            cls = "rank_invariant_but_uncalibrated"

        rows.append({
            "target": target, "mapping_class": MAPPING.mapping_class.iloc[j],
            "core_panel_role": MAPPING.core_panel_role.iloc[j],
            "spearman": rho_obs, "noncognate_spearman": rho_nc,
            "specificity_delta": delta, "r2": r2,
            "calibration_slope": slope, "sd_ratio": sd_ratio,
            "donor_median": donor_median, "donor_p10": donor_p10, "donor_min": donor_min,
            "classification": cls,
        })

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "P3_43TARGET_CLASSIFIED_MAIN_RESULT.csv", index=False)

    # summary
    summ = {
        "version": "p3_transport_expanded43_v2_classified",
        "direction": "GSE334503(discovery) -> GSE335494(validation)",
        "ridge_alpha": 1.0,
        "noncognate_seed": 20260902,
        "n_targets": int(len(out)),
        "classification_counts": out.classification.value_counts().to_dict(),
        "primary_direct_counts": out[out.core_panel_role.eq("primary")].classification.value_counts().to_dict(),
        "rule_ref": "TRANSPORTABILITY_RULE_FORMALIZATION_20260905.md §2.1",
    }
    (OUT / "P3_43TARGET_CLASSIFIED_SUMMARY.json").write_text(
        __import__("json").dumps(summ, indent=2), encoding="utf-8")

    print("=== classification counts (all 43) ===")
    print(out.classification.value_counts().to_string())
    print("\n=== primary_direct only ===")
    print(out[out.core_panel_role.eq("primary")].classification.value_counts().to_string())
    print("\n=== full table (rounded) ===")
    cols = ["target", "mapping_class", "spearman", "noncognate_spearman", "specificity_delta",
            "r2", "calibration_slope", "sd_ratio", "donor_median", "donor_min", "classification"]
    print(out[cols].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
