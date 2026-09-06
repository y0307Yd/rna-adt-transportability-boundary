"""Convergence bias analysis (v2): are non-converged targets systematically different?

Alignment: CLR matrix (31851 x 265) columns <-> tnc_feature_names.json (265 unique names).
vc has 273 rows = 265 unique + 8 duplicated features (each appearing twice).
For duplicated features, converged = True if ANY row converged (and we flag ambiguity).
"""
from pathlib import Path
import numpy as np
import pandas as pd
import json
from scipy.stats import mannwhitneyu

BASE = Path(r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\technical_biological_decomposition_data_v1_20260903\GSE245108_tnc_parsed")
OUT = BASE / "convergence_bias_20260905"
OUT.mkdir(parents=True, exist_ok=True)

z = np.load(BASE / "tnc_adt_clr.npz", allow_pickle=True)
key = [k for k in z.files if k != "allow_pickle"][0]
clr = z[key].astype(np.float64)  # (31851, 265)
del z

fn = json.load(open(BASE / "tnc_feature_names.json", encoding="utf-8"))  # 265
meta = json.load(open(BASE / "tnc_meta_full.json", encoding="utf-8"))  # 31851
donors = np.array([m["donor"] for m in meta])
clusters = np.array([m["cluster"] for m in meta])

vc = pd.read_csv(BASE / "tnc_variance_components_all.csv")

# build per-unique-feature converged + frac_donor + frac_rep (resolve duplicates)
conv_map = {}
fracd_map = {}
fracr_map = {}
fracc_map = {}
dup_flag = {}
for f in set(fn):
    sub = vc[vc.feature == f]
    conv_map[f] = bool(sub.converged.any())
    fracd_map[f] = float(sub.frac_donor.median())
    fracr_map[f] = float(sub.frac_rep.median())
    fracc_map[f] = float(sub.frac_concentration.median())
    dup_flag[f] = len(sub) > 1

print(f"clr shape {clr.shape}, fn {len(fn)}, unique features {len(set(fn))}")

uniq_donors = np.unique(donors)
uniq_clusters = np.unique(clusters)

rows = []
for i, f in enumerate(fn):
    c = clr[:, i]
    conv = conv_map[f]
    mean_val = float(c.mean())
    sd_val = float(c.std())
    cv = float(sd_val / abs(mean_val)) if abs(mean_val) > 1e-9 else np.nan
    # donor-level heterogeneity (of CLR signal)
    dmean = np.array([c[donors == d].mean() for d in uniq_donors])
    dsd = float(dmean.std())
    dcv = float(dmean.std() / abs(dmean.mean())) if abs(dmean.mean()) > 1e-9 else np.nan
    # cell-type specificity
    cmean = np.array([c[clusters == cl].mean() for cl in uniq_clusters])
    specificity = float(cmean.max() / cmean.mean()) if cmean.mean() > 0 else np.nan
    frac_top = float(cmean.max() / cmean.sum()) if cmean.sum() > 0 else np.nan
    rows.append({
        "feature": f, "converged": conv, "duplicated": dup_flag[f],
        "mean_clr": mean_val, "sd_clr": sd_val, "cv": cv,
        "donor_mean_sd": dsd, "donor_mean_cv": dcv,
        "specificity": specificity, "frac_signal_top_cluster": frac_top,
        "frac_donor": fracd_map[f], "frac_rep": fracr_map[f],
        "frac_concentration": fracc_map[f],
    })

feats = pd.DataFrame(rows)

tests = []
for col in ["mean_clr", "sd_clr", "cv", "donor_mean_sd", "donor_mean_cv",
            "specificity", "frac_signal_top_cluster"]:
    a = feats.loc[feats.converged, col].dropna()
    b = feats.loc[~feats.converged, col].dropna()
    if len(a) >= 5 and len(b) >= 5:
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        tests.append({"feature": col,
                      "converged_median": float(a.median()),
                      "nonconverged_median": float(b.median()),
                      "mannwhitney_p": float(p)})
    else:
        tests.append({"feature": col, "converged_median": np.nan,
                      "nonconverged_median": np.nan, "mannwhitney_p": np.nan})

tdf = pd.DataFrame(tests)
feats.to_csv(OUT / "convergence_bias_features.csv", index=False)
tdf.to_csv(OUT / "convergence_bias_tests.csv", index=False)

print("=== n converged / non-converged ===")
print(feats.converged.value_counts().to_string())
print()
print("=== Mann-Whitney: converged vs non-converged feature medians ===")
print(tdf.round(4).to_string(index=False))
print()
print("=== point-biserial corr(converged, x) ===")
for col in ["frac_donor", "frac_rep", "frac_concentration", "mean_clr", "cv", "specificity", "donor_mean_sd"]:
    sub = feats[[col, "converged"]].dropna()
    if len(sub) > 10 and sub[col].nunique() > 1:
        r = float(np.corrcoef(sub[col], sub["converged"].astype(float))[0, 1])
        print(f"  corr(converged, {col:20s}) = {r:+.3f}")

# save summary json
summary = {
    "n_converged": int(feats.converged.sum()),
    "n_nonconverged": int((~feats.converged).sum()),
    "tests": tdf.to_dict("records"),
}
(OUT / "convergence_bias_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
