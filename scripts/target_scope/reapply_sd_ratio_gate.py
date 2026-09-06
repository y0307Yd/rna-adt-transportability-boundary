# -*- coding: utf-8 -*-
"""
reapply_sd_ratio_gate.py
========================
Re-apply the SD-ratio primary gate (0.5-2.0, formerly 0.05-20) to the canonical
43-target classification, the GSE314416 cross-source classification, and the
51-target master boundary, WITHOUT re-running any underlying transport model.

This is the reproducible implementation of the Issue-4 reviewer request:
  "numerical transport thresholds lack external justification (SD ratio 0.05-20
   is extremely wide)".  Decision (user-confirmed 2026-09-06):
    - primary gate:  0.5 <= sd_ratio <= 2.0   (contributes to classification)
    - exploratory:   0.05 <= sd_ratio <= 20   (threshold-sensitivity only)

Read-only inputs; writes to a NEW versioned output directory (no baseline
overwrite).
"""
from pathlib import Path
import numpy as np
import pandas as pd
import json, re
from collections import Counter

BASE = Path(r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap")
OUT = BASE / "sd_ratio_gate_reeval_v1_20260906"
OUT.mkdir(parents=True, exist_ok=True)

SD_LO_OLD, SD_HI_OLD = 0.05, 20.0
SD_LO_NEW, SD_HI_NEW = 0.5, 2.0

# ---------------------------------------------------------------- 43-target canonical
D43 = pd.read_csv(BASE / r"p3_target_scope_20260902\p3_transport_expanded43_v2_classified\P3_43TARGET_CLASSIFIED_MAIN_RESULT.csv")


def classify43(rho, delta, r2, slope, sd, lo, hi):
    rank = (rho >= 0.20) and (delta >= 0.10)
    numeric = rank and (r2 > 0) and (0.5 <= slope <= 2.0) and (lo <= sd <= hi)
    if numeric:
        return "numerically_transportable"
    if rank:
        return "rank_invariant_but_uncalibrated"
    if rho >= 0.10:
        return "weak_context_dependent"
    return "non_specific_or_non_transportable"


old43 = D43.apply(lambda r: classify43(r.spearman, r.specificity_delta, r.r2, r.calibration_slope, r.sd_ratio, SD_LO_OLD, SD_HI_OLD), axis=1)
new43 = D43.apply(lambda r: classify43(r.spearman, r.specificity_delta, r.r2, r.calibration_slope, r.sd_ratio, SD_LO_NEW, SD_HI_NEW), axis=1)
D43_new = D43.copy()
D43_new["classification_old"] = old43
D43_new["classification_new"] = new43
D43_new.to_csv(OUT / "P3_43TARGET_CLASSIFIED_SD_RATIO_REAPPLIED.csv", index=False)

# ---------------------------------------------------------------- GSE314416 cross-source
M = pd.read_csv(BASE / r"p3_third_complete_panel_search_20260903\GSE314416\transport_v5_quality_gate\GSE314416_V4_METRICS.csv")
p = M[M.n_participants.notna()].copy()
obs = p[p.control.eq("cognate")].copy()
neg = p[p.control.eq("noncognate")][["source", "target", "spearman"]].rename(columns={"spearman": "noncognate_spearman"})
x = obs.merge(neg, on=["source", "target"], how="left")
x["specificity_delta"] = x.spearman - x.noncognate_spearman


def classify_gs(df, lo, hi):
    df = df.copy()
    df["variance_gate"] = df.prediction_observed_sd_ratio.between(lo, hi)
    df["rank_gate"] = (df.spearman >= 0.20) & (df.specificity_delta >= 0.10)
    df["numeric_gate"] = (df.r2 > 0) & df.calibration_slope.between(0.5, 2.0) & df.variance_gate
    c = df.groupby(["target", "analysis_role"]).agg(pass_rank=("rank_gate", "sum"), pass_numeric=("numeric_gate", "sum")).reset_index()
    c["cls"] = np.select(
        [(c.pass_numeric == 2) & (c.pass_rank == 2), c.pass_rank == 2, c.pass_rank == 1],
        ["numerically_transportable_both_sources", "rank_invariant_both_sources", "source_dependent_rank_signal"],
        default="weak_or_nontransportable")
    return c


c_old = classify_gs(x, SD_LO_OLD, SD_HI_OLD)
c_new = classify_gs(x, SD_LO_NEW, SD_HI_NEW)
c_new.to_csv(OUT / "GSE314416_CROSS_SOURCE_CLASSIFICATION_SD_RATIO_REAPPLIED.csv", index=False)


def canon(name):
    n = str(name).strip()
    base = n.split("(")[0].strip()
    base = re.split(r"\s+", base)[0].strip()
    return base


new_cls_map = {canon(t): cls for t, cls in zip(c_new.target, c_new.cls)}

# ---------------------------------------------------------------- 51-target master boundary
B = pd.read_csv(BASE / r"p2_p3_unified_target_boundary_v3_20260904\unified_43target_boundary_v3_20260904.csv")


def icc_float(v):
    if v is None:
        return None
    if isinstance(v, float) and np.isnan(v):
        return None
    s = str(v).strip()
    if s in ("", "nan", "NaN", "None"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def derive_tier(cls, tech_class, icc):
    low = (tech_class == "low") or (icc is not None and icc < 0.3)
    if cls == "numerically_transportable_both_sources":
        return "NUMERIC_BUT_MEASUREMENT_RISK" if (low or tech_class == "unknown") else "CONDITIONAL_numeric_transport"
    if cls == "rank_invariant_both_sources":
        if low:
            return "FAIL_measurement"
        if tech_class == "unknown":
            return "CONDITIONAL_rank_invariant_measurement_unknown"
        return "CONDITIONAL_rank_invariant"
    if cls == "weak_or_nontransportable":
        return "FAIL_universal"
    if low:
        return "FAIL_measurement"
    if tech_class == "unknown":
        return "INSUFFICIENT_evidence"
    return "INSUFFICIENT_transport_evidence"


def eff_tier(cls, final_tier, icc, tech_class):
    mf = (icc is not None and icc < 0.3) or tech_class == "low"
    if mf:
        return "NUMERIC_BUT_MEASUREMENT_RISK" if cls == "numerically_transportable_both_sources" else "FAIL_measurement"
    return final_tier


rows = []
for _, r in B.iterrows():
    tgt = r["target"]
    ckey = canon(tgt)
    old_cls = r["p3_classification"]
    new_cls = new_cls_map.get(ckey, old_cls) if pd.notna(old_cls) else new_cls_map.get(ckey, "")
    if pd.isna(new_cls):
        new_cls = ""
    icc = icc_float(r["kotliarov_ICC"])
    raw_tc = r["tech_repeatability_class"]
    tcls = str(raw_tc).strip() if pd.notna(raw_tc) and str(raw_tc).strip() not in ("", "nan") else "unknown"
    raw_p2 = r["p2_failure_tier"]
    p2_tier = str(raw_p2).strip() if pd.notna(raw_p2) else ""
    derived = derive_tier(new_cls, tcls, icc)
    final_tier = p2_tier if p2_tier else derived
    eff = eff_tier(new_cls, final_tier, icc, tcls)
    rows.append({"target": tgt, "p3_classification_old": old_cls, "p3_classification_new": new_cls,
                 "effective_boundary_tier_old": r["effective_boundary_tier"], "effective_boundary_tier_new": eff})

NR = pd.DataFrame(rows)
NR.to_csv(OUT / "MASTER_BOUNDARY_TIER_REAPPLIED.csv", index=False)

# ---------------------------------------------------------------- summary
old_cnt = Counter(NR.effective_boundary_tier_old.astype(str))
new_cnt = Counter(NR.effective_boundary_tier_new.astype(str))
changed = NR[NR.effective_boundary_tier_old != NR.effective_boundary_tier_new]
hard_old = old_cnt.get("FAIL_universal", 0) + old_cnt.get("FAIL_measurement", 0)
hard_new = new_cnt.get("FAIL_universal", 0) + new_cnt.get("FAIL_measurement", 0)

summary = {
    "version": "sd_ratio_gate_reeval_v1_20260906",
    "decision": "primary gate 0.5-2.0; exploratory 0.05-20 (sensitivity only)",
    "gse334503_to_gse335494_43target": {
        "old": {k: int(v) for k, v in old43.value_counts().items()},
        "new": {k: int(v) for k, v in new43.value_counts().items()},
        "changed": [(D43.target.iloc[i], old43.iloc[i], new43.iloc[i], round(D43.sd_ratio.iloc[i], 4))
                    for i in range(len(D43)) if old43.iloc[i] != new43.iloc[i]],
    },
    "gse314416_cross_source": {
        "old": {k: int(v) for k, v in c_old.cls.value_counts().items()},
        "new": {k: int(v) for k, v in c_new.cls.value_counts().items()},
        "changed": [(row.target, row.cls_old, row.cls_new) for _, row in
                    c_old[["target", "cls"]].merge(c_new[["target", "cls"]], on="target", suffixes=("_old", "_new"))
                    .query("cls_old != cls_new").iterrows()],
        "primary_direct_new": {k: int(v) for k, v in
                               c_new[c_new.analysis_role.eq("primary_direct")].cls.value_counts().items()},
    },
    "master_boundary_effective_tier": {
        "old": {k: int(v) for k, v in old_cnt.items()},
        "new": {k: int(v) for k, v in new_cnt.items()},
        "hard_fail_old": int(hard_old),
        "hard_fail_new": int(hard_new),
        "changed_targets": [(r.target, r.effective_boundary_tier_old, r.effective_boundary_tier_new)
                            for _, r in changed.iterrows()],
    },
}
(OUT / "SD_RATIO_GATE_REEVAL_SUMMARY.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

print(json.dumps(summary, indent=2, ensure_ascii=False))
print("\nWROTE to", OUT)
