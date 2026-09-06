# -*- coding: utf-8 -*-
"""Update A_target_direction_cohort_master_audit.csv for the SD-ratio gate
reapplication (0.5-2.0).  Applies to BOTH identical copies (repo + evidence pkg).

Changes (from reapply_sd_ratio_gate.py output):
  CD16 : p3_classification numerically_transportable_both_sources -> rank_invariant_both_sources
         effective_boundary_tier CONDITIONAL_numeric_transport -> CONDITIONAL_rank_invariant
  CD3  : same
  CD314: p3_classification numerically_transportable_both_sources -> rank_invariant_both_sources
         canonical_classification numerically_transportable -> rank_invariant_but_uncalibrated
         effective_boundary_tier NUMERIC_BUT_MEASUREMENT_RISK -> FAIL_measurement
"""
import pandas as pd

PATHS = [
    r"D:\rna-adt-transportability-boundary\supplementary\A_target_direction_cohort_master_audit.csv",
    r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\supplementary_evidence_package_v1_20260906\A_target_direction_cohort_master_audit.csv",
]

UPD = {
    "CD16": {
        "p3_classification": "rank_invariant_both_sources",
        "effective_boundary_tier": "CONDITIONAL_rank_invariant",
    },
    "CD3": {
        "p3_classification": "rank_invariant_both_sources",
        "effective_boundary_tier": "CONDITIONAL_rank_invariant",
    },
    "CD314": {
        "p3_classification": "rank_invariant_both_sources",
        "canonical_classification": "rank_invariant_but_uncalibrated",
        "effective_boundary_tier": "FAIL_measurement",
    },
}

for path in PATHS:
    d = pd.read_csv(path)
    for tgt, updates in UPD.items():
        mask = d.target == tgt
        assert mask.sum() == 1, f"{path}: {tgt} matched {mask.sum()} rows"
        for col, val in updates.items():
            old = d.loc[mask, col].iloc[0]
            d.loc[mask, col] = val
            print(f"{path.split(chr(92))[-1][:12]}  {tgt:6s} {col:26s} {old} -> {val}")
    d.to_csv(path, index=False)
    print(f"  WROTE {path}")

# final verification
d = pd.read_csv(PATHS[0])
print("\n=== verification (repo) ===")
print(d[d.target.isin(["CD16", "CD3", "CD314"])][["target", "p3_classification", "canonical_classification", "effective_boundary_tier"]].to_string(index=False))
print("\neffective_boundary_tier counts:")
print(d.effective_boundary_tier.value_counts().to_string())
