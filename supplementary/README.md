# Supplementary Evidence Package — v1 (2026-09-06)

Machine-readable tables backing every quantitative claim in the manuscript.
Each letter maps to the "Supplementary materials comprise (A)–(G)" sentence in the
Data and code availability section.

## File index

| Letter | File | Content | Manuscript section |
|---|---|---|---|
| A | `A_target_direction_cohort_master_audit.csv` | 51-target master boundary: transport class, measurement flag, mapping class, Kotliarov ICC, effective boundary tier | 3.2 / 3.4 |
| A2 | `A2_43target_bidirectional_evidence_matrix.csv` | 43-target bidirectional (forward/reverse) evidence: observed vs shuffled median, delta, bootstrap CI, sign-flip P, BH q, classification | 2.7 / 3.2 |
| B | `B_fold_level_model_table.csv` | LOPO fold-level metadata: held-out target, training targets, selected features, scaler mean/scale, ridge alpha/solver, observed vs predicted | 2.7 (LODO/LOPO) |
| C1 | `C1_spatial_sign_flip_bh_qvalues.csv` | Spatial sign-flip tests: analysis, target, timepoint, median, raw P, grade, BH q, passes-q05 | 3.6 (spatial) |
| C2 | `C2_spatial_block_rho_bootstrap_ci.csv` | Spatial block bootstrap CI: target, timepoint, n sections, observed median, CI lo/hi | 3.6 (spatial) |
| D | `D_gse245108_sample_design_matrix.csv` | GSE245108 sample-level design: donor × concentration × rep, n cells | 3.5 (technical variance) |
| D2 | `D2_gse245108_variance_components_all.csv` | Full 265-target variance components (donor/cell_type/concentration/rep/residual + fractions + converged flag) | 3.5 |
| D3 | `D3_convergence_bias_features.csv` | Per-feature convergence-bias covariates (mean/sd/cv, donor heterogeneity, specificity, fractions) | 3.5 |
| D4 | `D4_convergence_bias_tests.csv` | Converged vs non-converged Mann-Whitney tests | 3.5 |
| E | `E_frozen_evidence_timeline_audit.md` | Frozen-evidence timeline: which rules preceded which results (local mtime evidence) | 2.13 |
| G | `G_six_ICI_cohort_full_stats.csv` | Six ICI cohort full stats: AUROC + 95% CI, AUPRC, Brier, calibration intercept/slope, permutation P | 3.7 (clinical) |
| H1 | `H1_GSE334503_celltype_adjusted.csv` | Cell-type adjusted transport (GSE334503): pooled/within/ct-adjusted/baseline/donor p10/min | 3.2 / Figure 13 |
| H2 | `H2_GSE335494_celltype_adjusted.csv` | Cell-type adjusted transport (GSE335494) | 3.2 / Figure 13 |

## Key analysis parameters (extracted from source scripts)

- **43-target transport** (`run_expanded43_transport.py`): Ridge alpha=1.0; feature = mapped cognate RNA only; source-fitted StandardScaler; seed 20260902; evaluation = pooled Spearman, donor-median correlation, R², MAE, RMSE, calibration slope+intercept, predicted-vs-observed SD ratio; shuffle seed = 20260902 + target index.
- **Spatial validation** (`run_gse289326_spatial_protein_validation.py`): seed 20260820; N_PERM=1000; 6 patients × 12 sections (pre/post); technical replicates GSM8789211/GSM8789214 excluded; CD3E anchor threshold = quantile 0.75 (≥25 anchors); permutation P = (1 + Σ[|null| ≥ |rho|]) / (N_PERM + 1); grades = PAIRED-SPATIAL-PROTEIN GO/EXPLORATORY/NO-GO (sign-flip exact test + ≥4/6 agreement + |r| ≥ 0.10); CD103 locked as structurally unavailable.
- **Spatial cd3e-conditional bootstrap** (`manifest.json`): unit = section; statistic = median block partial rho; n_bootstrap = 10000; seed 20260901.
- **GSE245108 variance decomposition** (`variance_decomp_step4_all.py`): per-target model `ADT_CLR ~ 1 + (1|donor) + (1|cell_type) + (1|concentration) + (1|rep)`; statsmodels MixedLM, REML, method lbfgs, maxiter 300; variance = f.params["<group> Var"] × f.scale; CLR = log1p(counts) − row mean.
- **Cell-type adjusted** (`run_celltype_adjusted.py`): pooled = Spearman(log1p(cognate RNA), CLR(ADT)) over all cells; within = per-cell-type Spearman (n≥100), median across types; adjusted = Spearman after removing cell-type mean; baseline = cell-type one-hot mean ADT; donor = median/p10/min of per-donor Spearman (n≥50 cells); "Ambiguous" cells excluded from cell-type-aware analyses.
- **Clinical frozen score** (`validate_frozen_score.py` etc.): score = mean(log1p) of {CD8A, HLA-DRA, PDCD1, ITGAE}; permutation = 10000 response shuffles; Brier/calibration via logistic link; Gide P uses (cnt+1)/(nperm+1); seed 42.
