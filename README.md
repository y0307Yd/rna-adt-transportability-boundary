# Cross-cohort transportability of RNA-to-protein prediction is target- and context-specific

**Reproducible analysis of cross-cohort RNA–ADT transportability and applicability boundaries**

This repository contains the analysis code and machine-readable evidence tables that
back the quantitative claims in the manuscript
*"Cross-cohort transportability of RNA-to-protein prediction is target- and context-specific"*
(current manuscript: `manuscript_final_v8i_submission_ready.docx`).

## What this study is about

We ask a deceptively simple question: *when a model predicts protein abundance from
RNA in held-out cells, does that prediction transport across cohorts and contexts as if
it were a measurement?* The answer, for 51 antibody targets evaluated across 5 CITE-seq
cohorts and 6 clinical ICI cohorts, is **no — and the boundary is target- and
context-specific**.

We replace the binary "works / does not work" framing with an **evidence hierarchy**
for external validity:

1. **Numerical transport** — predicted and measured protein agree in scale, direction, and variance.
2. **Rank-only transport** — only the rank order is preserved; scale/calibration is lost.
3. **Context-dependent** — transport holds in some contexts but not others.
4. **Weak / non-transportable** — no usable signal.
5. **Measurement / hard failure** — the target fails as a measured protein surrogate in *every evaluated setting*.

The central message: **a model that predicts protein well in held-out cells is not
necessarily a model that measures protein across cohorts.** Transportability is an
empirical property of the *target–context–measurement system*, not merely an algorithm
property.

## Three principal results

- **Within-source donor generalization is not cross-cohort transport.** The strict
  five-component gate was attainable in 4/22 source-internal donor-held-out
  target-cohort settings, with 0/22 shuffled controls passing, but no target passed
  the same strict numerical gate in the principal expanded external transfer.
- **High rank correlation is not numerical interchangeability.** In the principal
  43-target external analysis, 29 targets retained rank information without adequate
  numerical calibration, 3 were weakly context-dependent and 11 were non-specific or
  non-transportable under the strict gate.
- **Technical repeatability does not guarantee transportability.** Of eight targets
  classified as failures across all evaluated transport settings with available
  orthogonal ICC values, four had ICC >= 0.50. Those four failures cannot be dismissed
  as low repeatability alone; twelve separate targets failed at the measurement layer.

## Repository structure

```
rna-adt-transportability-boundary/
├── README.md
├── LICENSE                      # MIT
├── requirements.txt             # Python dependencies
├── .gitignore
├── supplementary/               # Machine-readable evidence tables (A–H, including B3) + index README
└── scripts/                     # Core reproducible analysis scripts
    ├── target_scope/            # 43/51-target transport, direction-dependence, LOPO fold metadata
    ├── celltype/                # Cell-type-adjusted transport (Figure 13)
    ├── spatial/                 # GSE289326 IMC spatial protein validation
    ├── technical_variance/      # GSE245108 titration variance decomposition + identifiability
    ├── clinical/                # Frozen four-gene score validation (IMvigor210, Gide2019, etc.)
    ├── external_lopo/           # External LOPO dataset construction (GSE305118/334503/335494)
    └── figures/                 # Key figure generation scripts
```

## Supplementary evidence package

The `supplementary/` directory contains machine-readable tables, each mapped to a
manuscript section (see `supplementary/README.md` for the full index). These are the
tables referenced by the *"Supplementary materials comprise (A)–(G)"* sentence in the
Data and code availability section:

| Letter | File | Manuscript section |
|---|---|---|
| A | `A_target_direction_cohort_master_audit.csv` | 3.2 / 3.4 |
| A2 | `A2_43target_bidirectional_evidence_matrix.csv` | 2.7 / 3.2 |
| B | `B_fold_level_model_table.csv` | 2.7 (LODO/LOPO) |
| B3 | `B3_source_lodo_target_metrics.csv` | 2.6 / 3.1 (source-internal strict-gate positive control) |
| B3 donor | `B3_source_lodo_donor_target_metrics.csv` | 2.6 / 3.1 (121 held-out donor-target rows) |
| B3 mapping | `B3_null_b_noncognate_mapping.csv` | 2.6 (fixed Null-B mapping) |
| B3 audit | `B3_legacy_fold_consistency_audit.csv` | 2.6 (legacy consistency) |
| B3 manifest | `B3_input_manifest_sha256.csv` | 2.6 (input hashes) |
| C1 | `C1_spatial_sign_flip_bh_qvalues.csv` | 3.6 (spatial) |
| C2 | `C2_spatial_block_rho_bootstrap_ci.csv` | 3.6 (spatial) |
| D | `D_gse245108_sample_design_matrix.csv` | 3.5 (technical variance) |
| D2 | `D2_gse245108_variance_components_all.csv` | 3.5 |
| D3 | `D3_convergence_bias_features.csv` | 3.5 |
| D4 | `D4_convergence_bias_tests.csv` | 3.5 |
| E | `E_frozen_evidence_timeline_audit.md` | 2.13 |
| G | `G_six_ICI_cohort_full_stats.csv` | 3.7 (clinical) |
| H1 | `H1_GSE334503_celltype_adjusted.csv` | 3.2 / Figure 13 |
| H2 | `H2_GSE335494_celltype_adjusted.csv` | 3.2 / Figure 13 |

## Methodological boundary statement (important)

All analysis rules (target panel, transport thresholds, classification gates, frozen
clinical score) are **analysis-defined and retrospectively locked**, not pre-registered.
A local versioned contract (LODO/LOPO fold specification, see `supplementary/E_`) fixed
the evaluation rules *before* the corresponding results were generated, but no public
pre-registration or immutable release hash exists. We therefore describe boundaries as
**observed / empirical**, not as prospectively validated claims.

The frozen four-gene clinical score (CD8A + HLA-DRA + PDCD1 + ITGAE) is a
**deliberately frozen falsification test**, not a clinical application. It is *not* a
validated predictive biomarker: it lacks multivariable adjustment, independent
validation of incremental value, calibration, and treatment-interaction evidence.
Across six ICI cohorts (n = 27/49/65/298/73/172), only one cohort (Gide2019, n = 73)
showed nominal significance, which did not survive as a four-gene-specific discriminator
after considering that broad immune-infiltration markers performed similarly.

Spatial validation (12 IMC sections, 6 patients, GSE289326) returned **no target passing
the family-level gate** (minimum BH-q = 0.109, all > 0.05). It is reported as a negative
orthogonal check, not as cell-level transport evidence.

## Reproducing the analysis

1. Install dependencies: `pip install -r requirements.txt`
2. Download the source datasets from GEO (accessions listed in the manuscript and in
   `supplementary/README.md`).
3. Run scripts under `scripts/` in the order implied by their directory names.

> **Note on paths.** Scripts contain absolute paths from the original analysis
> environment (Windows). They are provided for transparency and audit, not as a
> turnkey pipeline. The `supplementary/` tables already contain the final numerical
> results; the scripts document exactly how each was produced.

## Key analysis parameters (extracted from source)

- **43-target transport** (`scripts/target_scope/run_expanded43_transport.py`):
  Ridge α = 1.0; feature = mapped cognate RNA only; source-fitted StandardScaler;
  seed 20260902; shuffle seed = 20260902 + target index.
- **Spatial validation** (`scripts/spatial/run_gse289326_spatial_protein_validation.py`):
  seed 20260820; N_PERM = 1000; 6 patients × 12 sections (pre/post); technical
  replicates excluded; CD3E anchor threshold = quantile 0.75; permutation
  P = (1 + Σ[|null| ≥ |ρ|]) / (N_PERM + 1).
- **GSE245108 variance decomposition** (`scripts/technical_variance/gse245108/`):
  per-target mixed model `ADT_CLR ~ 1 + (1|donor) + (1|cell_type) + (1|concentration)
  + (1|rep)`; statsmodels MixedLM, REML, lbfgs, maxiter 300; variance =
  f.params["<group> Var"] × f.scale; CLR = log1p(counts) − row mean.
- **Cell-type adjusted** (`scripts/celltype/run_celltype_adjusted.py`): pooled /
  within-cell-type / cell-type-adjusted / baseline / per-donor Spearman; "Ambiguous"
  cells excluded from cell-type-aware analyses.
- **Clinical frozen score**: score = mean(log1p) of {CD8A, HLA-DRA, PDCD1, ITGAE};
  10000 response shuffles; Brier/calibration via logistic link; seed 42.

## License

Code is released under the MIT License (see `LICENSE`). Source datasets are governed by
their respective GEO / authors' terms.

## Contact

For questions about the analysis, open an issue on this repository.
