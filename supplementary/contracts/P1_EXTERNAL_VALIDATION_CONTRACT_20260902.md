# P1 external LOPO validation contract

## Objective

Prospectively test the frozen LOPO target-prioritization rule in external RNA-ADT paired cohorts without re-selecting features or tuning thresholds.

## Frozen rule

- Features: mean cognate RNA-ADT coupling and mean within-cohort LODO correlation.
- Model: ridge regression, alpha=1, source-only standardization.
- Held-out unit: target protein; each external cohort is evaluated without using its target labels to select features.
- Primary metrics: Spearman correlation and MAE between predicted and observed target transportability.
- Secondary output: target-level transport class with unavailable ADTs explicitly excluded.

## Candidate cohorts

- GSE316782 NSCLC: 9 donors, 11-target panel after naming correction.
- GSE317605 biliary cancer: 18 donors, 9-target panel; HLA-DR and TIGIT unavailable.

These cohorts are external stress-test resources. They are not assumed to be complete-panel confirmation until raw paired matrices and donor metadata pass the input gate.

## Required controls

Training-mean, donor-only, batch-only and composition-only baselines where metadata permit. No baseline result may overwrite existing transport outputs.

## Output policy

All P1 files must remain under `D:/ai_multimodal_cholesterol_study_outputs/revision-roadmap/p1_external_lopo_validation/`. Existing baseline directories are immutable.
