import pandas as pd
import numpy as np

G = pd.read_csv(r"D:\rna-adt-transportability-boundary\supplementary\G_six_ICI_cohort_full_stats.csv")

p = G["permutation_p"].values
order = np.argsort(p)
n = len(p)
q = np.empty(n)
for i, idx in enumerate(order):
    q[idx] = p[idx] * n / (i + 1)
# enforce monotonic (BH step-up): take min over suffix
for i in range(n - 2, -1, -1):
    q[order[i]] = min(q[order[i]], q[order[i + 1]])

G["family_adjusted_p_bh"] = q
G["passes_family_q05"] = q < 0.05

# response coding column (R=CR/PR responder, NR=SD/PD non-responder, NA=NE not evaluable)
coding = {
    "GSE78220": "R=CR/PR, NR=SD/PD (binaryResponse); permutation stratified within cohort",
    "GSE91061": "R=CR/PR, NR=SD/PD (binaryResponse); permutation stratified within cohort",
    "GSE93157": "R=CR/PR, NR=SD/PD (binaryResponse); permutation stratified within cohort",
    "IMvigor210": "R=CR/PR, NR=SD/PD, NA=NE (Best.Confirmed.Overall.Response); permutation stratified within cohort",
    "Gide2019": "R=CR/PR, NR=SD/PD (RECIST); permutation stratified within cohort",
    "Braun2020": "R=CR/PR, NR=SD/PD (RECIST); permutation stratified within cohort",
}
G["response_coding"] = G["cohort"].map(coding)
G["permutation_stratification"] = "within-cohort response label shuffle (10,000 permutations, seed 42)"

cols = ["cohort","cancer","drug","n_evaluable","n_responders",
        "auroc","auroc_ci_lo","auroc_ci_hi","auprc","brier",
        "calibration_intercept","calibration_slope",
        "permutation_p","family_adjusted_p_bh","passes_family_q05",
        "response_coding","permutation_stratification"]

G = G[cols]
G.to_csv(r"D:\rna-adt-transportability-boundary\supplementary\G_six_ICI_cohort_full_stats.csv", index=False)

print(G.to_string(index=False))
print()
print("=== family-adjusted q (BH) ===")
print("only Gide passes q<0.05:", bool((G['passes_family_q05']).sum()==1 and G.loc[G.cohort=='Gide2019','passes_family_q05'].iloc[0]))
