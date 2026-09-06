#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Gide 2019 (PRJEB23709) - HLA-DRA a priori + frozen 4-gene score validation
冻结四基因 score = mean(log1p(normalized_counts)) over {CD8A, HLA-DRA, PDCD1, ITGAE}
HLA-DRA 单分析物 = log1p(HLA-DRA)  (a priori hypothesis, pre-registered intent)
不重训、不调阈值；AUROC + 置换检验(10000)
"""
import sys, csv, json
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

d = r"D:/ai_multimodal_cholesterol_study_outputs/gide2019"
GENES = ["CD8A", "HLA-DRA", "PDCD1", "ITGAE"]

# ---- 1. read counts matrix ----
with open(d + "/cancercell_normalized_counts_genenames.txt", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f, delimiter="\t")
    header = next(rd)
    all_rows = list(rd)

samples = header[3:]  # 91 sample columns
gene_symbol = {r[1]: i for i, r in enumerate(all_rows)}

# ---- 2. read response annotation (coldata_example.rda via pyreadr) ----
import bz2
import pyreadr
coldata = pyreadr.read_r(d + "/coldata_example.rda")["coldata_example"]
# rownames like "Sample_10_PD1_PRE" -> map to column "PD1_10_PRE"
def rowname_to_col(name):
    # Sample_{N}_{trt}_{time} -> {trt}_{N}_{time}
    parts = name.split("_")
    # parts[0]=Sample, parts[1]=N, parts[2]=trt(PRE-final), ...
    # e.g. Sample_10_PD1_PRE -> PD1_10_PRE
    num = parts[1]
    # treatment may be "PD1" or "ipiPD1"
    if parts[2] == "PD1":
        trt = "PD1"
        time = parts[3]
    else:
        trt = "ipiPD1"
        time = parts[3]
    return f"{trt}_{num}_{time}"

coldata["colname"] = [rowname_to_col(rn) for rn in coldata.index]
coldata["treatment"] = [("PD1" if c.startswith("PD1_") else "ipiPD1") for c in coldata["colname"]]

# ---- 3. build expression matrix for the 73 PRE samples with response ----
sample_to_col = {s: i for i, s in enumerate(samples)}
expr = {}
for g in GENES:
    if g not in gene_symbol:
        print(f"FATAL: {g} missing")
        sys.exit(1)
    r = all_rows[gene_symbol[g]]
    vals = [float(x) for x in r[3:]]
    expr[g] = vals

meta = []
for _, row in coldata.iterrows():
    col = row["colname"]
    if col not in sample_to_col:
        print(f"WARN: col {col} not in counts, skip")
        continue
    ci = sample_to_col[col]
    rec = {
        "sample": col,
        "treatment": row["treatment"],
        "response": row["Response"],  # "R" or "NR"
        "CD8A": expr["CD8A"][ci],
        "HLA_DRA": expr["HLA-DRA"][ci],
        "PDCD1": expr["PDCD1"][ci],
        "ITGAE": expr["ITGAE"][ci],
    }
    meta.append(rec)

print(f"Total PRE samples with response: {len(meta)}")
from collections import Counter
print("by treatment:", Counter(m["treatment"] for m in meta))
print("by response:", Counter(m["response"] for m in meta))
trt_resp = Counter((m["treatment"], m["response"]) for m in meta)
print("treatment x response:", dict(trt_resp))

# ---- 4. compute scores ----
for m in meta:
    g = [m["CD8A"], m["HLA_DRA"], m["PDCD1"], m["ITGAE"]]
    m["score_log1p"] = float(np.mean(np.log1p(g)))
    m["HLA_DRA_log1p"] = float(np.log1p(m["HLA_DRA"]))

y = np.array([1 if m["response"] == "R" else 0 for m in meta])
score = np.array([m["score_log1p"] for m in meta])
hla = np.array([m["HLA_DRA_log1p"] for m in meta])

# ---- 5. AUROC + permutation ----
def auc(y, s):
    if len(np.unique(y)) < 2:
        return np.nan
    order = np.argsort(s)
    y_sorted = y[order]
    n_pos = y_sorted.sum()
    n_neg = len(y_sorted) - n_pos
    ranks = np.arange(1, len(y_sorted) + 1)
    # handle ties: average rank
    # simple Mann-Whitney via pairwise
    pos = s[y == 1]
    neg = s[y == 0]
    # U statistic
    n1, n2 = len(pos), len(neg)
    if n1 == 0 or n2 == 0:
        return np.nan
    # count pos > neg + 0.5 ties
    # use ranksum
    from scipy.stats import mannwhitneyu
    u, _ = mannwhitneyu(pos, neg, alternative="two-sided")
    return u / (n1 * n2)

def auc_perm(y, s, nperm=10000, seed=42):
    rng = np.random.default_rng(seed)
    obs = auc(y, s)
    cnt = 0
    for _ in range(nperm):
        yp = rng.permutation(y)
        if auc(yp, s) >= obs:
            cnt += 1
    p = (cnt + 1) / (nperm + 1)
    return obs, p

def run(name, s, subset=None):
    yy, ss = y, s
    if subset is not None:
        yy = y[subset]
        ss = s[subset]
    a, p = auc_perm(yy, ss)
    print(f"[{name}] n={len(yy)} R={int(yy.sum())} NR={int(len(yy)-yy.sum())} AUROC={a:.4f} perm_p={p:.4f}")
    return a, p

results = {}
print("\n===== FULL cohort (all 73, PD1 + ipiPD1) =====")
results["full_4gene"] = run("frozen 4-gene score", score)
results["full_HLA"] = run("HLA-DRA single analyte", hla)

print("\n===== PD1 monotherapy only (41) =====")
pd1_mask = np.array([m["treatment"] == "PD1" for m in meta])
results["pd1_4gene"] = run("frozen 4-gene (PD1)", score, pd1_mask)
results["pd1_HLA"] = run("HLA-DRA (PD1)", hla, pd1_mask)

print("\n===== ipiPD1 combo only (32) =====")
combo_mask = ~pd1_mask
results["combo_4gene"] = run("frozen 4-gene (combo)", score, combo_mask)
results["combo_HLA"] = run("HLA-DRA (combo)", hla, combo_mask)

# ---- 6. also per-gene t-test (responder vs non-responder) ----
from scipy.stats import ttest_ind
print("\n===== per-gene t-test (log1p normalized counts, R vs NR, full cohort) =====")
for g in ["CD8A", "HLA_DRA", "PDCD1", "ITGAE"]:
    gv = np.log1p(np.array([m[g] for m in meta]))
    r = gv[y == 1]; nr = gv[y == 0]
    t, p = ttest_ind(r, nr)
    print(f"{g}: R mean={r.mean():.3f} NR mean={nr.mean():.3f} t={t:+.2f} p={p:.4f}")

# ---- save ----
out_rows = []
for m in meta:
    out_rows.append({
        "sample": m["sample"], "treatment": m["treatment"],
        "response": m["response"], "score_log1p": m["score_log1p"],
        "HLA_DRA_log1p": m["HLA_DRA_log1p"],
        "CD8A": m["CD8A"], "HLA_DRA": m["HLA_DRA"],
        "PDCD1": m["PDCD1"], "ITGAE": m["ITGAE"],
    })
with open(d + "/gide2019_frozen_score_with_response.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=out_rows[0].keys())
    w.writeheader()
    w.writerows(out_rows)

summary = {
    "cohort": "Gide 2019 (PRJEB23709, Cancer Cell 2019)",
    "cancer": "metastatic melanoma",
    "drug": "anti-PD-1 (pembro/nivo) monotherapy + anti-PD-1/anti-CTLA-4 combo",
    "n_PRE_with_response": len(meta),
    "n_PD1": int(pd1_mask.sum()), "n_combo": int(combo_mask.sum()),
    "response_criteria": "RECIST 1.1",
    "score_definition": "frozen mean log1p normalized_counts of {CD8A,HLA-DRA,PDCD1,ITGAE}",
    "HLA_DRA_hypothesis": "a priori: HLA-DRA single analyte predicts response",
    "results": results,
}
with open(d + "/gide2019_frozen_score_summary.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)
print("\nSaved CSV + JSON.")
