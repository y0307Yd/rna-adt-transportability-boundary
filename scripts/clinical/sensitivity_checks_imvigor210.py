# -*- coding: utf-8 -*-
"""Sensitivity checks for IMvigor210 frozen-score result."""
import csv, sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
rng = np.random.default_rng(20260904)

D = r"D:/ai_multimodal_cholesterol_study_outputs/imvigor210/csv"

with open(D + "/fData_IMvigor210.csv", encoding="utf-8", errors="replace") as f:
    fr = list(csv.reader(f))
f_rows = fr[1:]
sym_to_len = {r[1]: float(r[3]) if r[3] else 0.0 for r in f_rows}
sym_to_idx = {r[1]: i for i, r in enumerate(f_rows)}

with open(D + "/exmat_censored_IMvigor210.csv", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f)
    header = next(rd)
    sample_ids = [c.strip('"') for c in header[1:]]
    counts = np.array([[float(x) for x in row[1:]] for row in rd], dtype=np.float64)

with open(D + "/pData_IMvigor210.csv", encoding="utf-8", errors="replace") as f:
    pr = list(csv.reader(f))
p_rows = pr[1:]
binaryResponse = [r[2] for r in p_rows]
resp = np.full(len(p_rows), np.nan)
for i, b in enumerate(binaryResponse):
    if b == "CR/PR": resp[i] = 1.0
    elif b == "SD/PD": resp[i] = 0.0

def auroc(s, y):
    m = ~np.isnan(y); s, yy = s[m], y[m]
    n1 = np.sum(yy==1); n0 = np.sum(yy==0)
    if n1==0 or n0==0: return np.nan
    r = np.argsort(s); ranks = np.empty_like(r, dtype=float)
    i=0; n=len(r)
    while i<n:
        j=i
        while j+1<n and s[r[j+1]]==s[r[i]]: j+=1
        ranks[r[i:j+1]] = (i+j)/2.0+1.0
        i=j+1
    R1 = ranks[yy==1].sum(); U1 = R1 - n1*(n1+1)/2.0
    return U1/(n1*n0)

def perm_p(s, y, n_iter=10000):
    m=~np.isnan(y); s,yy=s[m],y[m]
    obs=auroc(s,yy); cnt=0
    for _ in range(n_iter):
        if auroc(s, rng.permutation(yy))>=obs: cnt+=1
    return obs, (cnt+1)/(n_iter+1)

# CPM + log1p
def cpm_log1p(cnt):
    lib = cnt.sum(axis=0, keepdims=True); lib[lib==0]=1.0
    return np.log1p(cnt/lib*1e6)

# upper-quartile (UQ) normalization, log1p (robust alternative)
def uq_log1p(cnt):
    q75 = np.percentile(cnt, 75, axis=0); q75[q75==0]=1.0
    return np.log1p(cnt / q75[None,:] * 1e3)

genes = ["CD8A","HLA-DRA","PDCD1","ITGAE"]
expr_cpm = cpm_log1p(counts)
expr_uq  = uq_log1p(counts)

print("=== single-gene AUROC (CPM log1p) ===")
for g in genes:
    s = expr_cpm[sym_to_idx[g], :]
    au, p = perm_p(s, resp)
    print(f"  {g:8s} AUROC={au:.4f}  perm_p={p:.4f}")

print("\n=== 4-gene frozen score (mean) variants ===")
for label, expr in [("CPM", expr_cpm), ("UQ", expr_uq)]:
    idxs = [sym_to_idx[g] for g in genes]
    sc = expr[idxs, :].mean(axis=0)
    au, p = perm_p(sc, resp)
    print(f"  score({label}) AUROC={au:.4f}  perm_p={p:.4f}")

# leave-one-out: does any single gene drive it?
print("\n=== leave-one-out score AUROC (CPM) ===")
for drop in genes:
    idxs = [sym_to_idx[g] for g in genes if g != drop]
    sc = expr_cpm[idxs, :].mean(axis=0)
    au, p = perm_p(sc, resp)
    print(f"  drop {drop:8s} -> AUROC={au:.4f}  perm_p={p:.4f}")

# simple 2-group t-test of each gene (effect size context)
from scipy import stats as st
print("\n=== per-gene log1p(CPM) responder vs non-responder (t-test) ===")
m = ~np.isnan(resp)
for g in genes:
    s = expr_cpm[sym_to_idx[g], :]
    r = s[m & (resp==1)]; nr = s[m & (resp==0)]
    t, pv = st.ttest_ind(r, nr, equal_var=False)
    print(f"  {g:8s} meanR={r.mean():.3f} meanNR={nr.mean():.3f}  t={t:+.2f}  p={pv:.4f}")
