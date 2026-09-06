#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
稳健性/混杂检查：Gide 2019 上四基因 score 显著是否只是"免疫浸润"混杂
1) 对照基因集 AUROC 背景分布（免疫浸润标志物 vs 随机基因）
2) 原始 counts 重算（归一化方法敏感性）
3) leave-one-out
4) 文库大小/总表达混杂检查
"""
import sys, csv, json
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
d = r"D:/ai_multimodal_cholesterol_study_outputs/gide2019"
from scipy.stats import mannwhitneyu, ttest_ind

GENES = ["CD8A", "HLA-DRA", "PDCD1", "ITGAE"]
# immune infiltration / T-cell / cytotoxic marker set (positive controls)
IMMUNE = ["CD3D","CD3E","CD2","PTPRC","GZMB","PRF1","CD8B","CXCR3","CD27","LAG3",
          "CTLA4","TIGIT","HAVCR2","IFNG","CXCL9","CXCL10","CCL5","GZMA","NKG7","CD96"]
# housekeeping / random background
HOUSEKEEPING = ["ACTB","GAPDH","TUBB","B2M","RPL13A","RPLP0","HPRT1","PGK1","VIM","YWHAZ"]

with open(d + "/cancercell_normalized_counts_genenames.txt", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f, delimiter="\t")
    header = next(rd)
    all_rows = list(rd)
samples = header[3:]
gene_symbol = {r[1]: i for i, r in enumerate(all_rows)}

# response annotation
import pyreadr
coldata = pyreadr.read_r(d + "/coldata_example.rda")["coldata_example"]
def rowname_to_col(name):
    parts = name.split("_")
    num = parts[1]
    trt = "PD1" if parts[2] == "PD1" else "ipiPD1"
    time = parts[3]
    return f"{trt}_{num}_{time}"
sample_to_col = {s: i for i, s in enumerate(samples)}
meta = []
for rn, row in coldata.iterrows():
    col = rowname_to_col(rn)
    if col in sample_to_col:
        meta.append((col, row["Response"], ("PD1" if col.startswith("PD1_") else "ipiPD1")))
y = np.array([1 if r == "R" else 0 for _, r, _ in meta])

def get_vals(g):
    if g not in gene_symbol:
        return None
    r = all_rows[gene_symbol[g]]
    return np.array([float(x) for x in r[3:]])

def log1p_of_col(colname_list):
    idx = [sample_to_col[c] for c, _, _ in meta]
    out = {}
    for c, r, t in meta:
        pass

# build sample-aligned expression for all genes of interest
sample_idx = [sample_to_col[c] for c, _, _ in meta]
def gene_matrix(gs):
    M = []
    for g in gs:
        v = get_vals(g)
        M.append(np.log1p(v[sample_idx]) if v is not None else None)
    return gs, M

def auc(y, s):
    if len(np.unique(y)) < 2: return np.nan
    pos = s[y==1]; neg = s[y==0]
    if len(pos)==0 or len(neg)==0: return np.nan
    u,_ = mannwhitneyu(pos, neg, alternative="two-sided")
    return u/(len(pos)*len(neg))

def perm_p(y, s, nperm=2000, seed=1):
    rng = np.random.default_rng(seed)
    obs = auc(y, s)
    cnt = 0
    for _ in range(nperm):
        if auc(rng.permutation(y), s) >= obs: cnt += 1
    return obs, (cnt+1)/(nperm+1)

print("=== 1) 对照基因集 AUROC 背景分布 (全队列 n=73) ===")
def report_set(name, gs):
    gs2, M = gene_matrix(gs)
    print(f"\n[{name}]")
    for g, v in zip(gs2, M):
        if v is None:
            print(f"  {g}: MISSING"); continue
        a, p = perm_p(y, v, nperm=1000)
        # also t-test direction
        t, pt = ttest_ind(v[y==1], v[y==0])
        print(f"  {g}: AUROC={a:.3f} perm_p={p:.3f} t={t:+.2f} p_t={pt:.4f}")

report_set("免疫浸润/T细胞/细胞毒 (阳性对照)", IMMUNE)
report_set("管家基因 (阴性对照)", HOUSEKEEPING)

print("\n=== 2) 四基因 leave-one-out (全队列) ===")
gs4, M4 = gene_matrix(GENES)
base = np.mean([v for v in M4], axis=0)
a, p = perm_p(y, base, nperm=5000)
print(f"4-gene base: AUROC={a:.4f} perm_p={p:.4f}")
for i, g in enumerate(GENES):
    dropped = np.mean([M4[j] for j in range(4) if j != i], axis=0)
    a2, p2 = perm_p(y, dropped, nperm=5000)
    print(f"  drop {g}: AUROC={a2:.4f} perm_p={p2:.4f}")

print("\n=== 3) 文库大小/总表达混杂检查 ===")
# total expression per sample (sum of all genes, proxy for library size / tumor purity)
all_vals = np.array([[float(x) for x in r[3:]] for r in all_rows])
lib_size = all_vals.sum(axis=0)[sample_idx]
lib_log = np.log1p(lib_size)
a_lib, p_lib = perm_p(y, lib_log, nperm=2000)
print(f"library size log1p: AUROC={a_lib:.4f} perm_p={p_lib:.4f}")
# correlation between score and library size
from scipy.stats import spearmanr, pearsonr
rho, pr = spearmanr(base, lib_log)
print(f"Spearman(score, lib_size) = {rho:.3f} (p={pr:.2e})")
rho2, pr2 = spearmanr(np.log1p(M4[1]), lib_log)  # HLA-DRA
print(f"Spearman(HLA-DRA, lib_size) = {rho2:.3f} (p={pr2:.2e})")

print("\n=== 4) 原始 counts 重算 (PD1-IPIPD1_counts.txt) ===")
with open(d + "/PD1-IPIPD1_counts.txt", encoding="utf-8", errors="replace") as f:
    rd2 = csv.reader(f, delimiter="\t")
    hdr2 = next(rd2)
    raw_rows = list(rd2)
raw_samples = hdr2[1:]  # first col = ID
raw_gene = {r[0]: i for i, r in enumerate(raw_rows)}  # key by ENSG ID
# need gene symbol mapping from cancercell file (col0=ID, col1=Gene)
id_to_sym = {r[0]: r[1] for r in all_rows}
raw_sym_to_id = {}
for r in raw_rows:
    eid = r[0]
    sym = id_to_sym.get(eid, None)
    if sym: raw_sym_to_id[sym] = eid
raw_sample_to_col = {s: i for i, s in enumerate(raw_samples)}
print("raw samples n:", len(raw_samples), "raw genes n:", len(raw_rows))
# check overlap
raw_meta = [(c, r, t) for c, r, t in meta if c in raw_sample_to_col]
print("raw overlap samples:", len(raw_meta))

def raw_vals(g):
    if g not in raw_sym_to_id: return None
    eid = raw_sym_to_id[g]
    if eid not in raw_gene: return None
    r = raw_rows[raw_gene[eid]]
    return np.array([float(x) for x in r[1:]])

ry = np.array([1 if r == "R" else 0 for _, r, _ in raw_meta])
ridx = [raw_sample_to_col[c] for c, _, _ in raw_meta]
raw_score = None
raw_hla = None
for g in GENES:
    v = raw_vals(g)
    if v is None:
        print(f"  raw {g}: MISSING"); continue
    lv = np.log1p(v[ridx])
    a3, p3 = perm_p(ry, lv, nperm=1000)
    print(f"  raw {g} log1p: AUROC={a3:.4f} perm_p={p3:.4f}")
    if g == "HLA-DRA": raw_hla = lv
    if raw_score is None:
        raw_score = lv.copy()
    else:
        raw_score += lv
raw_score /= 4
a4, p4 = perm_p(ry, raw_score, nperm=5000)
print(f"  raw 4-gene mean log1p counts: AUROC={a4:.4f} perm_p={p4:.4f}")
if raw_hla is not None:
    a5, p5 = perm_p(ry, raw_hla, nperm=5000)
    print(f"  raw HLA-DRA log1p: AUROC={a5:.4f} perm_p={p5:.4f}")
