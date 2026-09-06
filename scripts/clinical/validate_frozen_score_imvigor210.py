# -*- coding: utf-8 -*-
"""
IMvigor210 (Mariathasan 2018, Nature 554:544) frozen four-gene score re-validation.

Purpose:
  Re-validate the FROZEN four-gene RNA score against ICI response in a large
  independent cohort (n=348 bladder cancer, anti-PD-L1 atezolizumab), WITHOUT
  any re-training or threshold optimization (pre-registered style).

Frozen score definition (unchanged from prior P5 boundary):
  score = mean( log1p( normalized expression ) ) over {CD8A, HLA-DRA, PDCD1, ITGAE}

Normalization (pre-registered, two variants for sensitivity, NOT optimized on outcome):
  1. CPM:  counts / library_size * 1e6, then log1p   (primary)
  2. FPKM: counts / (len_kb * library_size_millions), then log1p  (sensitivity; uses gene length)

Outcome (binaryResponse, RECIST-derived, no re-derivation):
  responder = CR/PR, non-responder = SD/PD; NA (NE) samples excluded from AUROC.

Metrics (pre-registered):
  AUROC (Mann-Whitney), AUPRC, Brier score (score standardized to [0,1] via logistic
  calibration intercept/slope), permutation p (response shuffled, 10000 iters).
"""
import csv, sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

D = r"D:/ai_multimodal_cholesterol_study_outputs/imvigor210/csv"
OUT = r"D:/ai_multimodal_cholesterol_study_outputs/imvigor210"
rng = np.random.default_rng(20260904)

TARGETS = {"CD8A": 925, "HLA-DRA": 3122, "PDCD1": 5133, "ITGAE": 3682}  # symbol -> entrez

# ---- load fData: symbol -> (row_index 0-based, length) ----
with open(D + "/fData_IMvigor210.csv", encoding="utf-8", errors="replace") as f:
    fr = list(csv.reader(f))
f_hdr = fr[0]
f_rows = fr[1:]  # 31085 rows, order == exmat gene order
sym_to_len = {}
for r in f_rows:
    sym_to_len[r[1]] = float(r[3]) if r[3] else 0.0

# locate four-gene row indices (0-based) in exmat gene order
gene_rows = {}
for sym, ent in TARGETS.items():
    idx = None
    for i, r in enumerate(f_rows):
        if r[1] == sym:
            idx = i
            break
    assert idx is not None, sym
    gene_rows[sym] = idx
    print(f"gene {sym:8s} entrez={ent}  row={idx}  len={sym_to_len[sym]:.0f}")

# ---- load counts (31085 x 348) ----
print("loading counts ...")
counts = []
sample_ids = None
with open(D + "/exmat_censored_IMvigor210.csv", encoding="utf-8", errors="replace") as f:
    reader = csv.reader(f)
    header = next(reader)
    sample_ids = [c.strip('"') for c in header[1:]]
    for row in reader:
        counts.append([float(x) for x in row[1:]])
counts = np.array(counts, dtype=np.float64)
print("counts shape:", counts.shape, "samples:", len(sample_ids))

# ---- load pData response ----
with open(D + "/pData_IMvigor210.csv", encoding="utf-8", errors="replace") as f:
    pr = list(csv.reader(f))
p_hdr = pr[0]
p_rows = pr[1:]
p_sids = [r[0] for r in p_rows]
binaryResponse = [r[2] for r in p_rows]
assert p_sids == sample_ids, "sample order mismatch"

# build response vector: 1=CR/PR, 0=SD/PD, NaN=NA (excluded)
resp = np.full(len(p_rows), np.nan)
for i, b in enumerate(binaryResponse):
    if b == "CR/PR":
        resp[i] = 1.0
    elif b == "SD/PD":
        resp[i] = 0.0
# NA stays nan
n_total = len(resp)
n_R = int(np.nansum(resp))
n_NR = int(np.sum(resp == 0))
n_NA = int(np.sum(np.isnan(resp)))
print(f"n_total={n_total}  responders={n_R}  nonresponders={n_NR}  NA={n_NA}")

# ---- normalization helpers ----
def norm_cpm(cnt):
    lib = cnt.sum(axis=0, keepdims=True)
    lib[lib == 0] = 1.0
    return np.log1p(cnt / lib * 1e6)

def norm_fpkm(cnt):
    # len_kb per gene, library in millions
    lengths = np.array([sym_to_len.get(r[1], 1.0) for r in f_rows], dtype=np.float64)
    len_kb = lengths / 1000.0
    lib = cnt.sum(axis=0, keepdims=True)
    lib[lib == 0] = 1.0
    lib_m = lib / 1e6
    fpkm = cnt / (len_kb[:, None] * lib_m)
    return np.log1p(fpkm)

# ---- score builder ----
def frozen_score(expr_matrix):
    # expr_matrix: n_genes x n_samples, normalized
    idxs = [gene_rows[s] for s in ["CD8A", "HLA-DRA", "PDCD1", "ITGAE"]]
    sub = expr_matrix[idxs, :]  # 4 x n_samples
    return sub.mean(axis=0)

# ---- metrics ----
def auroc(scores, y):
    m = ~np.isnan(y)
    s, yy = scores[m], y[m]
    n1 = np.sum(yy == 1); n0 = np.sum(yy == 0)
    if n1 == 0 or n0 == 0:
        return np.nan
    # Mann-Whitney U
    r = np.argsort(s)
    # rank (1-based) with average for ties
    ranks = np.empty_like(r, dtype=np.float64)
    i = 0
    n = len(r)
    while i < n:
        j = i
        while j + 1 < n and s[r[j + 1]] == s[r[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        ranks[r[i:j + 1]] = avg
        i = j + 1
    R1 = ranks[yy == 1].sum()
    U1 = R1 - n1 * (n1 + 1) / 2.0
    return U1 / (n1 * n0)

def auprc(scores, y):
    m = ~np.isnan(y)
    s, yy = scores[m], y[m]
    order = np.argsort(-s)  # descending
    yy = yy[order]
    tp = np.cumsum(yy == 1).astype(np.float64)
    fp = np.cumsum(yy == 0).astype(np.float64)
    n_pos = np.sum(yy == 1)
    if n_pos == 0:
        return np.nan
    recall = tp / n_pos
    precision = tp / np.maximum(tp + fp, 1e-12)
    # integrate via trapezoid over recall, monotonic precision handled by taking max suffix
    prec = np.maximum.accumulate(precision[::-1])[::-1]
    # prepend (0, prec[0])
    rec = np.concatenate(([0.0], recall))
    prec = np.concatenate(([prec[0]], prec))
    return np.trapezoid(prec, rec)

def brier(scores, y):
    m = ~np.isnan(y)
    s, yy = scores[m], y[m]
    # logistic calibration on standardized score -> prob
    z = (s - s.mean()) / (s.std() + 1e-12)
    # fit intercept + slope via simple logistic (Newton)
    X = np.column_stack([np.ones_like(z), z])
    w = np.zeros(2)
    for _ in range(100):
        eta = X @ w
        p = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
        grad = X.T @ (yy - p)
        W = p * (1 - p)
        H = (X * W[:, None]).T @ X
        try:
            step = np.linalg.solve(H + 1e-6 * np.eye(2), grad)
        except np.linalg.LinAlgError:
            break
        w += step
        if np.abs(step).max() < 1e-8:
            break
    p = 1.0 / (1.0 + np.exp(-np.clip(X @ w, -30, 30)))
    brier = np.mean((yy - p) ** 2)
    return brier, w[0], w[1]  # intercept, slope

def permutation_p(scores, y, n_iter=10000):
    m = ~np.isnan(y)
    s, yy = scores[m], y[m]
    obs = auroc(s, yy)
    cnt = 0
    for _ in range(n_iter):
        yp = rng.permutation(yy)
        a = auroc(s, yp)
        if a >= obs:
            cnt += 1
    p = (cnt + 1) / (n_iter + 1)
    return obs, p

# ---- run ----
results = {}
for label, expr in [("CPM", norm_cpm(counts)), ("FPKM", norm_fpkm(counts))]:
    sc = frozen_score(expr)
    au = auroc(sc, resp)
    ap = auprc(sc, resp)
    br, icpt, slope = brier(sc, resp)
    obs, p = permutation_p(sc, resp, n_iter=10000)
    results[label] = dict(auroc=au, auprc=ap, brier=br, intercept=icpt,
                          slope=slope, perm_p=p, n_eval=int(n_R + n_NR))
    print(f"\n=== {label} normalization ===")
    print(f"  n_evaluable = {int(n_R+n_NR)} (R={n_R}, NR={n_NR}, NA={n_NA})")
    print(f"  AUROC  = {au:.4f}")
    print(f"  AUPRC  = {ap:.4f}")
    print(f"  Brier  = {br:.4f}  (calib intercept={icpt:.3f}, slope={slope:.3f})")
    print(f"  permutation p (10000) = {p:.4f}")

# ---- save scores + response for record ----
sc_cpm = frozen_score(norm_cpm(counts))
with open(OUT + "/imvigor210_frozen_score_with_response.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["sample_id", "binaryResponse", "BOR", "frozen_four_gene_score_cpm"])
    for i in range(len(sample_ids)):
        w.writerow([sample_ids[i], binaryResponse[i], p_rows[i][1], f"{sc_cpm[i]:.8f}"])

# summary json
import json
summary = {
    "cohort": "IMvigor210 (Mariathasan 2018, Nature 554:544)",
    "trial": "NCT02108652/NCT02951767, atezolizumab (anti-PD-L1), metastatic urothelial carcinoma",
    "n_total": n_total, "n_responder_CRPR": n_R, "n_nonresponder_SDPD": n_NR, "n_NA": n_NA,
    "frozen_score": "mean(log1p(norm)) over {CD8A, HLA-DRA, PDCD1, ITGAE}",
    "genes": {"CD8A": 925, "HLA-DRA": 3122, "PDCD1": 5133, "ITGAE": 3682},
    "results": results,
    "notes": "No re-training. No threshold optimization. NA(NE) excluded. Score is continuous; AUROC rank-based.",
}
with open(OUT + "/imvigor210_frozen_score_summary.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)
print("\nSaved: imvigor210_frozen_score_with_response.csv, imvigor210_frozen_score_summary.json")
