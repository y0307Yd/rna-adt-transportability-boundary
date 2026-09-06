# -*- coding: utf-8 -*-
"""
真实数据 PCA 特征值谱 → 估计有效独立靶标数 K（把发现 B1 的 K 从假设升级为实测）

数据：GSE245108 31851 细胞 × 265 ADT 靶标 CLR 矩阵。

方法：
  1. 相关矩阵（265×265）特征值分解。
  2. Neff = (Σλ)²/Σλ²  （有效独立维度）
  3. 平行分析（parallel analysis）：随机打乱每靶标（破坏靶标间相关但保留边缘分布），
     重复 B 次取 95% 分位特征值谱，真实特征值超过该谱的个数 = 显著因子数 K。
  4. 两条谱：
     (a) 全数据（含细胞类型维度）
     (b) 细胞类型残差（先回归掉 78 种 cell_type，再看类型内因子结构）

输出：
  identifiability_sim_v2_20260904/real_pca_spectrum.json
  identifiability_sim_v2_20260904/REAL_PCA_SPECTRUM_REPORT_20260904.md
"""
import os, json, sys, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings('ignore')
sys.stdout.reconfigure(encoding="utf-8")

BASE = r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\technical_biological_decomposition_data_v1_20260903"
OUT = r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\identifiability_sim_v2_20260904"
os.makedirs(OUT, exist_ok=True)

D = os.path.join(BASE, "GSE245108_tnc_parsed")
z = np.load(os.path.join(D, "tnc_adt_clr.npz"))
CLR = z["clr"].astype(np.float64)   # (31851, 265)
z.close()
names = json.load(open(os.path.join(D, "tnc_feature_names.json"), encoding="utf-8"))

meta = json.load(open(os.path.join(D, "tnc_meta_full.json"), encoding="utf-8"))
# meta 是 list，每元素含 cluster 字段（真实细胞类型注释），长度=CLR 行数
print("meta 类型:", type(meta), "长度:", len(meta) if isinstance(meta, list) else '?', flush=True)

N, P = CLR.shape
print(f"\nCLR: {N} cells x {P} targets", flush=True)

def spectrum(X):
    """相关矩阵特征值谱（降序）。"""
    Xc = X - X.mean(axis=0, keepdims=True)
    sd = Xc.std(axis=0, keepdims=True); sd[sd == 0] = 1.0
    Z = Xc / sd
    C = (Z.T @ Z) / (Z.shape[0] - 1)   # P x P 相关矩阵
    eig = np.linalg.eigvalsh(C)[::-1]  # 降序
    return eig

def parallel_analysis(X, B=20, seed=42):
    """平行分析：随机置换每列（破坏列间相关），取 95% 分位特征值谱。"""
    rng = np.random.default_rng(seed)
    Xc = X - X.mean(axis=0, keepdims=True)
    sd = Xc.std(axis=0, keepdims=True); sd[sd == 0] = 1.0
    Z = Xc / sd
    P = Z.shape[1]; N = Z.shape[0]
    spectra = np.zeros((B, P))
    for b in range(B):
        Zp = Z.copy()
        for j in range(P):
            Zp[:, j] = rng.permutation(Zp[:, j])
        C = (Zp.T @ Zp) / (N - 1)
        spectra[b] = np.linalg.eigvalsh(C)[::-1]
    return np.percentile(spectra, 95, axis=0)

def analyze(X, label, do_parallel=True):
    eig = spectrum(X)
    total = eig.sum()
    Neff = float((eig.sum() ** 2) / (eig ** 2).sum()) if eig.sum() > 0 else 1.0
    # 累计解释方差
    cum = np.cumsum(eig) / total
    k50 = int(np.searchsorted(cum, 0.50) + 1)
    k80 = int(np.searchsorted(cum, 0.80) + 1)
    k90 = int(np.searchsorted(cum, 0.90) + 1)
    kaiser = int((eig > eig.mean()).sum())   # Kaiser 准则（特征值>平均）
    print(f"\n=== {label} ===", flush=True)
    print(f"  P={P}, Neff={Neff:.2f}, 独立性比例={Neff/P:.3f}", flush=True)
    print(f"  累计50%需要 {k50} 成分, 80%需 {k80}, 90%需 {k90}", flush=True)
    print(f"  Kaiser(>mean) 成分数: {kaiser}", flush=True)
    print(f"  前12特征值: {np.round(eig[:12], 3)}", flush=True)
    K_par = None
    if do_parallel:
        pa = parallel_analysis(X)
        K_par = int((eig > pa).sum())
        print(f"  平行分析(95%) 显著因子数 K={K_par}", flush=True)
    return dict(label=label, P=P, Neff=round(Neff,2),
                independence_fraction=round(Neff/P,4),
                k50=k50, k80=k80, k90=k90, kaiser=kaiser, K_parallel=K_par,
                top_eigenvalues=[round(float(x),3) for x in eig[:12]])

results = {}

# (a) 全数据谱
results["full"] = analyze(CLR, "全数据（含细胞类型维度）")

# (b) 细胞类型残差谱 —— 从 tnc_meta_full.json 取 cluster
cell_type = None
if isinstance(meta, list) and len(meta) == N:
    cell_type = np.array([m.get("cluster", "NA") for m in meta])
    print(f"\n[cell_type] cluster n_unique={len(np.unique(cell_type))}", flush=True)

if cell_type is not None:
    # 回归掉 cell_type（one-hot 最小二乘残差）
    uniq = np.unique(cell_type)
    Dmat = np.zeros((N, len(uniq)))
    for i, u in enumerate(uniq):
        Dmat[cell_type == u, i] = 1.0
    # 最小二乘投影残差
    beta = np.linalg.lstsq(Dmat, CLR, rcond=None)[0]
    resid = CLR - Dmat @ beta
    results["celltype_resid"] = analyze(resid, "细胞类型残差（类型内因子结构）")
else:
    print("\n[警告] 未找到 cell_type 标签，跳过残差谱", flush=True)

with open(os.path.join(OUT, "real_pca_spectrum.json"), "w", encoding="utf-8") as f:
    json.dump(results, f, indent=2, ensure_ascii=False)

print("\n\n=== 汇总 ===", flush=True)
for k, v in results.items():
    print(f"  {k}: Neff={v['Neff']}  K_parallel={v['K_parallel']}  kaiser={v['kaiser']}  k90={v['k90']}", flush=True)
print("\nOutputs ->", OUT, flush=True)
