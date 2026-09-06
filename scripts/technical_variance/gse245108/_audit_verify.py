#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
审计验证：确认 statsmodels MixedLM 方差分量的正确提取公式
对比 f.params["<组> Var"]、f.cov_re 对角线、f.scale 三者的关系
"""
import numpy as np, json, os, sys
import pandas as pd
import statsmodels.formula.api as smf
import warnings
warnings.filterwarnings('ignore')
sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\technical_biological_decomposition_data_v1_20260903\GSE245108_tnc_parsed"

clr = np.load(os.path.join(OUT_DIR, 'tnc_adt_clr.npz'))['clr']
meta = json.load(open(os.path.join(OUT_DIR, 'tnc_meta_full.json')))
feat_names = json.load(open(os.path.join(OUT_DIR, 'tnc_feature_names.json')))
print('cells', clr.shape[0], 'features', clr.shape[1])

df = pd.DataFrame({
    'donor': [m['donor'] for m in meta],
    'cell_type': [m['cluster'] for m in meta],
    'concentration': [m['concentration'] for m in meta],
    'rep': [m['rep'] for m in meta],
})

# 检查 meta 字段类型和取值
print('\nconcentration unique:', sorted(set(df['concentration'])))
print('rep unique:', sorted(set(df['rep'])))
print('donor unique:', sorted(set(df['donor'])))
print('n cell_type:', df['cell_type'].nunique())

target = 'Hu.HLA.DR'
idx = feat_names.index(target)
d = df.copy()
d['y'] = clr[:, idx].astype(np.float64)

m = smf.mixedlm("y ~ 1", d, groups=d['donor'], re_formula="1",
                vc_formula={'cell_type':'0+C(cell_type)',
                            'concentration':'0+C(concentration)',
                            'rep':'0+C(rep)'})
f = m.fit(reml=True, method='lbfgs', maxiter=300)

print('\n===== f.params (完整) =====')
print(f.params)
print('\n===== f.scale =====')
print(f.scale)
print('\n===== f.cov_re (主 groups=donor) =====')
print(f.cov_re)

# vc 的协方差矩阵
print('\n===== f.random_effects (keys) =====')
print(list(f.random_effects.keys()))

# 对比：params 里随机效应参数 vs 平方 vs 乘 scale
print('\n===== 提取公式对比 =====')
p = f.params
for k in ['Group Var', 'cell_type Var', 'concentration Var', 'rep Var']:
    if k in p.index:
        v = float(p[k])
        print(f'{k}: params={v:.6f} | params^2={v**2:.6f} | params*scale={v*f.scale:.6f}')

# cov_re 对角线（如果可访问）
try:
    print('\ncov_re diag:', np.diag(f.cov_re))
except Exception as e:
    print('cov_re diag err:', e)

# vc 的各组方差
try:
    for k, re in f.random_effects.items():
        if hasattr(re, 'values'):
            pass
except Exception as e:
    pass

# 用 f.summary() 看随机效应
print('\n===== f.summary 随机效应部分 =====')
sm = f.summary()
# 打印包含 Var 的行
for line in str(sm).split('\n'):
    if 'Var' in line or 'Group' in line or 'RE' in line or 'vc' in line.lower():
        print(line)
