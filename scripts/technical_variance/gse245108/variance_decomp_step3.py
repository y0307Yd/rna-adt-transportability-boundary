#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE245108 TNC 滴定数据 step3 (v2): 分层方差分解 (正确提取方差分量)

模型 (每靶标, ADT CLR 值):
  ADT_CLR ~ 1 + (1|donor) + (1|cell_type) + (1|concentration) + (1|rep)

方差分量提取 (已用 cov_re 精确验证):
  实际方差 = f.params["<组> Var"] * f.scale
    donor Var         = params["Group Var"] * scale
    cell_type Var     = params["cell_type Var"] * scale
    concentration Var = params["concentration Var"] * scale
    rep Var           = params["rep Var"] * scale
    residual          = scale

随机效应含义:
  donor          = 生物学供体方差
  cell_type      = 细胞类型组成方差 (控制项)
  concentration  = 抗体染色浓度噪声 (0.25x~4x, 本数据集新增正交杠杆)
  rep            = 技术重复方差 (rep1/rep2)
"""
import numpy as np
import json
import os
import pandas as pd
import statsmodels.formula.api as smf
import warnings
warnings.filterwarnings('ignore')

OUT_DIR = r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\technical_biological_decomposition_data_v1_20260903\GSE245108_tnc_parsed"

def main():
    clr = np.load(os.path.join(OUT_DIR, 'tnc_adt_clr.npz'))['clr']
    meta = json.load(open(os.path.join(OUT_DIR, 'tnc_meta_full.json')))
    feat_names = json.load(open(os.path.join(OUT_DIR, 'tnc_feature_names.json')))

    n_cell, n_feat = clr.shape
    print('cells', n_cell, 'features', n_feat)

    df = pd.DataFrame({
        'donor': [m['donor'] for m in meta],
        'cell_type': [m['cluster'] for m in meta],
        'concentration': [m['concentration'] for m in meta],
        'rep': [m['rep'] for m in meta],
    })

    key_targets = ['Hu.CD103', 'Hu.CD8', 'Hu.CD279', 'Hu.HLA.DR',
                   'Hu.CD274', 'Hu.CD3_UCHT1', 'Hu.CD4_RPA.T4',
                   'Hu.CD45RA', 'Hu.CD45RO', 'Hu.CD56', 'Hu.CD19',
                   'Hu.CD16', 'Hu.CD45_2D1', 'Hu.CD38_HIT2', 'Hu.CD69']

    results = {}
    for feat in key_targets:
        if feat not in feat_names:
            print('SKIP missing', feat)
            continue
        idx = feat_names.index(feat)
        d = df.copy()
        d['y'] = clr[:, idx].astype(np.float64)
        try:
            m = smf.mixedlm("y ~ 1", d, groups=d['donor'], re_formula="1",
                            vc_formula={'cell_type':'0+C(cell_type)',
                                        'concentration':'0+C(concentration)',
                                        'rep':'0+C(rep)'})
            f = m.fit(reml=True, method='lbfgs', maxiter=300)
            scale = float(f.scale)
            p = f.params
            var = {
                'donor': float(p.get('Group Var', 0)) * scale,
                'cell_type': float(p.get('cell_type Var', 0)) * scale,
                'concentration': float(p.get('concentration Var', 0)) * scale,
                'rep': float(p.get('rep Var', 0)) * scale,
                'residual': scale,
            }
            # 归一化比例
            total = sum(var.values())
            frac = {k: v/total for k, v in var.items()} if total > 0 else {}
            results[feat] = {'var': var, 'frac': frac}
            print(f'OK {feat}: frac donor={frac.get("donor",0):.4f} '
                  f'cell_type={frac.get("cell_type",0):.4f} '
                  f'conc={frac.get("concentration",0):.4f} '
                  f'rep={frac.get("rep",0):.4f} '
                  f'resid={frac.get("residual",0):.4f}')
        except Exception as e:
            print(f'FAIL {feat}: {e}')
            results[feat] = {'error': str(e)}

    json.dump(results, open(os.path.join(OUT_DIR, 'tnc_variance_components_key.json'), 'w'), indent=2)
    print('\nsaved key results to tnc_variance_components_key.json')

if __name__ == '__main__':
    main()
