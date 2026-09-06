#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE245108 TNC 滴定数据 step4 (v3): 全 265 靶标方差分解 — 增量写入 + 断点续跑
"""
import numpy as np
import json
import os
import pandas as pd
import statsmodels.formula.api as smf
import warnings
warnings.filterwarnings('ignore')

OUT_DIR = r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\technical_biological_decomposition_data_v1_20260903\GSE245108_tnc_parsed"
CSV_PATH = os.path.join(OUT_DIR, 'tnc_variance_components_all.csv')

COLUMNS = ['feature','donor','cell_type','concentration','rep','residual',
           'frac_donor','frac_cell_type','frac_concentration','frac_rep','frac_residual',
           'converged','error']

def main():
    clr = np.load(os.path.join(OUT_DIR, 'tnc_adt_clr.npz'))['clr']
    meta = json.load(open(os.path.join(OUT_DIR, 'tnc_meta_full.json')))
    feat_names = json.load(open(os.path.join(OUT_DIR, 'tnc_feature_names.json')))
    n_cell, n_feat = clr.shape
    print('cells', n_cell, 'features', n_feat, flush=True)

    df = pd.DataFrame({
        'donor': [m['donor'] for m in meta],
        'cell_type': [m['cluster'] for m in meta],
        'concentration': [m['concentration'] for m in meta],
        'rep': [m['rep'] for m in meta],
    })

    # 断点续跑: 读取已完成 feature
    done = set()
    if os.path.exists(CSV_PATH):
        existing = pd.read_csv(CSV_PATH)
        done = set(existing['feature'])
        print(f'resume: {len(done)} already done', flush=True)

    # 打开追加写入
    header_needed = not os.path.exists(CSV_PATH)
    import csv
    fout = open(CSV_PATH, 'a', newline='', encoding='utf-8')
    writer = csv.DictWriter(fout, fieldnames=COLUMNS)
    if header_needed:
        writer.writeheader()

    for i, feat in enumerate(feat_names):
        if feat in done:
            continue
        d = df.copy()
        d['y'] = clr[:, i].astype(np.float64)
        row = {'feature': feat}
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
            total = sum(var.values())
            row.update(var)
            row.update({f'frac_{k}': (v/total if total > 0 else 0) for k, v in var.items()})
            row['converged'] = bool(f.converged)
            row['error'] = ''
        except Exception as e:
            row.update({'donor': np.nan, 'cell_type': np.nan, 'concentration': np.nan,
                        'rep': np.nan, 'residual': np.nan,
                        'frac_donor': np.nan, 'frac_cell_type': np.nan,
                        'frac_concentration': np.nan, 'frac_rep': np.nan, 'frac_residual': np.nan,
                        'converged': False, 'error': str(e)[:200]})
        writer.writerow(row)
        fout.flush()
        if (i+1) % 10 == 0:
            print(f'  {i+1}/{n_feat} done (csv {len(done)+i+1-len(done)} new)', flush=True)

    fout.close()
    print('\nDONE. saved tnc_variance_components_all.csv', flush=True)

    # 汇总
    out = pd.read_csv(CSV_PATH)
    ok = out[out['error'] == '']
    print(f'\nconverged+ok: {len(ok)}/{len(out)}')
    print('frac_concentration 中位数:', round(ok['frac_concentration'].median(), 4))
    print('frac_rep 中位数:', round(ok['frac_rep'].median(), 4))
    print('frac_donor 中位数:', round(ok['frac_donor'].median(), 4))

if __name__ == '__main__':
    main()
