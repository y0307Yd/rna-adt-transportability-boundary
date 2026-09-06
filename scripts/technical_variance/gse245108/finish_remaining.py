#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""只处理 todo_remaining 的靶标，前台快速收尾"""
import numpy as np, json, os, csv
import pandas as pd
import statsmodels.formula.api as smf
import warnings
warnings.filterwarnings('ignore')

OUT = r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\technical_biological_decomposition_data_v1_20260903\GSE245108_tnc_parsed"
CSV_PATH = os.path.join(OUT, 'tnc_variance_components_all.csv')
COLUMNS = ['feature','donor','cell_type','concentration','rep','residual',
           'frac_donor','frac_cell_type','frac_concentration','frac_rep','frac_residual',
           'converged','error']

clr = np.load(os.path.join(OUT, 'tnc_adt_clr.npz'))['clr']
meta = json.load(open(os.path.join(OUT, 'tnc_meta_full.json')))
feat_names = json.load(open(os.path.join(OUT, 'tnc_feature_names.json')))

df = pd.DataFrame({
    'donor': [m['donor'] for m in meta],
    'cell_type': [m['cluster'] for m in meta],
    'concentration': [m['concentration'] for m in meta],
    'rep': [m['rep'] for m in meta],
})

done = set()
if os.path.exists(CSV_PATH):
    done = set(pd.read_csv(CSV_PATH)['feature'])

todo = [f for f in feat_names if f not in done]
print('todo:', len(todo), flush=True)

fout = open(CSV_PATH, 'a', newline='', encoding='utf-8')
writer = csv.DictWriter(fout, fieldnames=COLUMNS)

for k, feat in enumerate(todo):
    i = feat_names.index(feat)
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
        row.update({f'frac_{kk}': (v/total if total > 0 else 0) for kk, v in var.items()})
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
    print(f'  {k+1}/{len(todo)} {feat} done', flush=True)

fout.close()
print('DONE remaining', len(todo), flush=True)
