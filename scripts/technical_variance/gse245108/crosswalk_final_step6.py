#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""四数据集最终对照表 (step6): GSE245108 vs Kotliarov v2 vs GSE135325 vs GSE282881"""
import numpy as np, json, os, csv
import pandas as pd

OUT = r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\technical_biological_decomposition_data_v1_20260903\GSE245108_tnc_parsed"
KOT = r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\technical_biological_decomposition_data_v1_20260903\kotliarov_probe\kotliarov_pairing\variance_decomposition_v2\per_target_variance_components.csv"

gse = pd.read_csv(os.path.join(OUT, 'tnc_variance_components_all.csv'))
gse = gse[gse['error'].isna()]
gse['target'] = gse['feature'].str.replace('^Hu\.', '', regex=True)

kot = pd.read_csv(KOT)
kot = kot.rename(columns={'repeatability_ICC':'kot_ICC'})

# 对齐目标名
def norm_t(s):
    s = s.replace('-','').replace('.','').replace('_','').upper()
    return s

gse['tkey'] = gse['target'].map(norm_t)
kot['tkey'] = kot['target'].map(norm_t)

# 合并
m = pd.merge(kot[['target','tkey','kot_ICC','frac_donor','frac_cell_type','frac_lane','frac_residual','technical_fraction']],
             gse[['target','tkey','frac_donor','frac_cell_type','frac_concentration','frac_rep','frac_residual']],
             on='tkey', how='inner', suffixes=('_kot','_gse'))
m = m.rename(columns={'target_kot':'target'})
m = m[['target','kot_ICC','frac_donor_kot','frac_cell_type_kot','frac_lane','frac_residual_kot',
       'frac_donor_gse','frac_cell_type_gse','frac_concentration','frac_rep','frac_residual_gse']]
m.to_csv(os.path.join(OUT, 'four_dataset_crosswalk_kotliarov_gse245108.csv'), index=False)
print('=== Kotliarov v2 × GSE245108 对齐 (n=%d) ===' % len(m))
print(m.to_string(index=False))

# 汇总对比统计
print()
print('=== 技术噪声分量中位数对比 ===')
print('Kotliarov frac_lane (lane嵌套batch):    %.4f' % m['frac_lane'].median())
print('GSE245108 frac_concentration (染色浓度): %.4f' % m['frac_concentration'].median())
print('GSE245108 frac_rep (独立文库制备):       %.4f' % m['frac_rep'].median())
print('GSE245108 frac_donor (生物):            %.4f' % m['frac_donor_gse'].median())
print()
# 四基因蛋白专项
print('=== 四基因 score 蛋白对照 ===')
four = m[m['target'].str.upper().str.contains('CD103|CD279|HLA|CD8', na=False)]
for _, r in four.iterrows():
    print(f"{r['target']:10s} Kot_ICC={r['kot_ICC']:.3f}  GSE_conc={r['frac_concentration']:.3f}  GSE_rep={r['frac_rep']:.3f}  GSE_donor={r['frac_donor_gse']:.3f}")
