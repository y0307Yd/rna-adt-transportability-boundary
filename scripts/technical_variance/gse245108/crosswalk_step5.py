#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
四数据集技术/生物方差分解横向对照 (step5)

整合:
  1. GSE245108 TNC 滴定 (本数据集) — donor/cell_type/concentration/rep 分层
  2. Kotliarov 2020 — donor/cell_type/lane 分层 (donor×batch 混淆)
  3. GSE135325 — 双平台 experiment 效应
  4. GSE282881 — split-pool 技术重复 CV

产出:
  - four_dataset_crosswalk.csv  (按靶标对齐的方差分量对照)
  - 汇总统计
"""
import numpy as np
import json
import os
import pandas as pd

OUT_DIR = r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\technical_biological_decomposition_data_v1_20260903\GSE245108_tnc_parsed"
KOTLIAROV_DIR = r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\technical_biological_decomposition_data_v1_20260903\kotliarov_probe"

def main():
    # 1. GSE245108 全量结果
    all_csv = os.path.join(OUT_DIR, 'tnc_variance_components_all.csv')
    if not os.path.exists(all_csv):
        print('GSE245108 all CSV not ready yet:', all_csv)
        return
    gse = pd.read_csv(all_csv)
    # 规范化 feature 名 (去掉 Hu. 前缀)
    gse['target'] = gse['feature'].str.replace('^Hu\.', '', regex=True)

    # 2. Kotliarov 结果 (若存在)
    kot = None
    kot_csv = os.path.join(KOTLIAROV_DIR, 'kotliarov_pairing', 'variance_decomposition_v2', 'per_target_variance_components.csv')
    if os.path.exists(kot_csv):
        kot = pd.read_csv(kot_csv)
        print('loaded Kotliarov:', kot.shape)
        print('Kotliarov columns:', list(kot.columns))

    # 关键靶标对照
    key = ['CD103', 'CD8', 'CD279', 'HLA-DR', 'CD274', 'CD3', 'CD4',
           'CD45RA', 'CD45RO', 'CD56', 'CD19', 'CD16', 'CD38', 'CD69',
           'CD11b', 'CD11c', 'CD25', 'CD27', 'CD28', 'CD57', 'CD185',
           'CD196', 'CD244', 'CD123', 'CD303', 'CD223', 'CD137', 'CD33',
           'CD141', 'CD127', 'CD194', 'CD314', 'CD161', 'CD21', 'CD1c',
           'CD183', 'CD20']

    print('\n===== GSE245108 关键靶标方差分解 =====')
    rows = []
    for t in key:
        sub = gse[gse['target'] == t]
        if len(sub) == 0:
            # 尝试模糊匹配
            sub = gse[gse['target'].str.contains(t.replace('-',''), na=False, regex=False)]
        if len(sub) == 0:
            continue
        sub = sub.iloc[0]
        rows.append({
            'target': t,
            'frac_donor': round(sub['frac_donor'], 4),
            'frac_cell_type': round(sub['frac_cell_type'], 4),
            'frac_concentration': round(sub['frac_concentration'], 4),
            'frac_rep': round(sub['frac_rep'], 4),
            'frac_residual': round(sub['frac_residual'], 4),
        })
    df_key = pd.DataFrame(rows)
    print(df_key.to_string(index=False))

    # 汇总全靶标 concentration 分布
    ok = gse[gse['error'].isna()]
    print(f'\n===== GSE245108 全 265 靶标 concentration 分布 =====')
    print('n ok:', len(ok), '/', len(gse))
    print(ok['frac_concentration'].describe().round(4).to_string())
    print('\nrep 分布:')
    print(ok['frac_rep'].describe().round(4).to_string())
    print('\ndonor 分布:')
    print(ok['frac_donor'].describe().round(4).to_string())

    # 保存关键靶标对照
    df_key.to_csv(os.path.join(OUT_DIR, 'tnc_key_targets_summary.csv'), index=False)
    print('\nsaved tnc_key_targets_summary.csv')

if __name__ == '__main__':
    main()
