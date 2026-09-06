#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE245108 TNC 滴定数据解析脚本 (step1)
目标：把 4 个 TNC filtered_feature_bc_matrix.h5 解析为统一配对数据：
  - ADT 计数矩阵 (细胞 x 265 抗体 + 10 HTO hashtag)
  - 每细胞的浓度(0.25/0.5/1/1.5/2/4) + donor(WM22/WF32/BF32/BM32) + rep(1/2) 标签
  - 细胞类型注释(来自 Level3R-Titrated-Titration-Annotation.txt 的 Cluster)
输出：
  - tnc_adt_counts.npz   (ADT counts, cells x features)
  - tnc_meta.json        (cell barcodes + concentration/donor/rep/population labels)
  - tnc_feature_names.json

设计依据(已核实)：
  TNC-1-2-1/2 = donor1+donor2 (WM22+WF32), rep1/rep2
  TNC-3-4-1/2 = donor3+donor4 (BF32+BM32), rep1/rep2
  HTO hashtag 名: TNC_<donoridx>_<conc>  其中 conc in {0_25x,0_5x,1x,2x,4x}
  (注意: 1.5 浓度出现在 annotation 但 HTO 命名未见 1_5x, 需核实)
"""
import h5py
import numpy as np
import json
import os

DATA_DIR = r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\technical_biological_decomposition_data_v1_20260903\GSE245108_titration"
META_DIR = r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\technical_biological_decomposition_data_v1_20260903\GSE245108_meta"
OUT_DIR = r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\technical_biological_decomposition_data_v1_20260903\GSE245108_tnc_parsed"
os.makedirs(OUT_DIR, exist_ok=True)

TNC_FILES = [
    ("GSE245108_TNC-1-2-1_filtered_feature_bc_matrix.h5", "rep1"),
    ("GSE245108_TNC-1-2-2_filtered_feature_bc_matrix.h5", "rep2"),
    ("GSE245108_TNC-3-4-1_filtered_feature_bc_matrix.h5", "rep1"),
    ("GSE245108_TNC-3-4-2_filtered_feature_bc_matrix.h5", "rep2"),
]
# 全局 donor 索引 -> 供体名 (来自 series_matrix: WM22 白男, WF32 白女, BF32 黑女, BM32 黑男)
DONOR_MAP = {'1':'WM22', '2':'WF32', '3':'BF32', '4':'BM32'}

def get_str(g):
    d = g[()]
    return [x.decode() if isinstance(x, bytes) else str(x) for x in d]

def read_csr_matrix(f):
    # 10x h5 用 CSC 存储: shape=(n_feat, n_barcode), indptr 长度 = n_barcode+1
    data = f['matrix/data'][()]
    indices = f['matrix/indices'][()]
    indptr = f['matrix/indptr'][()]
    shape = tuple(f['matrix/shape'][()])
    from scipy.sparse import csc_matrix
    return csc_matrix((data, indices, indptr), shape=shape)

def load_one(path, rep):
    f = h5py.File(path, 'r')
    ftype = get_str(f['matrix/features/feature_type'])
    fname = get_str(f['matrix/features/name'])
    barcodes = get_str(f['matrix/barcodes'])
    n_feat = len(ftype)
    n_cell = len(barcodes)
    mat = read_csr_matrix(f)  # (n_feat x n_cell)
    f.close()

    gene_idx = [i for i, t in enumerate(ftype) if t == 'Gene Expression']
    ab_idx = [i for i, t in enumerate(ftype) if t == 'Antibody Capture']
    ab_names = [fname[i] for i in ab_idx]

    # ADT 抗体 (非 HTO, 非 Isotype, 非 Ig Fc)
    hto_names = [n for n in ab_names if n.startswith('TNC_')]
    isotype_names = [n for n in ab_names if n.startswith('Isotype') or n == 'Hu.IgG.Fc']
    adt_names = [n for n in ab_names if not n.startswith('TNC_') and not n.startswith('Isotype') and n != 'Hu.IgG.Fc']

    hto_idx = [ab_idx[i] for i, n in enumerate(ab_names) if n in hto_names]
    adt_idx = [ab_idx[i] for i, n in enumerate(ab_names) if n in adt_names]

    # 提取矩阵 (转置为 cells x features)
    # mat 是 (n_feat x n_cell), 取子集再转置
    adt_counts = mat[adt_idx, :].T.toarray().astype(np.float32) if len(adt_idx) else np.zeros((n_cell,0))
    hto_counts = mat[hto_idx, :].T.toarray().astype(np.float32) if len(hto_idx) else np.zeros((n_cell,0))

    # HTO 去卷积: 每细胞取 HTO 计数最大的 hashtag 作为标签
    hto_names_local = [ab_names[i] for i, n in enumerate(ab_names) if n in hto_names]
    hto_assign = np.argmax(hto_counts, axis=1) if hto_counts.shape[1] else np.zeros(n_cell, dtype=int)

    meta = []
    for c in range(n_cell):
        hto = hto_names_local[hto_assign[c]] if len(hto_names_local) else 'NA'
        # 解析 TNC_<donor>_<conc>
        conc = 'NA'; donor = 'NA'
        if hto.startswith('TNC_'):
            parts = hto.split('_')  # ['TNC','1','0','25x'] etc
            if len(parts) >= 3:
                don_idx = parts[1]
                conc_part = parts[2]
                if len(parts) >= 4:
                    conc_part = parts[2] + '_' + parts[3]
                conc_map = {'0_25x':'0.25','0_5x':'0.5','1x':'1','2x':'2','4x':'4'}
                conc = conc_map.get(conc_part, conc_part)
                try:
                    donor = DONOR_MAP.get(don_idx, 'NA')
                except Exception:
                    donor = 'NA'
        meta.append({'barcode': barcodes[c], 'rep': rep, 'donor': donor, 'concentration': conc, 'hto': hto})

    return adt_names, adt_counts, meta

def main():
    all_adt_names = None
    all_counts = []
    all_meta = []

    for fname, rep in TNC_FILES:
        p = os.path.join(DATA_DIR, fname)
        exp = {'GSE245108_TNC-1-2-1_filtered_feature_bc_matrix.h5':36400985,
               'GSE245108_TNC-1-2-2_filtered_feature_bc_matrix.h5':35981413,
               'GSE245108_TNC-3-4-1_filtered_feature_bc_matrix.h5':44471480,
               'GSE245108_TNC-3-4-2_filtered_feature_bc_matrix.h5':43014084}[fname]
        if not os.path.exists(p) or os.path.getsize(p) != exp:
            print('SKIP (incomplete/missing)', fname, os.path.getsize(p) if os.path.exists(p) else 'NA')
            continue
        print('loading', fname, '...')
        adt_names, adt_counts, meta = load_one(p, rep)
        if all_adt_names is None:
            all_adt_names = adt_names
        else:
            assert all_adt_names == adt_names, 'ADT name mismatch across files!'
        all_counts.append(adt_counts)
        all_meta.extend(meta)
        print('  cells:', len(meta), '| adt features:', len(adt_names))

    counts = np.vstack(all_counts)
    print('\nTOTAL cells:', counts.shape[0], '| ADT features:', counts.shape[1])

    np.savez_compressed(os.path.join(OUT_DIR, 'tnc_adt_counts.npz'), counts=counts)
    json.dump(all_meta, open(os.path.join(OUT_DIR, 'tnc_meta.json'), 'w'))
    json.dump(all_adt_names, open(os.path.join(OUT_DIR, 'tnc_feature_names.json'), 'w'))

    # 打印标签分布
    from collections import Counter
    print('\nconcentration counts:', dict(Counter(m['concentration'] for m in all_meta)))
    print('donor counts:', dict(Counter(m['donor'] for m in all_meta)))
    print('rep counts:', dict(Counter(m['rep'] for m in all_meta)))

if __name__ == '__main__':
    main()
