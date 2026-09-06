#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE245108 TNC 滴定数据 step2: 合并细胞类型注释 + 归一化 + 构建方差分解输入

目标：
  1. 把 Level3R-Titrated-Titration-Annotation.txt 的 Cluster 按 barcode 前缀合并进 meta
  2. ADT 计数做 CLR (centered log-ratio) 归一化 (log1p 后中心化, 标准 ADT 处理)
  3. 输出最终配对数据: tnc_adt_clr.npz + tnc_meta_full.json

关键设计事实 (已核实)：
  - TNC-1-2-1/2, TNC-3-4-1/2 是滴定数据 (5 浓度 0.25/0.5/1/2/4)
  - 4 donor (WM22/WF32/BF32/BM32) x 5 浓度 x 2 rep 完全交叉
  - Cluster 注释来自 RNA 转录组 (label transfer + scTriangulate)
"""
import pandas as pd
import numpy as np
import json
import os

DATA_DIR = r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\technical_biological_decomposition_data_v1_20260903\GSE245108_meta"
OUT_DIR = r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\technical_biological_decomposition_data_v1_20260903\GSE245108_tnc_parsed"

def main():
    # 加载 step1 产出
    counts = np.load(os.path.join(OUT_DIR, 'tnc_adt_counts.npz'))['counts']
    meta = json.load(open(os.path.join(OUT_DIR, 'tnc_meta.json')))
    feat_names = json.load(open(os.path.join(OUT_DIR, 'tnc_feature_names.json')))

    print('loaded counts', counts.shape, 'meta', len(meta))

    # 加载 annotation, 提取滴定 TNC 的 barcode -> cluster 映射
    ann = pd.read_csv(os.path.join(DATA_DIR, 'GSE245108_Level3R-Titrated-Titration-Annotation.txt'), sep='\t')
    tit_tnc = ann[ann['UID'].str.contains('TNC-', na=False)].copy()
    # UID = "<barcode-1>.TNC-<sample>", 提取 barcode 部分
    tit_tnc['barcode'] = tit_tnc['UID'].str.split('.').str[0]
    barcode2cluster = dict(zip(tit_tnc['barcode'], tit_tnc['Cluster']))
    print('annotation barcodes:', len(barcode2cluster))

    # 合并 cluster 到 meta
    n_matched = 0
    for m in meta:
        c = barcode2cluster.get(m['barcode'])
        m['cluster'] = c if c else 'NA'
        if c:
            n_matched += 1
    print('cluster matched:', n_matched, '/', len(meta))

    # 过滤掉无 cluster 注释的细胞
    keep = [i for i, m in enumerate(meta) if m['cluster'] != 'NA']
    counts = counts[keep]
    meta = [meta[i] for i in keep]
    print('after cluster filter:', counts.shape, 'cells')

    # CLR 归一化: log1p 后每细胞中心化 (除以每细胞的几何均值)
    # 标准做法: clr(x) = log(x+1) - mean(log(x+1))
    logc = np.log1p(counts)
    clr = logc - logc.mean(axis=1, keepdims=True)

    # 保存
    np.savez_compressed(os.path.join(OUT_DIR, 'tnc_adt_clr.npz'), clr=clr.astype(np.float32))
    json.dump(meta, open(os.path.join(OUT_DIR, 'tnc_meta_full.json'), 'w'))
    json.dump(feat_names, open(os.path.join(OUT_DIR, 'tnc_feature_names.json'), 'w'))

    # 统计
    from collections import Counter
    print('\nfinal concentration counts:', dict(Counter(m['concentration'] for m in meta)))
    print('final donor counts:', dict(Counter(m['donor'] for m in meta)))
    print('final rep counts:', dict(Counter(m['rep'] for m in meta)))
    print('n clusters:', len(set(m['cluster'] for m in meta)))
    # 每 cluster 细胞数 (top 15)
    cc = Counter(m['cluster'] for m in meta)
    print('\ntop 15 clusters:')
    for k, v in cc.most_common(15):
        print('  ', k, v)

if __name__ == '__main__':
    main()
