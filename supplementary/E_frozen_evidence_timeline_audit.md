# 冻结证据时间线审计（Frozen-evidence timeline audit）

> 版本：v1_20260905
> 目的：回应复审 Critical「frozen / locked / prespecified 是否在结果前完成」。
> 结论：本地文件时间戳可证部分规则先于结果；但**无公开预注册、无不可变 release、无哈希链**，故全文不得使用 frozen/locked/prespecified 的强含义，统一降级为 analysis-defined 或 retrospectively-locked-before-stated-extension。

## 一、核心判定

| 断言 | 证据 | 判定 |
|---|---|---|
| LOPO 外部分析规则先于结果 | `P1_EXTERNAL_VALIDATION_CONTRACT_20260902.md`（09-02 07:57）早于六队列报告（09-02 19:22）与 `FROZEN_LOPO_EXTERNAL_VALIDATION_BY_COHORT.csv`（09-03 18:52） | 本地时间戳可证，但仅本地版本化文件 |
| 四基因 score 定义先于后续队列验证 | `build_gse78220_score_table.py`（09-01 18:12）含 CD8A/HLA-DRA/PDCD1/ITGAE 定义，早于 IMvigor210（09-04 12:58）、Gide（09-04 13:10）、Braun（09-04 13:32） | 可证"score 定义"先于三大队列，但原始三 ICI 队列与 score 定义几乎同步 |
| 43 靶标 locked panel 先于结果 | `P3_SHARED43_LOCKED_PANEL.csv`（09-02 22:10）早于 `P3_EXPANDED43_TRANSPORT_METRICS.csv`（09-02 22:17） | 仅 7 分钟间隔，接近同时，不能声称结果前预注册 |
| 阈值先于结果 | 阈值敏感性分析 09-04 17:13 才产出 | **阈值本身未在结果前锁定**，属于事后敏感性验证 |
| 公开预注册 / 不可变 release / 哈希链 | 无 | **不存在** |

## 二、时间线（关键里程碑，按文件 mtime 升序）

| 时间 | 事件 | 文件 |
|---|---|---|
| 09-01 18:12 | 四基因 score 定义最早落盘（CD8A/HLA-DRA/PDCD1/ITGAE） | `clinical_gse206325_exploratory/v3_robustness/build_gse78220_score_table.py` |
| 09-02 07:57 | LOPO 外部分析冻结契约 | `p1_external_lopo_validation/P1_EXTERNAL_VALIDATION_CONTRACT_20260902.md` |
| 09-02 19:22 | 六队列外部 transport 报告（结果） | `p1_external_lopo_validation/P1_SIX_COHORT_FINAL_REPORT_20260902.md` |
| 09-02 22:10 | 43 靶标 locked panel | `p3_target_scope_20260902/shared43_datasets_v1/P3_SHARED43_LOCKED_PANEL.csv` |
| 09-02 22:17 | 43 靶标 transport 结果 | `p3_transport_expanded43_v1/P3_EXPANDED43_TRANSPORT_METRICS.csv` |
| 09-03 18:52 | 六队列冻结预测结果 | `external_lopo_predictor_validation_20260903/FROZEN_LOPO_EXTERNAL_VALIDATION_BY_COHORT.csv` |
| 09-03 18:56 | 临床四基因 score 冻结验证（GSE93157） | `clinical_external_validation_v4_20260903/run_gse93157_frozen_validation.py` |
| 09-04 12:58 | IMvigor210 冻结 score 验证 | `imvigor210/validate_frozen_score.py` |
| 09-04 13:10 | Gide 2019 冻结 score 验证 | `gide2019/validate_gide_frozen_score.py` |
| 09-04 13:32 | Braun 2020 冻结 score 验证 | `braun2020/validate_braun_frozen_score.py` |
| 09-04 17:13 | 阈值敏感性分析（事后） | `GSE314416/threshold_sensitivity_v5.py` |

## 三、措辞替换规范（全文统一执行）

| 原文措辞（应弃用） | 替换为 |
|---|---|
| frozen analysis / frozen score / frozen predictor | analysis-defined analysis / score locked before the stated extension |
| prespecified / prospectively locked | retrospectively locked before the stated extension（并注明本地时间戳依据） |
| locked thresholds（无先验依据时） | investigator-defined thresholds with sensitivity analysis（阈值本身事后验证，不称预注册） |
| universal failure / hard failure（绝对含义） | failure across all evaluated settings / observed hard failure |

## 四、可验证的"先于结果"证据（可写入正文/补充）

1. **LOPO 契约**：`P1_EXTERNAL_VALIDATION_CONTRACT_20260902.md` 明确"without re-selecting features or tuning thresholds"、"ridge alpha=1 source-only standardization"，其 mtime 早于六队列结果约 12 小时——可称"analysis rules were fixed in a local versioned contract before external cohort results were generated"。
2. **四基因 score 定义**：在 Gide/IMvigor210/Braun 三大队列验证之前已存在（09-01），可称"the four-gene score was defined before the three largest independent cohorts were evaluated"。
3. **不能称**：任何阈值、任何 43 靶标分类是"结果前预注册"——因为 43 panel 与结果仅隔 7 分钟，阈值敏感性是事后补充。

## 五、若需升级为"可验证冻结"的选项（供作者决策）

1. 将关键规则文件 + SHA-256 上传至 Zenodo/OSF 生成不可变 DOI release，并把 DOI 写入正文 Data availability。
2. 或保持现状，全文降级措辞为"analysis-defined / retrospectively locked"，并附本时间线审计表作为 Supplement。

> 本审计由本地文件 mtime 生成；mtime 可被修改，因此仅作为"作者内部时间线"证据，不构成第三方可验证的预注册。
