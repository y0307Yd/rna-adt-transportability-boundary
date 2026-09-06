from pathlib import Path
import json, pandas as pd

ROAD=Path(r'D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap')
OUT=ROAD/'integrated_evidence_v5_20260903'; OUT.mkdir(exist_ok=True)
MAP=pd.read_csv(ROAD/'p3_target_scope_20260902'/'P3_SHARED_TARGET_RNA_MAPPING_20260902.csv')
f6=ROAD/'p1_external_lopo_validation'/'P1_SIX_COHORT_TARGET_CLASSIFICATION_20260902.csv'
P1=pd.read_csv(f6 if f6.exists() else ROAD/'p1_external_lopo_validation'/'P1_FIVE_COHORT_TARGET_CLASSIFICATION_20260902.csv')
V5=pd.read_csv(ROAD/'p3_third_complete_panel_search_20260903'/'GSE314416'/'transport_v5_quality_gate'/'GSE314416_V5_CROSS_SOURCE_CLASSIFICATION.csv')
P2=pd.read_csv(ROAD/'p2_technical_composition_20260902'/'GSE282881_ORIGINAL11_TECHNICAL_REPEATABILITY_AUDIT.csv')
P2=P2[['study_target','available','median_cv_across_27_clusters','p90_cv_across_27_clusters','median_range_ratio']].rename(columns={'study_target':'anchor_target','available':'technical_available'})
def canon(x): return ''.join(ch for ch in str(x).upper() if ch.isalnum())
alias={'CD279':'PD1','CD274B7H1PDL1':'PDL1','CD274':'PDL1','PD1':'PD1','PDL1':'PDL1','CD103':'CD103','CD8A':'CD8A','CD4':'CD4','HLADR':'HLADR','TIGIT':'TIGIT'}
def key(x): return alias.get(canon(x),canon(x))
P1['key']=P1.target.map(key); V5['key']=V5.target.map(key); P2['key']=P2.anchor_target.map(key)
rows=[]
for _,m in MAP.iterrows():
    t=m.gse334503_name; k=key(t); a=P1[P1.key.eq(k)]; b=V5[V5.key.eq(k)]; q=P2[P2.key.eq(k)]
    rows.append({'target':t,'rna_gene':m.rna_gene,'mapping_class':m.mapping_class,'analysis_role':'primary_direct' if bool(m.primary_direct_mapping_set) else 'sensitivity_non_direct','anchor_cohorts':int(a.cohorts_evaluated.iloc[0]) if len(a) else 0,'anchor_complete_cohorts':int(a.complete_panel_cohorts.iloc[0]) if len(a) else 0,'anchor_class':a.evidence_class.iloc[0] if len(a) else 'not_evaluated','third_cohort_class':b.cross_source_classification.iloc[0] if len(b) else 'not_evaluated','third_median_spearman':float(b.median_spearman.iloc[0]) if len(b) else None,'third_median_r2':float(b.median_r2.iloc[0]) if len(b) else None,'technical_repeatability_available':bool(q.technical_available.iloc[0]) if len(q) else False,'technical_median_cv':float(q.median_cv_across_27_clusters.iloc[0]) if len(q) and pd.notna(q.median_cv_across_27_clusters.iloc[0]) else None,'claim_boundary':'multi-cohort/context estimate; not universal; clinical evidence exploratory'})
df=pd.DataFrame(rows); df.to_csv(OUT/'INTEGRATED_EVIDENCE_MATRIX_V5_LATEST.csv',index=False)
s={'status':'LATEST_EVIDENCE_INTEGRATION_COMPLETE','targets':len(df),'source_versions':{'anchor':'P1 six-cohort final classification or five-cohort fallback','expanded_third':'GSE314416 v5 quality gate','technical':'P2 final technical repeatability audit','excluded':['GSE314416 baseline empty ADT','expanded v1/v2 invalid shuffle evidence','old integrated matrices']},'third_cohort_note':'GSE314416 is healthy PBMC follow-up, 17 independent participants; not a third cancer cohort','counts':df.third_cohort_class.value_counts().to_dict()}
(OUT/'INTEGRATED_EVIDENCE_MATRIX_V5_SUMMARY.json').write_text(json.dumps(s,indent=2,ensure_ascii=False),encoding='utf8'); print(json.dumps(s,indent=2,ensure_ascii=False))
