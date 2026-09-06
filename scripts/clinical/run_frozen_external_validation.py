from pathlib import Path
import json, hashlib, numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, roc_auc_score

OUT=Path(__file__).parent; OUT.mkdir(exist_ok=True)
TRAIN=Path(r'C:\Users\Y\Documents\Codex\2026-08-13\jia\work\ai_multimodal_cholesterol_study\results\ai_extension\failure_sources\failure_source_matrix.csv')
ROOT=Path(r'D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\p1_external_lopo_validation')
FILES={
'GSE305118_AML':ROOT/'GSE305118/v1_transport_bc_to_gse305118/target_metrics.csv',
'GSE309625_CRLM':ROOT/'GSE309625/v1_transport_bc_to_gse309625/target_metrics.csv',
'GSE316782_NSCLC':ROOT/'GSE316782/v1_transport_bc_to_gse316782_tumour/target_metrics.csv',
'GSE317605_BILIARY':ROOT/'GSE317605/v1_transport_bc_to_gse317605/target_metrics.csv',
'GSE334503_EXTERNAL':ROOT/'GSE334503/v1_transport_bc_to_gse334503_memorysafe_v2/target_metrics.csv',
'GSE335494_MELANOMA':ROOT/'GSE335494/v1_transport_bc_to_gse335494/target_metrics.csv'}
features=['target_cognate_rna_adt_rho','target_within_cohort_loo_rho']; outcome='ai_cross_median_donor_rho'
raw=pd.read_csv(TRAIN); train=raw.groupby('protein',as_index=False).mean(numeric_only=True)
sc=StandardScaler().fit(train[features]); model=Ridge(alpha=1.0).fit(sc.transform(train[features]),train[outcome])
pred=pd.DataFrame({'target':train.protein,'frozen_prediction':model.predict(sc.transform(train[features]))})
def canon(x):
 s=''.join(c for c in str(x).upper() if c.isalnum()); return {'CD279':'PD1','CD274':'PDL1','HLADR':'HLADR','CD103INTEGRINE':'CD103'}.get(s,s)
pred['key']=pred.target.map(canon); rows=[]
for cohort,f in FILES.items():
 d=pd.read_csv(f); d=d[d['mode'].astype(str).str.lower().eq('observed')] if 'mode' in d else d
 d['key']=d.target.map(canon); m=d.merge(pred[['key','frozen_prediction']],on='key',how='inner')
 y=m.median_within_donor_spearman.astype(float); score=m.frozen_prediction.astype(float); binary=(y>0.05).astype(int)
 rho=float(spearmanr(score,y).statistic) if len(m)>=3 else np.nan
 rng=np.random.default_rng(20260903); null=np.array([spearmanr(score,rng.permutation(y)).statistic for _ in range(10000)]) if len(m)>=3 else np.array([])
 p=float((1+(np.abs(null)>=abs(rho)).sum())/(len(null)+1)) if len(null) else np.nan
 auc=float(roc_auc_score(binary,score)) if binary.nunique()==2 else np.nan
 rows.append({'cohort':cohort,'targets':len(m),'spearman':rho,'permutation_p':p,'mae':float(mean_absolute_error(y,score)),'auc_transport_flag':auc,'positive_targets':int(binary.sum())})
 m.assign(cohort=cohort,observed_external=y).to_csv(OUT/f'{cohort}_FROZEN_PREDICTIONS.csv',index=False)
res=pd.DataFrame(rows); res.to_csv(OUT/'FROZEN_LOPO_EXTERNAL_VALIDATION_BY_COHORT.csv',index=False)
manifest={'status':'FROZEN_EXTERNAL_OUTCOME_VALIDATION_COMPLETE','training_targets':len(train),'features':features,'alpha':1.0,'external_cohorts':list(FILES),'seed':20260903,'permutations':10000,'coefficient':model.coef_.tolist(),'intercept':float(model.intercept_),'boundary':'validates external outcome generalization for previously seen targets; does not validate prediction for unseen proteins','training_sha256':hashlib.sha256(TRAIN.read_bytes()).hexdigest()}
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8'); print(res.to_string(index=False)); print(json.dumps(manifest,indent=2))
