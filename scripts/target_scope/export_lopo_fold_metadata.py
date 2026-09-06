import json
from pathlib import Path
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge

root=Path(r'C:/Users/Y/Documents/Codex/2026-08-13/jia/work/ai_multimodal_cholesterol_study')
raw=root/'results/ai_extension/failure_sources/failure_source_matrix.csv'
out=Path(r'D:/ai_multimodal_cholesterol_study_outputs/revision-roadmap/lopo_lodo_audit/fold_level'); out.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(raw).groupby('protein',as_index=False).mean(numeric_only=True)
feats=['target_cognate_rna_adt_rho','target_within_cohort_loo_rho']; ycol='ai_cross_median_donor_rho'; rows=[]
for i,r in d.iterrows():
    tr=d.index!=i; scaler=StandardScaler().fit(d.loc[tr,feats]); model=Ridge(alpha=1.0).fit(scaler.transform(d.loc[tr,feats]),d.loc[tr,ycol]); pred=model.predict(scaler.transform(d.loc[[i],feats]))[0]
    rows.append({'fold_id':int(i),'held_out_target':r.protein,'training_targets':';'.join(d.loc[tr,'protein'].tolist()),'selected_features':';'.join(feats),'scaler_mean':json.dumps(scaler.mean_.tolist()),'scaler_scale':json.dumps(scaler.scale_.tolist()),'ridge_alpha':1.0,'ridge_solver':'auto','observed_value':float(r[ycol]),'predicted_value':float(pred)})
pd.DataFrame(rows).to_csv(out/'lopo_fold_level_metadata.csv',index=False)
(out/'manifest.json').write_text(json.dumps({'source_matrix':str(raw),'n_folds':len(rows),'features':feats,'outcome':ycol,'ridge_alpha':1.0,'scaler':'StandardScaler fit within training fold'},indent=2),encoding='utf-8')
print('folds',len(rows))
