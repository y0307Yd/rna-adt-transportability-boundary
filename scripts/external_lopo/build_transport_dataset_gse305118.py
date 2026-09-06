from pathlib import Path
import gzip, json, re
import numpy as np, pandas as pd
from scipy.io import mmread
from scipy import sparse

ROOT=Path(__file__).parent.parent/'raw_extracted'; OUT=Path(__file__).parent/'transport_dataset'; OUT.mkdir(parents=True,exist_ok=True)
targets=['CD8A','CD4','CD14','CD19','CD56','HLA-DR','PD-1','PD-L1','CTLA4','TIGIT','CD103']
aliases={'CD8A':['CD8A','CD8'],'CD4':['CD4'],'CD14':['CD14'],'CD19':['CD19'],'CD56':['CD56','NCAM1'],'HLA-DR':['HLA-DR','HLA-DRA'],'PD-1':['CD279','PD-1'],'PD-L1':['CD274','PD-L1'],'CTLA4':['CD152','CTLA4'],'TIGIT':['TIGIT'],'CD103':['CD103','ITGAE']}
blocks_r=[]; blocks_a=[]; metas=[]; genes=None
for bf in sorted(ROOT.glob('*_barcodes.tsv.gz')):
    stem=bf.name[:-len('_barcodes.tsv.gz')]; ff=ROOT/(stem+'_features.tsv.gz'); mf=ROOT/(stem+'_matrix.mtx.gz')
    if not ff.exists() or not mf.exists(): continue
    rows=[x.rstrip('\n').split('\t') for x in gzip.open(ff,'rt')]
    gnames=[r[1] if len(r)>1 else r[0] for r in rows]; types=[r[2] if len(r)>2 else '' for r in rows]
    mat=mmread(mf).tocsr(); barcodes=[x.rstrip('\n') for x in gzip.open(bf,'rt')]
    if len(barcodes)!=mat.shape[1]: raise ValueError(stem+' barcode mismatch')
    gi=[i for i,t in enumerate(types) if t.lower() in ('gene expression','gene expression')]
    ai=[i for i,t in enumerate(types) if 'antibody' in t.lower() or 'adt' in t.lower()]
    if not gi: gi=list(range(mat.shape[0]))
    if not ai: ai=[i for i,n in enumerate(gnames) if any(a.lower()==n.lower() for aa in aliases.values() for a in aa)]
    if genes is None: genes=[gnames[i] for i in gi]
    if [gnames[i] for i in gi]!=genes: raise ValueError('GEX genes differ: '+stem)
    R=mat[gi,:].T.tocsr(); blocks_r.append(R)
    norm={re.sub(r'^(hu|human)[._-]*','',gnames[i],flags=re.I).lower():i for i in ai}; idx=[]
    for t in targets:
        hit=next((norm[a.lower()] for a in aliases[t] if a.lower() in norm),None)
        if hit is None: hit=next((i for k,i in norm.items() if any(k.startswith(a.lower()+'_') or k.startswith(a.lower()+'.') for a in aliases[t])),None)
        if hit is None and t=='HLA-DR': hit=next((i for k,i in norm.items() if k.replace('.','-') in ('hla-dr','hla.dra')),None)
        if hit is None: raise ValueError('missing '+t+' in '+stem)
        idx.append(hit)
    Y=np.asarray(mat[idx,:].todense(),dtype=np.float32).T; Y=Y-Y.mean(axis=1,keepdims=True); blocks_a.append(Y)
    m=re.match(r'GSM\d+_(Patient[^_]+)_([^_]+)',stem); patient=m.group(1) if m else 'unknown'; tissue=m.group(2) if m else 'unknown'; gsm=stem.split('_')[0]
    metas.append(pd.DataFrame({'barcode':barcodes,'patient_id':patient,'tissue':tissue,'gsm':gsm,'donor_id':patient,'library':stem}))
R=sparse.vstack(blocks_r,format='csr'); Y=np.vstack(blocks_a); M=pd.concat(metas,ignore_index=True)
sparse.save_npz(OUT/'rna_raw_cells_by_genes.npz',R); np.save(OUT/'adt_clr_cells_by_targets.npy',Y); np.save(OUT/'adt_raw_cells_by_targets.npy',Y)
pd.DataFrame({'gene':genes}).to_csv(OUT/'rna_genes.csv',index=False); pd.DataFrame({'target':targets}).to_csv(OUT/'adt_targets.csv',index=False); M.to_csv(OUT/'cell_metadata.csv',index=False)
(OUT/'build_manifest.json').write_text(json.dumps({'cells':int(len(M)),'genes':len(genes),'patients':int(M.patient_id.nunique()),'patient_tissue_units':int(M.library.nunique()),'targets':targets},indent=2),encoding='utf-8')
print(json.dumps({'shape':R.shape,'adt_shape':Y.shape,'patients':M.patient_id.nunique(),'units':M.library.nunique()},indent=2))
