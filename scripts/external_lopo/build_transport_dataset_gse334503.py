from pathlib import Path
import gzip, json
import numpy as np, pandas as pd
from scipy.io import mmread
from scipy import sparse

p=Path(__file__).parent; out=p/'transport_dataset'; out.mkdir(exist_ok=True)
targets=['CD8A','CD4','CD14','CD19','CD56','HLA-DR','PD-1','PD-L1','CTLA4','TIGIT','CD103']
aliases={'CD8A':'CD8','CD4':'CD4','CD14':'CD14','CD19':'CD19','CD56':'CD56','HLA-DR':'HLA-DR','PD-1':'CD279','PD-L1':'CD274','CTLA4':'CD152','TIGIT':'TIGIT','CD103':'CD103'}
blocks_r=[]; blocks_a=[]; metas=[]; genes=None
for b in range(1,6):
    d=pd.read_csv(p/f'GSE334503_Batch{b}_hto_demultiplex.csv'); keep=np.flatnonzero(d.demux_status.eq('singlet').to_numpy())
    gf=next(p.glob(f'*_GEX_Batch{b}_features.tsv.gz')); raw=[x.strip().split('\t') for x in gzip.open(gf,'rt')]; gnames=[x[1] if len(x)>1 else x[0] for x in raw]
    if genes is None: genes=gnames
    G=mmread(next(p.glob(f'*_GEX_Batch{b}_matrix.mtx.gz'))).tocsr()[:,keep].T.tocsr(); blocks_r.append(G)
    af=[x.strip() for x in gzip.open(next(p.glob(f'*_ADT_Batch{b}_features.tsv.gz')),'rt')]
    A=mmread(next(p.glob(f'*_ADT_Batch{b}_matrix.mtx.gz'))).tocsr()
    idx=[next(i for i,n in enumerate(af) if aliases[t].lower() in n.lower()) for t in targets]
    Y=np.asarray(A[idx,:][:,keep].todense(),dtype=np.float32).T
    Y=Y-Y.mean(axis=1,keepdims=True)
    blocks_a.append(Y); metas.append(d.iloc[keep][['batch','barcode','patient_id','visit']])
R=sparse.vstack(blocks_r,format='csr'); Y=np.vstack(blocks_a); M=pd.concat(metas,ignore_index=True); sparse.save_npz(out/'rna_raw_cells_by_genes.npz',R); np.save(out/'adt_clr_cells_by_targets.npy',Y); np.save(out/'adt_raw_cells_by_targets.npy',Y); pd.DataFrame({'gene':genes}).to_csv(out/'rna_genes.csv',index=False); pd.DataFrame({'target':targets}).to_csv(out/'adt_targets.csv',index=False); M['donor_id']=M.patient_id.astype(str); M.to_csv(out/'cell_metadata.csv',index=False); (out/'build_manifest.json').write_text(json.dumps({'cells':len(M),'genes':len(genes),'patients':int(M.patient_id.nunique()),'targets':targets},indent=2),encoding='utf-8'); print(R.shape,Y.shape,M.patient_id.nunique())
