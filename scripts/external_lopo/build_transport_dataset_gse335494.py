from pathlib import Path
import gzip,json,re
import numpy as np,pandas as pd
from scipy import sparse
from scipy.io import mmread
p=Path(__file__).parent.parent/'raw_extracted'; out=Path(__file__).parent/'transport_dataset'; out.mkdir(parents=True,exist_ok=True)
targets=['CD8A','CD4','CD14','CD19','CD56','HLA-DR','PD-1','PD-L1','TIGIT','CD103']
aliases={'CD8A':['CD8'],'CD4':['CD4'],'CD14':['CD14'],'CD19':['CD19'],'CD56':['CD56'],'HLA-DR':['HLA-DR'],'PD-1':['CD279'],'PD-L1':['CD274'],'TIGIT':['TIGIT'],'CD103':['CD103']}
def lines(f): return [x.rstrip().split('\t') for x in gzip.open(f,'rt') if x.strip()]
files=sorted(p.glob('*_features.tsv.gz')); gene_lists=[]
for f in files:
 r=lines(f); gene_lists.append([x[1] for x in r if len(x)>2 and x[2]=='Gene Expression'])
common=set(gene_lists[0])
for g in gene_lists[1:]: common &= set(g)
genes=[g for g in gene_lists[0] if g in common]; blocks=[]; ys=[]; metas=[]; exclusions=[]
for f in files:
 stem=f.name[:-len('_features.tsv.gz')]; r=lines(f); names=[x[1] for x in r]; types=[x[2] if len(x)>2 else '' for x in r]; lookup={n:i for i,n in enumerate(names)}; lookup.update({x[0]:i for i,x in enumerate(r)}); gi=[lookup[g] for g in genes]
 mat=mmread(p/(stem+'_matrix.mtx.gz')).tocsr(); bc=[x[0] for x in lines(p/(stem+'_barcodes.tsv.gz'))]
 if mat.shape[1]!=len(bc): raise ValueError(stem+' barcode mismatch')
 idx=[]; missing=[]
 for t in targets:
  hit=next((lookup[a] for a in aliases[t] if a in lookup and types[lookup[a]]=='Antibody Capture'),None)
  if hit is None: missing.append(t)
  else: idx.append(hit)
 if missing:
  exclusions.append({'sample':stem,'reason':'missing_locked_ADT','missing_targets':';'.join(missing),'cells':len(bc)})
  continue
 blocks.append(mat[gi,:].T.tocsr())
 A=np.asarray(mat[idx,:].todense(),dtype=np.float32).T; A=A-A.mean(axis=1,keepdims=True); ys.append(A)
 m=re.search(r'_(S\d+)-([A-Z]+)$',stem); donor=m.group(1) if m else stem; time=m.group(2) if m else 'unknown'; gsm=stem.split('_')[0]
 metas.append(pd.DataFrame({'barcode':[stem+':'+x for x in bc],'raw_barcode':bc,'donor_id':donor,'timepoint':time,'sample':stem,'gsm':gsm}))
R=sparse.vstack(blocks,format='csr'); Y=np.vstack(ys); M=pd.concat(metas,ignore_index=True)
sparse.save_npz(out/'rna_raw_cells_by_genes.npz',R); np.save(out/'adt_clr_cells_by_targets.npy',Y); np.save(out/'adt_raw_cells_by_targets.npy',Y); pd.DataFrame({'gene':genes}).to_csv(out/'rna_genes.csv',index=False); pd.DataFrame({'target':targets}).to_csv(out/'adt_targets.csv',index=False); M.to_csv(out/'cell_metadata.csv',index=False)
pd.DataFrame(exclusions).to_csv(out/'excluded_samples.csv',index=False)
manifest={'cells':len(M),'genes_common_all_samples':len(genes),'independent_patients':int(M.donor_id.nunique()),'sample_timepoints':int(M['sample'].nunique()),'excluded_samples':len(exclusions),'targets_evaluated':targets,'structurally_missing':['CTLA4'],'source_archive_sha256':'40B655EEAB5C939ED0762074C64F47B1034EBDB72D3D941BCDBAFDC472869417'}
(out/'build_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8'); print(json.dumps(manifest,indent=2))
