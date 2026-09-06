from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

src = Path(r'D:/ai_multimodal_cholesterol_study_outputs/v5_external_transport/validated/TRANSPORT_SHUFFLE_SUPERIORITY_MATRIX_20260901.csv')
out = Path(r'D:/ai_multimodal_cholesterol_study_outputs/revision-roadmap/figures')
out.mkdir(parents=True, exist_ok=True)
d = pd.read_csv(src)
d['analysis'] = d['analysis'].replace({'BC_to_GSE316782_NSCLC_tumour':'BC→NSCLC','BC_to_GSE316782_NSCLC_tumour_v2':'BC→NSCLC (full)','GSE316782_NSCLC_tumour_to_BC':'NSCLC→BC','BC_to_GSE317605_biliary_partial9':'BC→BTC (partial)','GSE317605_biliary_partial9_to_BC':'BTC→BC (partial)','bc_to_crlm_tumour_only_sensitivity':'BC→CRLM (sens.)','bc_to_dlbcl_partial5':'BC→DLBCL (partial)','bc_to_hcc':'BC→HCC'})
tab = d.pivot_table(index='target', columns='analysis', values='delta_observed_minus_shuffled')
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.linewidth':0.6})
fig, ax = plt.subplots(figsize=(12, 5))
im=ax.imshow(tab.to_numpy(),aspect='auto',cmap='RdBu_r',vmin=-0.2,vmax=0.4)
ax.set_xticks(range(len(tab.columns))); ax.set_xticklabels(tab.columns,rotation=45,ha='right',fontsize=7)
ax.set_yticks(range(len(tab.index))); ax.set_yticklabels(tab.index)
for i in range(tab.shape[0]):
 for j in range(tab.shape[1]):
  v=tab.iloc[i,j]
  if pd.notna(v): ax.text(j,i,f'{v:.2f}',ha='center',va='center',fontsize=6)
fig.colorbar(im,ax=ax,label='Observed - shuffled Spearman',fraction=.025,pad=.02)
ax.set_title('External RNA–ADT transportability')
ax.set_xlabel('Transport direction; BTC/DLBCL panels are partial')
ax.set_ylabel('Target')
fig.tight_layout()
fig.savefig(out / 'Figure1_external_transport_superiority.png', dpi=300)
fig.savefig(out / 'Figure1_external_transport_superiority.pdf')
tab.to_csv(out / 'Figure1_source_data.csv')
