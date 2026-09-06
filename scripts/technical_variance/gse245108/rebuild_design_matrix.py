import json
import pandas as pd

# rebuild sample-level design matrix WITH cell_type counts from tnc_meta_full.json
meta = json.load(open(r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\technical_biological_decomposition_data_v1_20260903\GSE245108_tnc_parsed\tnc_meta_full.json", encoding="utf-8"))
df = pd.DataFrame(meta)

# normalize concentration to the 0.25/0.5/1/2/4 scheme
def norm_conc(c):
    c = str(c).strip()
    if c in ("0.25","0.5","1","2","4"):
        return c
    return c

# group by donor/concentration/rep -> n cells + n cell types + top cell type
rows = []
for (donor, conc, rep), g in df.groupby(["donor","concentration","rep"]):
    n_cells = len(g)
    n_celltypes = g["cluster"].nunique()
    top_ct = g["cluster"].value_counts().idxmax()
    top_ct_n = g["cluster"].value_counts().max()
    rows.append({
        "donor": donor,
        "concentration": str(conc),
        "rep": rep,
        "n_cells": n_cells,
        "n_cell_types": n_celltypes,
        "top_cell_type": top_ct,
        "top_cell_type_n_cells": int(top_ct_n),
    })

D = pd.DataFrame(rows)
# sort donor by known order, concentration numeric, rep
order = ["WM22","WF32","BF32","BM32"]
D["_donor_ord"] = D["donor"].map({d:i for i,d in enumerate(order)})
D["_conc_ord"] = D["concentration"].astype(float)
D["_rep_ord"] = D["rep"].str.replace("rep","").astype(int)
D = D.sort_values(["_donor_ord","_conc_ord","_rep_ord"]).drop(columns=["_donor_ord","_conc_ord","_rep_ord"]).reset_index(drop=True)

D.to_csv(r"D:\rna-adt-transportability-boundary\supplementary\D_gse245108_sample_design_matrix.csv", index=False)
print(D.to_string(index=False))
print()
print("total cells:", D["n_cells"].sum())
print("n cell types total:", df["cluster"].nunique())
