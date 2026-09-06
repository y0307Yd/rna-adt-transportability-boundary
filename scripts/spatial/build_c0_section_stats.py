import pandas as pd

src = r"C:\Users\Y\Documents\Codex\2026-08-13\jia\work\ai_multimodal_cholesterol_study\results\ai_extension\specificity_reliability_lopo\cd3e_conditional_section_metrics.csv"
d = pd.read_csv(src)

# rename to clear, manuscript-facing column names
out = d.rename(columns={
    "gsm": "section_id",
    "patient": "patient",
    "timepoint": "timepoint",
    "target": "target",
    "target_class": "target_class",
    "n_spots": "n_spots",
    "partial_spearman": "partial_spearman",
    "block_partial_rho": "block_partial_rho",
    "block_permutation_p": "block_permutation_p",
    "n_blocks": "n_blocks",
})[[
    "section_id","patient","timepoint","target","target_class","n_spots",
    "partial_spearman","block_partial_rho","block_permutation_p","n_blocks"
]]

# exclude technical replicate sections (per methods: GSM8789211/GSM8789214)
drops = ["GSM8789211","GSM8789214"]
out = out[~out["section_id"].isin(drops)].reset_index(drop=True)

out.to_csv(r"D:\rna-adt-transportability-boundary\supplementary\C0_spatial_patient_section_statistics.csv", index=False)

print("rows:", len(out))
print("sections:", out["section_id"].nunique())
print("patients:", sorted(out["patient"].unique()))
print("targets:", sorted(out["target"].unique()))
print()
print(out.groupby(["patient","timepoint"])["section_id"].nunique())
print()
print("n_spots range:", out["n_spots"].min(), "-", out["n_spots"].max())
