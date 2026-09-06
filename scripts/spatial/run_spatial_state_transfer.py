from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.spatial import cKDTree
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler


STATES = {"CD8A": "CD8A", "HLA-DR": "HLA-DRA", "PD-1": "PDCD1", "CD103": "ITGAE"}
MARKERS = {
    "tcell": ["CD3D", "CD3E", "CD8A", "CD8B", "TRBC1", "TRBC2", "CCL5"],
    "myeloid": ["LST1", "TYROBP", "FCER1G", "C1QA", "C1QB", "C1QC", "CD68"],
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT14"],
}


def decode(x): return [v.decode() if isinstance(v, bytes) else str(v) for v in x]
def open_text(path):
    with open(path, "rb") as h: magic = h.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, encoding="utf-8")


def parse_titles(path):
    out, accession = {}, None
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("^SAMPLE = "): accession = line.split("=", 1)[1].strip()
        elif accession and line.startswith("!Sample_title = "): out[accession] = line.split("=", 1)[1].strip()
    return out


def patient_from_title(title, accession):
    m = re.search(r"patient\s+([^,]+)", title, re.I)
    return m.group(1).strip() if m else accession


def positions(path, valid):
    out = {}
    with open_text(path) as h:
        for row in csv.reader(h):
            if row and row[0] in valid and len(row) >= 6:
                out[row[0]] = (float(row[4]), float(row[5]))
    return out


def score(matrix, lookup, genes):
    idx = [lookup[g] for g in genes if g in lookup]
    return np.asarray(matrix[idx].mean(axis=0)).ravel() if idx else np.zeros(matrix.shape[1])


def train_models(source_dir: Path, top_genes=1000):
    x = sparse.load_npz(source_dir / "rna_raw_cells_by_genes.npz").tocsr().astype(float)
    y = np.load(source_dir / "adt_clr_cells_by_targets.npy").astype(float)
    meta = pd.read_csv(source_dir / "cell_metadata.csv")
    genes = pd.read_csv(source_dir / "rna_genes.csv").gene.astype(str).to_numpy()
    targets = pd.read_csv(source_dir / "adt_targets.csv").target.astype(str).tolist()
    lib = np.asarray(x.sum(axis=1)).ravel(); x = x.multiply((10000.0 / lib)[:, None]).tocsr(); x.data = np.log1p(x.data)
    det = np.asarray(x.getnnz(axis=0)).ravel(); mean = np.asarray(x.mean(axis=0)).ravel(); var = np.maximum(0, np.asarray(x.power(2).mean(axis=0)).ravel() - mean**2)
    elig = np.flatnonzero(det >= max(20, int(x.shape[0] * 0.01))); ranked = elig[np.argsort(var[elig])[-top_genes:]]
    fixed = np.flatnonzero(np.isin(genes, ["ABCG1", "CD8A", "CD4", "HLA-DRA", "PDCD1", "ITGAE", "FDPS", "EPCAM", "ERBB2"]))
    selected = np.unique(np.concatenate([ranked, fixed]))
    scaler = StandardScaler(with_mean=False); xs = scaler.fit_transform(x[:, selected])
    models = {}
    for marker in STATES:
        j = targets.index(marker); ys = (y[:, j] - y[:, j].mean()) / y[:, j].std()
        models[marker] = Ridge(alpha=10.0, solver="lsqr").fit(xs, ys)
    return {
        "models": models,
        "scaler": scaler,
        "selected_genes": genes[selected],
        "source_detection": (det[selected] / x.shape[0]),
    }


def selected_spatial_matrix(norm, names, selected_genes):
    lookup = {}
    for i, gene in enumerate(names):
        lookup.setdefault(gene, []).append(i)
    columns = []
    present = []
    for gene in selected_genes:
        idx = lookup.get(gene, [])
        if not idx:
            columns.append(sparse.csr_matrix((norm.shape[1], 1)))
            present.append(False)
        else:
            columns.append(sparse.csr_matrix(norm[idx].sum(axis=0).T))
            present.append(True)
    return sparse.hstack(columns, format="csr"), np.asarray(present)


def standardized_mean_difference(values, high, low):
    if high.sum() < 5 or low.sum() < 5:
        return float("nan")
    pooled_var = (
        (high.sum() - 1) * values[high].var(ddof=1)
        + (low.sum() - 1) * values[low].var(ddof=1)
    ) / (high.sum() + low.sum() - 2)
    if not np.isfinite(pooled_var) or pooled_var <= 0:
        return float("nan")
    return float((values[high].mean() - values[low].mean()) / np.sqrt(pooled_var))


def analyze(h5_path, pos_path, cohort, title, bundle):
    with h5py.File(h5_path, "r") as h:
        m = h["matrix"]; names = decode(m["features"]["name"][:]); barcodes = decode(m["barcodes"][:])
        raw = sparse.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=tuple(m["shape"][:]))
    pos = positions(pos_path, set(barcodes)); keep = np.asarray([i for i,b in enumerate(barcodes) if b in pos])
    barcodes = [barcodes[i] for i in keep]; raw = raw[:, keep]; lib = np.asarray(raw.sum(axis=0)).ravel(); lib[lib == 0] = 1
    norm = raw.multiply(1e4 / lib).log1p().tocsc(); lookup = {n:i for i,n in enumerate(names)}
    tcell, myeloid = score(norm, lookup, MARKERS["tcell"]), score(norm, lookup, MARKERS["myeloid"])
    epithelial = score(norm, lookup, MARKERS["epithelial"])
    xy = np.asarray([pos[b] for b in barcodes])
    neigh = cKDTree(xy).query(xy, k=min(7, len(xy)))[1][:, 1:]
    epithelial_neigh = epithelial[neigh].mean(axis=1)
    myeloid_neigh = myeloid[neigh].mean(axis=1)
    tcell_anchor = tcell >= np.median(tcell)

    spatial, present = selected_spatial_matrix(norm, names, bundle["selected_genes"])
    transformed = bundle["scaler"].transform(spatial)
    rows = []
    for marker, model in bundle["models"].items():
        pred = model.predict(transformed)
        q1, q3 = np.quantile(pred[tcell_anchor], [0.25, 0.75])
        high = tcell_anchor & (pred >= q3)
        low = tcell_anchor & (pred <= q1)
        accession = re.match(r"GSM\d+", h5_path.name).group(0)
        rows.append({
            "cohort": cohort,
            "accession": accession,
            "patient": patient_from_title(title, accession),
            "state": marker,
            "n_spots": len(barcodes),
            "n_tcell_anchor_spots": int(tcell_anchor.sum()),
            "selected_gene_coverage": float(present.mean()),
            "predicted_state_epithelial_neighborhood_smd": standardized_mean_difference(epithelial_neigh, high, low),
            "predicted_state_myeloid_neighborhood_smd": standardized_mean_difference(myeloid_neigh, high, low),
            "predicted_state_mean_high": float(pred[high].mean()),
            "predicted_state_mean_low": float(pred[low].mean()),
        })
    return rows


def main():
    p=argparse.ArgumentParser(); p.add_argument("source_dir",type=Path); p.add_argument("spatial_root",type=Path); p.add_argument("output_dir",type=Path); p.add_argument("--top-genes",type=int,default=1000); args=p.parse_args(); args.output_dir.mkdir(parents=True,exist_ok=True)
    bundle=train_models(args.source_dir,args.top_genes); jobs=[("GSE210616",args.spatial_root/"GSE210616_family.soft.txt",args.spatial_root/"gse210616_h5",args.spatial_root/"gse210616_spatial"),("GSE213688",args.spatial_root/"GSE213688_family.soft.txt",args.spatial_root/"gse213688_h5",args.spatial_root/"gse213688_spatial")]; rows=[]
    for cohort,soft,h5dir,posdir in jobs:
        titles=parse_titles(soft)
        for h5 in sorted(h5dir.glob("*.h5")):
            acc=re.match(r"GSM\d+",h5.name).group(0); stem=h5.name.replace("_filtered_feature_bc_matrix.h5",""); pos=posdir/f"{stem}_tissue_positions_list.csv.gz"
            if not pos.exists():
                raise FileNotFoundError(pos)
            rows.extend(analyze(h5,pos,cohort,titles.get(acc,""),bundle))
    df=pd.DataFrame(rows)
    df.to_csv(args.output_dir/"spatial_state_transfer_section.csv",index=False)
    patient=df.groupby(["cohort","patient","state"],as_index=False).agg(
        primary_smd=("predicted_state_epithelial_neighborhood_smd","mean"),
        myeloid_smd=("predicted_state_myeloid_neighborhood_smd","mean"),
        n_sections=("accession","nunique"),
        selected_gene_coverage=("selected_gene_coverage","median"),
    )
    patient.to_csv(args.output_dir/"spatial_state_transfer_patient.csv",index=False)
    cohort_rows=[]
    for (cohort,state), group in patient.groupby(["cohort","state"]):
        values=group.primary_smd.dropna().to_numpy()
        p_value=float(stats.wilcoxon(values).pvalue) if len(values) and np.any(values != 0) else float("nan")
        cohort_rows.append({"cohort":cohort,"state":state,"n_patients":len(values),"median_primary_smd":float(np.median(values)) if len(values) else float("nan"),"positive_fraction":float(np.mean(values>0)) if len(values) else float("nan"),"wilcoxon_p":p_value})
    cohort_summary=pd.DataFrame(cohort_rows)
    cohort_summary.to_csv(args.output_dir/"spatial_state_transfer_cohort.csv",index=False)
    patient_keys=patient[["cohort","patient"]].drop_duplicates()
    summary={"source":"GSE199219 T-cell CITE-seq","states":list(STATES),"spatial_cohorts":["GSE210616","GSE213688"],"sections":int(df[["cohort","accession"]].drop_duplicates().shape[0]),"patients":int(len(patient_keys)),"primary_endpoint":"epithelial-neighbourhood SMD for high versus low predicted state within T-cell-rich anchor spots","interpretation":"exploratory patient-level spatial transfer; predicted protein states are not measured spatial protein"}; (args.output_dir/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8"); print(json.dumps(summary,indent=2))


if __name__ == "__main__": main()
