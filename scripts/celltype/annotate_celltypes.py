"""Cell type annotation for GSE334503 / GSE335494 (43-target transport cohorts).

Marker-score based lineage annotation. No external reference needed.
Outputs: cell_type column appended to a copy of cell_metadata for each cohort.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse

BASE = Path(r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap")
OUT = BASE / "celltype_annotation_20260905"
OUT.mkdir(parents=True, exist_ok=True)

# lineage -> marker genes (RNA symbols, matched case-insensitively)
LINEAGES = {
    "CD4_T":  ["CD3D", "CD3E", "CD3G", "CD4"],
    "CD8_T":  ["CD3D", "CD3E", "CD8A", "CD8B"],
    "NK":     ["NKG7", "GNLY", "KLRD1", "KLRB1", "NCAM1"],
    "B":      ["MS4A1", "CD19", "CD79A", "CD79B", "IGHM", "MZB1", "TCL1A"],
    "Mono":   ["CD14", "LYZ", "S100A8", "S100A9", "ITGAM"],
    "cDC":    ["CLEC9A", "CD1C", "FCER1A", "IRF8"],
    "pDC":    ["CLEC4C", "LILRA4", "TCF4"],
    "Platelet": ["PPBP"],
}


def annotate(cohort):
    td = BASE / f"p1_external_lopo_validation/{cohort}/v1_reconstruction/transport_dataset"
    genes = pd.read_csv(td / "rna_genes.csv").gene.astype(str).tolist()
    gmap = {g.upper(): i for i, g in enumerate(genes)}

    # resolve marker indices
    marker_to_idx = {}
    all_markers = []
    for lin, ms in LINEAGES.items():
        for m in ms:
            if m.upper() in gmap and m.upper() not in marker_to_idx:
                marker_to_idx[m.upper()] = gmap[m.upper()]
                all_markers.append(m)

    print(f"[{cohort}] {len(all_markers)} unique markers resolved")

    # load RNA marker columns (sparse -> dense, only marker columns)
    R = sparse.load_npz(td / "rna_raw_cells_by_genes.npz")
    idx = [marker_to_idx[m.upper()] for m in all_markers]
    X = R[:, idx].toarray().astype(np.float32)
    del R
    n_cells = X.shape[0]

    # log1p + per-gene z-score (across cells)
    Xl = np.log1p(X)
    mu = Xl.mean(axis=0)
    sd = Xl.std(axis=0)
    sd[sd == 0] = 1.0
    Z = (Xl - mu) / sd
    del X, Xl

    # module score per lineage = mean z-score of its markers
    lin_names = list(LINEAGES.keys())
    scores = np.zeros((n_cells, len(lin_names)), dtype=np.float32)
    for i, lin in enumerate(lin_names):
        cols = [marker_to_idx[m.upper()] for m in LINEAGES[lin]]
        # position of these markers in the resolved marker order
        pos = [all_markers.index(m) for m in LINEAGES[lin]]
        scores[:, i] = Z[:, pos].mean(axis=1)

    # assignment: argmax with gap + minimum thresholds
    order = np.argsort(-scores, axis=1)
    top1 = scores[np.arange(n_cells), order[:, 0]]
    top2 = scores[np.arange(n_cells), order[:, 1]]
    gap = top1 - top2
    assign_idx = order[:, 0]

    MIN_SCORE = 0.5
    MIN_GAP = 0.25
    celltype = np.array(lin_names)[assign_idx].astype(object)
    ambiguous = (top1 < MIN_SCORE) | (gap < MIN_GAP)
    celltype[ambiguous] = "Ambiguous"

    meta = pd.read_csv(td / "cell_metadata.csv")
    meta["cell_type"] = celltype
    meta["top1_score"] = top1
    meta["top2_gap"] = gap

    out = OUT / f"{cohort}_cell_metadata_annotated.csv"
    meta.to_csv(out, index=False)
    print(f"[{cohort}] annotated {n_cells} cells -> {out.name}")
    print(meta.cell_type.value_counts().to_string())
    print()
    return meta


if __name__ == "__main__":
    for c in ["GSE334503", "GSE335494"]:
        annotate(c)
