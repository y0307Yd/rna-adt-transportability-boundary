# -*- coding: utf-8 -*-
"""重新生成 Figure 4 六队列版 600dpi TIFF + PDF（替换旧两队列 publication 版本）。"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGDIR = r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\figures"
BLUE = "#2166AC"; RED = "#B2182B"; PURPLE = "#7A3E9D"; GREY = "#888888"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.edgecolor": "#444444", "axes.linewidth": 0.6,
    "axes.titlesize": 11, "axes.titleweight": "bold",
})

cohorts = [
    ("GSE78220",   "Melanoma\n(anti-PD-1)",       27,  0.483, 0.552),
    ("GSE91061",   "Melanoma\n(anti-PD-1)",       49,  0.631, 0.109),
    ("GSE93157",   "Multi-cancer\n(anti-PD-1)",   65,  0.592, 0.122),
    ("IMvigor210", "Bladder\n(anti-PD-L1)",      298,  0.542, 0.147),
    ("Gide 2019",  "Melanoma\n(anti-PD-1\u00b1CTLA-4)", 73, 0.771, 0.0002),
    ("Braun 2020", "Kidney RCC\n(anti-PD-1)",    172,  0.516, 0.383),
]
labels = [c[0] for c in cohorts]
nvals  = [c[2] for c in cohorts]
aurocs = [c[3] for c in cohorts]
pvals  = [c[4] for c in cohorts]
x = np.arange(len(cohorts))
colors = [BLUE if c[0] != "Gide 2019" else RED for c in cohorts]

fig, axes = plt.subplots(2, 1, figsize=(8.2, 6.4),
                         gridspec_kw={"height_ratios": [1.35, 1.0]})

ax = axes[0]
ax.axhline(0.5, ls="--", color="grey", lw=1.0, zorder=0)
ax.bar(x, aurocs, 0.62, color=colors, edgecolor="#303030", lw=0.6, zorder=2)
for i, (au, nv, pv) in enumerate(zip(aurocs, nvals, pvals)):
    ax.text(i, au + 0.02, f"{au:.3f}", ha="center", va="bottom", fontsize=8, fontweight="bold")
    ax.text(i, -0.13, f"n={nv}\nP={pv:.3f}", ha="center", va="top", fontsize=6.8, color="#333333")
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylim(-0.22, 0.95)
ax.set_ylabel("AUROC", fontsize=9)
ax.set_title("Frozen four-gene score in six independent ICI cohorts (n=684 evaluable)", fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
ax.annotate("* infiltration proxy,\nnot a biomarker", xy=(4, 0.771), xytext=(4.35, 0.90),
            fontsize=7.5, color=RED, fontweight="bold",
            arrowprops=dict(arrowstyle="->", color=RED, lw=0.8))

ax = axes[1]
ax.axhline(0.05, ls="--", color=GREY, lw=1.0, zorder=0)
ax.text(len(cohorts)-0.3, 0.05, " P=0.05", va="bottom", ha="right", fontsize=7, color=GREY)
for i, pv in enumerate(pvals):
    c = RED if pv < 0.05 else PURPLE
    ax.scatter(i, pv, s=42, color=c, zorder=3, edgecolor="#303030", lw=0.5)
    ax.text(i, pv + 0.03, f"{pv:.3f}", ha="center", va="bottom", fontsize=7.5)
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylim(0, 0.72)
ax.set_ylabel("Label-permutation P", fontsize=9)
ax.set_title("Significance (all non-significant; Gide is an infiltration-proxy artefact)", fontsize=9)
ax.spines[["top", "right"]].set_visible(False)

fig.tight_layout()

# 600dpi TIFF + PDF（六队列版，覆盖旧 publication 版）
tiff_path = os.path.join(FIGDIR, "Figure4_publication.tiff")
pdf_path  = os.path.join(FIGDIR, "Figure4_publication.pdf")
png_path  = os.path.join(FIGDIR, "Figure4_publication.png")
fig.savefig(tiff_path, dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
fig.savefig(pdf_path,  bbox_inches="tight")
fig.savefig(png_path,  dpi=600, bbox_inches="tight")
plt.close(fig)
print("saved Figure4_publication.tiff/pdf/png @600dpi (six-cohort)")

# 同步 Figure4_source_data.csv（六队列完整统计口径）
rows = [["cohort","cancer","n_evaluable","auroc","perm_p"]]
for c in cohorts:
    cancer_label = c[1].replace("\n"," ").replace("\u00b1","+/-")
    rows.append([c[0], cancer_label, c[2], round(c[3],4), c[4]])
with open(os.path.join(FIGDIR, "Figure4_source_data.csv"), "w", newline="", encoding="utf-8") as f:
    csv.writer(f).writerows(rows)
print("saved Figure4_source_data.csv (six-cohort)")
