# -*- coding: utf-8 -*-
"""
重绘 publication figures（2026-09-04）—— 用 matplotlib 替代 R（Rscript NOT FOUND）。
覆盖：
  Figure 4 : 六队列 ICI score 外部评估（升级自旧两队列版）
  Figure 6 : 研究流程图（补技术噪声/可识别性/六队列证伪新节点）
  Figure 7 : 技术噪声分层方差分解（四数据集横向 + donor 双峰 + rep vs conc 散点）
  Figure 8 : 可识别性模拟（donor bias x design + 功效曲线 + CLR/CPM 归一化偏差）
  Figure 9 : PCA 特征值谱（scree + 累计方差 + Neff）
输出 PNG(300dpi) + SVG 至 figures 目录。
"""
import csv, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

FIGDIR = r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\figures"
os.makedirs(FIGDIR, exist_ok=True)

# 统一风格
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.edgecolor": "#444444",
    "axes.linewidth": 0.6,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
})

BLUE = "#2166AC"; RED = "#B2182B"; PURPLE = "#7A3E9D"
GREEN = "#2A9D8F"; ORANGE = "#E07A5F"; GREY = "#888888"

def save(fig, name, w, h):
    for ext, dpi in [("png", 300), ("svg", None)]:
        p = os.path.join(FIGDIR, name + "." + ext)
        if ext == "svg":
            fig.savefig(p, format="svg", bbox_inches="tight")
        else:
            fig.savefig(p, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print("saved", name)

# ============================================================
# Figure 4 : six-cohort ICI score external evaluation
# ============================================================
def fig4():
    cohorts = [
        ("GSE78220",   "Melanoma\n(anti-PD-1)",       27,  0.483, 0.552),
        ("GSE91061",   "Melanoma\n(anti-PD-1)",       49,  0.631, 0.109),
        ("GSE93157",   "Multi-cancer\n(anti-PD-1)",   65,  0.592, 0.122),
        ("IMvigor210", "Bladder\n(anti-PD-L1)",      298,  0.542, 0.147),
        ("Gide 2019",  "Melanoma\n(anti-PD-1\u00b1CTLA-4)", 73, 0.771, 0.0002),
        ("Braun 2020", "Kidney RCC\n(anti-PD-1)",    172,  0.516, 0.383),
    ]
    labels   = [c[0] for c in cohorts]
    nvals    = [c[2] for c in cohorts]
    aurocs   = [c[3] for c in cohorts]
    pvals    = [c[4] for c in cohorts]
    x = np.arange(len(cohorts))
    colors = [BLUE if c[0] != "Gide 2019" else RED for c in cohorts]

    fig, axes = plt.subplots(2, 1, figsize=(8.2, 6.4),
                             gridspec_kw={"height_ratios": [1.35, 1.0]})

    # Panel A: AUROC
    ax = axes[0]
    ax.axhline(0.5, ls="--", color="grey", lw=1.0, zorder=0)
    bars = ax.bar(x, aurocs, 0.62, color=colors, edgecolor="#303030", lw=0.6, zorder=2)
    for i, (au, nv, pv) in enumerate(zip(aurocs, nvals, pvals)):
        ax.text(i, au + 0.02, f"{au:.3f}", ha="center", va="bottom", fontsize=8, fontweight="bold")
        ax.text(i, -0.13, f"n={nv}\nP={pv:.3f}", ha="center", va="top", fontsize=6.8, color="#333333")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylim(-0.22, 0.95)
    ax.set_ylabel("AUROC", fontsize=9)
    ax.set_title("Frozen four-gene score in six independent ICI cohorts (n=684 evaluable)", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    # Gide annotation
    ax.annotate("* infiltration proxy,\nnot a biomarker", xy=(4, 0.771), xytext=(4.35, 0.90),
                fontsize=7.5, color=RED, fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=RED, lw=0.8))

    # Panel B: permutation P
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
    save(fig, "Figure4_ici_score_external_evaluation_v2", 8.2, 6.4)

    # source data csv
    rows = [["cohort","cancer","n_evaluable","auroc","perm_p"]]
    for c in cohorts:
        cancer_label = c[1].replace("\n"," ").replace("\u00b1","+/-")
        rows.append([c[0], cancer_label, c[2], round(c[3],4), c[4]])
    with open(os.path.join(FIGDIR, "Figure4_source_data.csv"), "w", newline="", encoding="utf-8") as f:
        csv.writer(f).writerows(rows)

# ============================================================
# Figure 6 : research framework / study flow (new nodes)
# ============================================================
def fig6():
    fig, ax = plt.subplots(figsize=(11, 5.2))
    ax.set_xlim(0, 12); ax.set_ylim(0, 4.4); ax.axis("off")

    # Row 1: main pipeline
    main = [
        (1.25, "RNA input\npatient / cell",        "#D9EAF7"),
        (3.55, "RNA-ADT mapping\nregularized ridge","#DDEAD7"),
        (5.85, "Held-out transport\ndonor / cohort", "#FCE8C3"),
        (8.15, "Evidence integration\nshuffle + spatial","#E8DFF5"),
        (10.45,"Clinical extension\nfrozen score", "#F7D6D0"),
    ]
    w = 1.85
    for x, label, fc in main:
        ax.add_patch(FancyBboxPatch((x-w/2, 3.1), w, 0.95,
                     boxstyle="round,pad=0.02,rounding_size=0.12",
                     fc=fc, ec="#303030", lw=0.7))
        ax.text(x, 3.575, label, ha="center", va="center", fontsize=7.2, fontweight="bold")
    for i in range(len(main)-1):
        ax.annotate("", xy=(main[i+1][0]-w/2-0.06, 3.575), xytext=(main[i][0]+w/2+0.06, 3.575),
                    arrowprops=dict(arrowstyle="-|>", color="#303030", lw=1.0))

    # Row 2: evidence/boundary layer
    row2 = [
        (2.2,  "Technical variance\ndecomposition\n(conc < lane < rep)", "#FDE7E7"),
        (6.0,  "Identifiability boundary\n(low-\u03c3\u00b2 targets not estimable)", "#E7F0FA"),
        (9.8,  "Six-cohort falsification\n(684 cases, no valid biomarker)", "#F1F8E9"),
    ]
    for x, label, fc in row2:
        ax.add_patch(FancyBboxPatch((x-1.55, 1.0), 3.1, 1.05,
                     boxstyle="round,pad=0.02,rounding_size=0.12",
                     fc=fc, ec="#303030", lw=0.7))
        ax.text(x, 1.525, label, ha="center", va="center", fontsize=7.0, fontweight="bold")

    # vertical connectors
    ax.annotate("", xy=(2.2, 2.05), xytext=(5.85, 3.1),
                arrowprops=dict(arrowstyle="-|>", color="#888888", lw=0.8, ls=":"))
    ax.annotate("", xy=(6.0, 2.05), xytext=(8.15, 3.1),
                arrowprops=dict(arrowstyle="-|>", color="#888888", lw=0.8, ls=":"))
    ax.annotate("", xy=(9.8, 2.05), xytext=(10.45, 3.1),
                arrowprops=dict(arrowstyle="-|>", color="#888888", lw=0.8, ls=":"))

    ax.text(6.0, 4.05, "MAIN LINE: cross-patient and cross-cancer RNA\u2192ADT transportability",
            ha="center", fontsize=9.5, fontweight="bold", color=BLUE)
    ax.text(6.0, 2.65, "EVIDENCE & BOUNDARY LAYER",
            ha="center", fontsize=8, fontweight="bold", color="#555555")
    save(fig, "Figure6_model_framework", 11, 5.2)

# ============================================================
# Figure 7 : technical-noise hierarchical decomposition
# ============================================================
def fig7():
    # Panel A data: median frac of technical layers across datasets
    layers = [
        ("Staining\nconcentration\n(GSE245108)", 0.0069, BLUE),
        ("Sequencing\nlane\n(Kotliarov 2020)",    0.0126, GREEN),
        ("Inter-experiment\nplatform\n(GSE135325)",0.0017, ORANGE),
        ("Library prep\nreplicate\n(GSE245108)",  0.1848, RED),
    ]
    names = [l[0] for l in layers]
    vals  = [l[1] for l in layers]
    cols  = [l[2] for l in layers]

    # donor bimodal bins
    bins = [("<0.01",113), ("0.01-0.05",10), ("0.05-0.1",1), ("0.1-0.2",3), ("0.2-0.3",44), ("0.3-0.4",101)]
    bin_labels = [b[0] for b in bins]
    bin_counts = [b[1] for b in bins]

    # rep vs conc scatter from full csv
    csv_path = r"D:\ai_multimodal_cholesterol_study_outputs\revision-roadmap\technical_biological_decomposition_data_v1_20260903\GSE245108_tnc_parsed\tnc_variance_components_all.csv"
    rep_frac, conc_frac, conv = [], [], []
    with open(csv_path) as f:
        r = csv.DictReader(f)
        for row in r:
            if row.get("error") not in (None, "", "NaN"):
                continue
            try:
                rf = float(row["frac_rep"]); cf = float(row["frac_concentration"])
                cv = row.get("converged","") == "True"
            except (ValueError, KeyError):
                continue
            rep_frac.append(rf); conc_frac.append(cf); conv.append(cv)
    rep_frac = np.array(rep_frac); conc_frac = np.array(conc_frac); conv = np.array(conv)

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.6))

    # A: technical layers
    ax = axes[0]
    y = np.arange(len(layers))[::-1]
    ax.barh(y, vals, 0.6, color=cols, edgecolor="#303030", lw=0.6)
    for yi, v in zip(y, vals):
        ax.text(v + 0.004, yi, f"{v:.3f}", va="center", fontsize=8, fontweight="bold")
    ax.set_yticks(y); ax.set_yticklabels(names, fontsize=7.2)
    ax.set_xlim(0, 0.26)
    ax.set_xlabel("Median variance fraction", fontsize=8)
    ax.set_title("Technical noise layers\n(library-prep replicate dominates)", fontsize=8.5)
    ax.spines[["top", "right"]].set_visible(False)

    # B: donor bimodal
    ax = axes[1]
    ax.bar(range(len(bins)), bin_counts, 0.7, color=BLUE, edgecolor="#303030", lw=0.6)
    for i, c in enumerate(bin_counts):
        ax.text(i, c + 2, str(c), ha="center", fontsize=7.5)
    ax.set_xticks(range(len(bins))); ax.set_xticklabels(bin_labels, fontsize=7, rotation=30)
    ax.set_ylabel("Number of targets", fontsize=8)
    ax.set_title("GSE245108 donor variance is bimodal\n(113 \u22480 lineage vs 101 \u22480.33 donor-specific)", fontsize=8.5)
    ax.spines[["top", "right"]].set_visible(False)

    # C: rep vs conc scatter
    ax = axes[2]
    ax.scatter(conc_frac[~conv], rep_frac[~conv], s=10, color="#cccccc", label="non-converged", alpha=0.6, zorder=1)
    ax.scatter(conc_frac[conv], rep_frac[conv], s=12, color=BLUE, alpha=0.55, label="converged", zorder=2)
    ax.plot([0, 0.55], [0, 0.55], ls="--", color=GREY, lw=0.8)
    ax.text(0.42, 0.34, "rep = conc", rotation=45, fontsize=7, color=GREY)
    ax.set_xlabel("frac_concentration", fontsize=8)
    ax.set_ylabel("frac_rep", fontsize=8)
    ax.set_title("Per-target: rep >> concentration\n(most points below diagonal)", fontsize=8.5)
    ax.set_xlim(0, 0.55); ax.set_ylim(0, 0.6)
    ax.legend(fontsize=6.5, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    save(fig, "Figure7_technical_noise_decomposition", 11.5, 3.6)

# ============================================================
# Figure 8 : identifiability simulation
# ============================================================
def fig8():
    designs = ["D_confound\n(Kotliarov-like)", "D_cross\n(GSE245108-like)", "D_ideal\n(novo-like)"]
    low  = [0.213, 0.725, 0.404]
    high = [0.101, 0.185, 0.036]

    power_n  = [4, 8, 16, 32]
    power_rmse = [0.261, 0.155, 0.135, 0.113]
    power_rel  = [0.773, 0.459, 0.408, 0.338]

    # CLR bias: low/high donor x raw/cpm/clr, donor & rep
    clr_donor = {"low": [0.493, 0.006, 0.010], "high": [0.012, -0.315, -0.311]}
    clr_rep   = {"low": [0.184, -0.178, -0.175], "high": [0.011, -0.179, -0.176]}

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8))

    # A: donor bias by design
    ax = axes[0]
    x = np.arange(3); wd = 0.36
    ax.bar(x - wd/2, low, wd, color=RED, label="low-\u03c3\u00b2 donor (failure-like)", edgecolor="#303030", lw=0.6)
    ax.bar(x + wd/2, high, wd, color=BLUE, label="high-\u03c3\u00b2 donor", edgecolor="#303030", lw=0.6)
    ax.axhline(0, color="#303030", lw=0.7)
    for xi, (l, h) in enumerate(zip(low, high)):
        ax.text(xi-wd/2, l+0.02, f"{l:.2f}", ha="center", fontsize=7)
        ax.text(xi+wd/2, h+0.02, f"{h:.2f}", ha="center", fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels(designs, fontsize=7)
    ax.set_ylabel("Donor variance bias", fontsize=8)
    ax.set_title("Identifiability: low-\u03c3\u00b2 donor\nvariance is unestimable (large +bias)", fontsize=8.5)
    ax.legend(fontsize=6.3, frameon=False, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)

    # B: power curve
    ax = axes[1]
    ax.plot(power_n, power_rmse, "o-", color=BLUE, lw=1.2, ms=5)
    for n, r in zip(power_n, power_rmse):
        ax.text(n, r+0.008, f"{r:.3f}", ha="center", fontsize=7)
    ax.set_xlabel("Number of donors", fontsize=8)
    ax.set_ylabel("Donor RMSE", fontsize=8, color=BLUE)
    ax.set_title("Power: \u226516-32 donors needed\nfor precise donor-variance split", fontsize=8.5)
    ax.set_xticks(power_n)
    ax2 = ax.twinx()
    ax2.plot(power_n, power_rel, "s--", color=RED, lw=1.0, ms=4)
    ax2.set_ylabel("Relative SD", fontsize=8, color=RED)
    ax2.set_ylim(0, 0.9)
    ax.spines[["top"]].set_visible(False); ax2.spines[["top"]].set_visible(False)

    # C: CLR/CPM normalization bias (donor & rep, high_donor focus)
    ax = axes[2]
    x = np.arange(3); wd = 0.3
    ax.bar(x - wd*1.05, clr_donor["high"], wd, color=BLUE, label="donor bias", edgecolor="#303030", lw=0.5)
    ax.bar(x + wd*0.05, clr_rep["high"], wd, color=ORANGE, label="rep bias", edgecolor="#303030", lw=0.5)
    ax.axhline(0, color="#303030", lw=0.7)
    for xi, (d, r) in enumerate(zip(clr_donor["high"], clr_rep["high"])):
        ax.text(xi-wd*1.05, d + (0.02 if d>=0 else -0.05), f"{d:.2f}", ha="center", fontsize=6.8)
        ax.text(xi+wd*0.05, r + (0.02 if r>=0 else -0.05), f"{r:.2f}", ha="center", fontsize=6.8)
    ax.set_xticks(x); ax.set_xticklabels(["raw", "CPM", "CLR"], fontsize=8)
    ax.set_ylabel("Bias", fontsize=8)
    ax.set_title("CLR/CPM removes global rep effect\n(rep underestimated \u2248 -0.18)", fontsize=8.5)
    ax.legend(fontsize=6.3, frameon=False, loc="lower left")
    ax.set_ylim(-0.45, 0.25)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    save(fig, "Figure8_identifiability_simulation", 11.5, 3.8)

# ============================================================
# Figure 9 : PCA eigenvalue spectrum
# ============================================================
def fig9():
    ev = np.array([51.854, 3.174, 2.71, 1.639, 1.478, 1.322, 1.314, 1.228, 1.199, 1.161, 1.142, 1.132])
    Neff = 24.18
    # cumulative variance over full spectrum (approximate using reported k50/k80/k90)
    k50, k80, k90 = 72, 156, 186
    idx = np.arange(1, len(ev)+1)

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.4))

    # A: scree (log)
    ax = axes[0]
    ax.plot(idx, ev, "o-", color=BLUE, lw=1.2, ms=5)
    ax.axhline(1.0, ls="--", color=GREY, lw=0.8)
    ax.text(1, 1.0, " Kaiser=1", va="bottom", fontsize=6.5, color=GREY)
    ax.set_yscale("log")
    ax.set_xticks(idx)
    ax.set_xlabel("Principal component", fontsize=8)
    ax.set_ylabel("Eigenvalue (log)", fontsize=8)
    ax.set_title("Scree: \u03bb\u2081=51.85 dominant\n(parallel-analysis K=9)", fontsize=8.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.annotate("\u03bb\u2081 = 51.85 (19.6%)", xy=(1, 51.85), xytext=(3.5, 20),
                fontsize=7.5, color=RED, fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=RED, lw=0.8))

    # B: effective dimensionality
    ax = axes[1]
    ax.bar([0, 1], [265, Neff], 0.5, color=[GREY, RED], edgecolor="#303030", lw=0.6)
    ax.text(0, 265+5, "265", ha="center", fontsize=8, fontweight="bold")
    ax.text(1, Neff+5, f"{Neff:.1f}", ha="center", fontsize=8, fontweight="bold", color=RED)
    ax.text(1, Neff/2, "9.1%", ha="center", fontsize=8, color="white", fontweight="bold")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Antibodies\n(measured)", "Effective\nindependent"], fontsize=7.5)
    ax.set_ylabel("Dimensionality", fontsize=8)
    ax.set_title(f"Effective independent targets: Neff={Neff:.1f}\n(\u2248 9 significant factors)", fontsize=8.5)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    save(fig, "Figure9_pca_spectrum", 8.6, 3.4)

if __name__ == "__main__":
    fig4()
    fig6()
    fig7()
    fig8()
    fig9()
    print("ALL FIGURES REGENERATED")
