#!/usr/bin/env python
"""Figures of the revision analysis. Reads the tables and results of steps 03 to 16.

fig1_antibody_titres               GMT of binding and neutralising antibodies by sex, arm and day
fig2_antibody_sex_difference       male / female GMR at V3D28 (day 252) by model and arm, and at each
                                   study day from the longitudinal mixed model
fig3_signature_response            D1 vs D0 change of the vaccine-response and ISG signatures per
                                   participant, and the female minus male difference by model
fig4_gene_level_response           female vs male log2 fold change per gene, per-sex DEG counts at
                                   full and equal group size, mean female minus male change over
                                   responsive genes against the permutation null
fig5_irf7_isg15                    IRF7 and ISG15 per participant at D0 and D1, sex difference in
                                   the fold change
fig6_wgcna_robustness              module-trait correlations with FDR, networks rebuilt without one
                                   participant
fig7_variance_and_cell_proportions variance partitioning, change in estimated cell proportions

Output: figures/<name>.pdf and .png
"""
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fig_style as fs
from paths import FIG, PART, RES, TAB

WGCNA = os.path.join(RES, "wgcna_robustness")

ARMS4 = ["HD SI", "HD LI", "LD SI", "LD LI"]
ASSAYS = ["Binding (ELISA S1)", "Neutralising (PRNT50)"]
ASSAY_COL = {ASSAYS[0]: "black", ASSAYS[1]: fs.POINT}
LLOQ = {ASSAYS[0]: 4.6, ASSAYS[1]: 16}
DAYS = [0, 14, 28, 42, 56, 70, 84, 168, 224, 252]
XP = {d: i for i, d in enumerate(DAYS)}
VACC = {"SI": [0, 28, 224], "LI": [0, 56, 224]}
GROUPS = [(a, v) for a in ["SI", "LI"] for v in ["T1", "T2", "T3"]]


def ptxt(p):
    return f"p = {p:.2f}" if p >= 0.1 else (f"p = {p:.3f}" if p >= 0.001 else "p < 0.001")


def group_axis(ax, xpos, y=-0.21):
    ax.set_xticks([xpos[k] for k in GROUPS]); ax.set_xticklabels([v for _, v in GROUPS])
    for a in ["SI", "LI"]:
        x = np.mean([xpos[(a, v)] for v in ["T1", "T2", "T3"]])
        ax.text(x, y, a, transform=ax.get_xaxis_transform(), ha="center", fontweight="bold")


fs.setup()

# ------------------------------------------------------------------ 1 antibody titres over time
gmt = pd.read_csv(os.path.join(TAB, "antibody_gmt.csv"))
pop = pd.read_csv(os.path.join(TAB, "antibody_population.csv"))
fig, axes = plt.subplots(2, 4, figsize=(7.2, 4.4), sharey="row", gridspec_kw={"hspace": 0.55, "wspace": 0.1})
for r, assay in enumerate(ASSAYS):
    ylim = (2, 2000) if r == 0 else (8, 2000)
    for c, arm in enumerate(ARMS4):
        ax = axes[r, c]
        for j, sex in enumerate(fs.SEXES):
            d = gmt[(gmt.assay == assay) & (gmt.arm == arm) & (gmt.Sex == sex)].sort_values("day")
            x = np.array([XP[k] for k in d.day]) + (j - 0.5) * 0.18
            ax.errorbar(x, d.GMT, yerr=[d.GMT - d.CI_low, d.CI_high - d.GMT], color=fs.SEX_COL[sex],
                        lw=0.9, elinewidth=0.5, marker="o", ms=2.2, capsize=0)
        ax.axhline(LLOQ[assay], color=fs.GREY, lw=0.5, ls=":")
        for v in VACC[arm.split()[1]]:
            ax.plot(XP[v], ylim[0] * 1.25, marker="^", color="black", ms=3.2, clip_on=False, zorder=5)
        ax.set_yscale("log"); ax.set_ylim(*ylim); ax.minorticks_off()
        ax.set_xticks(range(len(DAYS)))
        ax.set_xticklabels([str(k) if k in (0, 56, 84, 252) else "" for k in DAYS], fontsize=5.5)
        n = pop[pop.arm == arm].set_index("sex").mITT
        ax.set_title(f"{arm}\n{n['Female']} female, {n['Male']} male", fontsize=6.5)
        if c == 0:
            ax.set_ylabel(f"{assay}\nGMT (IU/mL)")
            ax.text(-0.45, 1.18, "ab"[r], transform=ax.transAxes, fontsize=8, fontweight="bold")
        if r == 1:
            ax.set_xlabel("Study day", labelpad=1)
h = fs.sex_handles()
h.append(plt.Line2D([], [], ls="", marker="^", color="black", ms=3.5, label="Vaccination"))
h.append(plt.Line2D([], [], color=fs.GREY, lw=0.5, ls=":", label="Lower limit of the assay"))
fig.legend(handles=h, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.02))
fs.save(fig, FIG, "fig1_antibody_titres")

# ------------------------------------------------------------------ 2 antibody sex difference
C = pd.read_csv(os.path.join(TAB, "antibody_sex_difference_by_day.csv"))
D = pd.read_csv(os.path.join(TAB, "antibody_sex_difference_mixed_model.csv"))
E = pd.read_csv(os.path.join(TAB, "antibody_model_tests.csv"))
fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw={"width_ratios": [1, 1.15], "wspace": 0.55})
ax = axes[0]
rows = [("unadjusted", "all vaccinated", "Unadjusted"), ("arm", "all vaccinated", "Arm"),
        ("arm + age", "all vaccinated", "Arm and age"), ("arm + age + D0", "all vaccinated", "Arm, age, D0 titre"),
        ("arm + age, Tobit", "all vaccinated", "Arm and age, censored")] + [("per arm", a, f"{a} only") for a in ARMS4]
for y, (mn, arm, lab) in enumerate(rows):
    for k, assay in enumerate(ASSAYS):
        s = C[(C.assay == assay) & (C.day == 252) & (C.model == mn) & (C.arm == arm)].iloc[0]
        yy = y + (k - 0.5) * 0.34
        ax.plot([s.CI_low, s.CI_high], [yy, yy], color=ASSAY_COL[assay], lw=0.7)
        ax.plot(s.GMR_male_vs_female, yy, "o", color=ASSAY_COL[assay], ms=2.5)
ax.axhline(4.5, color="#D9D9D9", lw=0.5)
ax.axvline(1, color=fs.GREY, lw=0.5, ls=":")
ax.set_xscale("log"); ax.set_xlim(0.2, 12); ax.minorticks_off()
ax.set_xticks([0.25, 0.5, 1, 2, 4, 8]); ax.set_xticklabels(["0.25", "0.5", "1", "2", "4", "8"])
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[2] for r in rows], fontsize=6)
ax.set_ylim(len(rows) - 0.4, -1.7); ax.tick_params(axis="y", length=0)
ax.set_xlabel("Male / female GMR (95% CI)")
ax.set_title("a  V3D28 (day 252), models", loc="left")
ax.legend(handles=[plt.Line2D([], [], color=ASSAY_COL[a], marker="o", ms=2.5, lw=0.7, label=a.split(" (")[0])
                   for a in ASSAYS], loc="upper right", fontsize=5.5)

ax = axes[1]
for k, assay in enumerate(ASSAYS):
    d = D[(D.assay == assay) & (D.model == "all time points")].sort_values("day")
    x = np.array([XP[v] for v in d.day]) + (k - 0.5) * 0.22
    ax.errorbar(x, d.GMR_male_vs_female, yerr=[d.GMR_male_vs_female - d.CI_low, d.CI_high - d.GMR_male_vs_female],
                color=ASSAY_COL[assay], lw=0, elinewidth=0.7, marker="o", ms=2.5, label=assay.split(" (")[0])
    pi = E[(E.assay == assay) & (E.test == "sex x time (F test, Satterthwaite)")].pvalue.iloc[0]
    ph = d[d.day == 252].p_holm.iloc[0]
    ax.text(0.02, 0.97 - 0.08 * k, f"{assay.split(' (')[0]}: sex x time {ptxt(pi)}, day 252 Holm {ptxt(ph)}",
            transform=ax.transAxes, fontsize=5.5, va="top")
ax.axhline(1, color=fs.GREY, lw=0.5, ls=":")
ax.set_yscale("log"); ax.set_ylim(0.25, 8); ax.minorticks_off()
ax.set_yticks([0.25, 0.5, 1, 2, 4]); ax.set_yticklabels(["0.25", "0.5", "1", "2", "4"])
ax.set_xticks(range(len(DAYS))); ax.set_xticklabels(DAYS, fontsize=5.5)
ax.set_xlabel("Study day"); ax.set_ylabel("Male / female GMR (95% CI)")
ax.set_title("b  Each study day, mixed model", loc="left")
ax.legend(loc="lower right", fontsize=5.5)
fs.save(fig, FIG, "fig2_antibody_sex_difference")

# ------------------------------------------------------------------ 3 signature response per participant
sc = pd.read_csv(os.path.join(PART, "participant_signature_scores.csv"))
res = pd.read_csv(os.path.join(TAB, "signature_sex_difference.csv"))
mix = pd.read_csv(os.path.join(TAB, "signature_sex_difference_mixed_model.csv"))
mixr = mix[mix.contrast == "response"].assign(model="Mixed model (participant random effect, age)")
res = pd.concat([res, mixr[["score", "arm", "vaccination", "model", "difference_F_minus_M", "CI_low", "CI_high"]]])
NAMES = {"vaccine_response_signature": "Vaccine-response signature", "ISG_signature": "ISG signature"}
MCOL = {"Unadjusted": "black", "Adjusted for age": fs.GREY, "Adjusted for D0 level": fs.POINT,
        "Adjusted for age and cell proportions": fs.ACCENT, "Mixed model (participant random effect, age)": fs.ACCENT_2}
pos = {g: i + 0.5 * (g[0] == "LI") for i, g in enumerate(GROUPS)}
fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.8), gridspec_kw={"width_ratios": [1.5, 1], "wspace": 0.35, "hspace": 0.62})
rng = np.random.default_rng(3)
for row, score in enumerate(NAMES):
    ax = axes[row, 0]
    for g, x in pos.items():
        d = sc[(sc.arm == g[0]) & (sc.vaccination == g[1])]
        for k, sex in enumerate(fs.SEXES):
            y = d.loc[d.sex == sex, score].values
            xx = x + (k - 0.5) * 0.36
            ax.scatter(xx + rng.uniform(-0.06, 0.06, len(y)), y, s=9, color=fs.SEX_COL[sex], edgecolor="black", lw=0.3, zorder=3)
            ax.hlines(np.median(y), xx - 0.14, xx + 0.14, color="black", lw=0.9, zorder=4)
    ax.axhline(0, color=fs.GREY, lw=0.5, ls=":")
    group_axis(ax, pos)
    ax.set_ylabel(f"{NAMES[score]}\nchange D1 vs D0 (log2)")
    ax.set_title(f"{'ac'[row]}  {NAMES[score]}, each participant", loc="left")
    if row == 0:
        ax.legend(handles=fs.sex_handles(markeredgecolor="black", markeredgewidth=0.3), loc="upper right", ncol=2)
    ax = axes[row, 1]
    rr = res[res.score == score]
    for i, g in enumerate(GROUPS):
        yc = i + 0.5 * (g[0] == "LI")
        for k, lab in enumerate(MCOL):
            r_ = rr[(rr.arm == g[0]) & (rr.vaccination == g[1]) & (rr.model == lab)].iloc[0]
            y = yc + (k - 2) * 0.14
            ax.plot([r_.CI_low, r_.CI_high], [y, y], color=MCOL[lab], lw=0.8)
            ax.plot(r_.difference_F_minus_M, y, "o", color=MCOL[lab], ms=2.3)
    ax.axvline(0, color=fs.GREY, lw=0.5, ls=":")
    ax.set_yticks([i + 0.5 * (g[0] == "LI") for i, g in enumerate(GROUPS)])
    ax.set_yticklabels([f"{a} {v}" for a, v in GROUPS]); ax.set_ylim(6.0, -0.6)
    ax.set_xlabel("Female minus male change (95% CI)")
    ax.set_title(f"{'bd'[row]}  Female minus male", loc="left")
fig.legend(handles=[plt.Line2D([], [], color=c, marker="o", ms=2.5, lw=0.8, label=l) for l, c in MCOL.items()],
           loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.06), fontsize=6)
fs.save(fig, FIG, "fig3_signature_response")

# ------------------------------------------------------------------ 4 gene-level response
gl = pd.read_csv(os.path.join(RES, "gene_level_deseq2", "gene_level_results.csv.gz"))
gl = gl[gl.model == "basic"]
A = pd.read_csv(os.path.join(RES, "deg_counts_permutation", "deg_counts_equal_n_all_subsets.csv"))
P = pd.read_csv(os.path.join(TAB, "sex_difference_permutation_null.csv"))
fig = plt.figure(figsize=(7.2, 4.2))
gs = fig.add_gridspec(2, 6, height_ratios=[0.8, 1.15], hspace=0.55, wspace=0.5)
gsb = gs[1, :].subgridspec(1, 2, wspace=0.35)
for c, (arm, vac) in enumerate(GROUPS):
    ax = fig.add_subplot(gs[0, c])
    b = gl[(gl.arm == arm) & (gl.vaccination == vac)]
    w = b.pivot(index="gene", columns="estimate", values="log2FC")
    q = b.pivot(index="gene", columns="estimate", values="padj")
    resp = q.index[(q.female_D1_vs_D0 < 0.05) | (q.male_D1_vs_D0 < 0.05)]
    sig = q.index[q.sex_difference_F_minus_M < 0.05].intersection(resp)
    ax.scatter(w.loc[resp, "male_D1_vs_D0"], w.loc[resp, "female_D1_vs_D0"], s=1, color="#BDBDBD", lw=0, rasterized=True)
    ax.scatter(w.loc[sig, "male_D1_vs_D0"], w.loc[sig, "female_D1_vs_D0"], s=5, color=fs.ACCENT, lw=0)
    ax.plot([-6, 6], [-6, 6], color="black", lw=0.5, ls="--")
    ax.set_xlim(-6, 6); ax.set_ylim(-6, 6); ax.set_aspect("equal")
    ax.set_xticks([-4, 0, 4]); ax.set_yticks([-4, 0, 4])
    r = np.corrcoef(w.loc[resp, "male_D1_vs_D0"], w.loc[resp, "female_D1_vs_D0"])[0, 1]
    ax.text(-5.6, 5.4, f"r = {r:.2f}", fontsize=5, va="top")
    ax.set_title(f"{arm} {vac}", fontsize=7, fontweight="bold")
    ax.set_xlabel("Male log2FC", fontsize=6)
    if c == 0:
        ax.set_ylabel("Female log2FC", fontsize=6)
        ax.text(-0.9, 1.35, "a  Change D1 vs D0 per gene, female against male", transform=ax.transAxes, fontsize=8)

xpos = {g: i + 0.6 * (g[0] == "LI") for i, g in enumerate(GROUPS)}
ax = fig.add_subplot(gsb[0])
for arm, vac in GROUPS:
    for j, sex in enumerate(fs.SEXES):
        x = xpos[(arm, vac)] + (j - 0.5) * 0.36
        d = A[(A.arm == arm) & (A.vaccination == vac) & (A.sex == sex)]
        full = d[d["sample"] == "full"].iloc[0]
        eq = d[d["sample"] == "equal n"]
        if len(eq):
            ax.scatter(x + rng.uniform(-0.07, 0.07, len(eq)), eq.n_DEG, s=2, color=fs.SEX_COL[sex], alpha=0.6, lw=0)
            ax.hlines(np.median(eq.n_DEG), x - 0.13, x + 0.13, color="black", lw=0.8, zorder=4)
        ax.scatter(x, full.n_DEG, s=20, marker="D", color=fs.SEX_COL[sex], edgecolor="black", lw=0.4, zorder=5)
        ax.text(x, -0.08, f"{full.n}", transform=ax.get_xaxis_transform(), ha="center", fontsize=5)
group_axis(ax, xpos, y=-0.28)
ax.tick_params(axis="x", pad=9)
ax.set_ylabel("Genes changed D1 vs D0\n(FDR < 0.05, each sex separately)")
ax.set_title("b  Genes changed in each sex", loc="left")
h = fs.sex_handles(markeredgecolor="black", markeredgewidth=0.3)
h += [plt.Line2D([], [], marker="D", color="white", markerfacecolor="white", markeredgecolor="black", ms=4, label="All participants"),
      plt.Line2D([], [], marker="o", color="white", markerfacecolor=fs.GREY, ms=2.5, label="Subsets, equal n per sex")]
ax.legend(handles=h, loc="upper left", ncol=2, fontsize=5)
ax.set_ylim(0, 9800)

PM = {"unadjusted": ("Unadjusted", "black"), "D0": ("Adjusted for D0 level of the gene", fs.POINT),
      "D0 + age": ("Adjusted for D0 level of the gene and age", fs.ACCENT_2)}
ax = fig.add_subplot(gsb[1])
for g in GROUPS:
    for j, mn in enumerate(PM):
        r_ = P[(P.arm == g[0]) & (P.vaccination == g[1]) & (P.model == mn) & (P.measure == "direction")].iloc[0]
        x = xpos[g] + (j - 1) * 0.25
        ax.plot([x, x], [r_.null_median, r_.null_q95], color="#BDBDBD", lw=3.5, solid_capstyle="butt", zorder=1)
        ax.scatter(x, r_.observed, s=10, color=PM[mn][1], zorder=3, lw=0)
ax.axhline(0, color="black", lw=0.5, ls=":")
group_axis(ax, xpos)
ax.set_ylabel("Female minus male change,\nmean over responsive genes (log2)")
ax.set_title("c  Female minus male, all responsive genes", loc="left")
h = [plt.Line2D([], [], color=c, marker="o", lw=0, ms=3, label=l) for l, c in PM.values()]
h.append(plt.Line2D([], [], color="#BDBDBD", lw=3.5, label="Sex labels permuted, median to 95th percentile"))
ax.legend(handles=h, loc="lower left", fontsize=5)
lo, hi = ax.get_ylim(); ax.set_ylim(lo * 1.9, hi * 1.1)
fs.save(fig, FIG, "fig4_gene_level_response")

# ------------------------------------------------------------------ 5 IRF7 and ISG15 per participant
pp = pd.read_csv(os.path.join(PART, "irf7_isg15_per_participant.csv"))
st = pd.read_csv(os.path.join(TAB, "irf7_isg15_fold_change.csv"))
GENES, ARMS, VACS = ["IRF7", "ISG15"], ["SI", "LI"], ["T1", "T2", "T3"]
fig, axes = plt.subplots(2, 6, figsize=(7.2, 3.6), sharey="row",
                         gridspec_kw={"wspace": 0.15, "hspace": 0.75})
for r, g in enumerate(GENES):
    for c, (arm, vac) in enumerate([(a, v) for a in ARMS for v in VACS]):
        ax = axes[r, c]
        d = pp[(pp.gene == g) & (pp.arm == arm) & (pp.vaccination == vac)]
        for sex in fs.SEXES[::-1]:
            for _, x in d[d.sex == sex].iterrows():
                ax.plot([0, 1], [x.D0, x.D1], color=fs.SEX_COL[sex], lw=0.7, alpha=0.9,
                        marker="o", ms=2.2, mec="black", mew=0.2)
        s = st[(st.gene == g) & (st.arm == arm) & (st.vaccination == vac)].iloc[0]
        p = s["pvalue_sex_difference_F_minus_M"]
        dd = s['log2FC_sex_difference_F_minus_M']
        plab = f"p = {p:.2f}" if p >= 0.01 else (f"p = {p:.3f}" if p >= 0.001 else "p < 0.001")
        ph = s["p_holm_12_sex_difference_F_minus_M"]
        if ph < 0.05:
            plab += f"\nHolm p = {ph:.3f}"
        ax.text(0.5, 1.0, f"F minus M = {dd:.2f}\n{plab}",
                transform=ax.transAxes, ha="center", va="bottom", fontsize=5.5)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["D0", "D1"]); ax.set_xlim(-0.3, 1.3)
        ax.set_xlabel(f"{arm} {vac}", fontweight="bold")
        if c == 0:
            ax.set_ylabel(f"{g}\nlog2(TPM+1)")
        if c == 3:
            ax.spines["left"].set_position(("outward", 0))
fig.legend(handles=[plt.Line2D([], [], color=c, lw=1, marker="o", ms=2.5, mec="black", mew=0.2, label=s)
                    for s, c in fs.SEX_COL.items()], loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.04))
fs.save(fig, FIG, "fig5_irf7_isg15")

# ------------------------------------------------------------------ 6 WGCNA robustness
mt = pd.read_csv(os.path.join(WGCNA, "module_trait_with_FDR.csv"))
num = lambda s: int(re.sub(r"\D", "", s))  # noqa: E731
HIGHLIGHT = {"SI": "M19", "LI": "M23"}
LABEL_Y = {"SI": (0.95, 0.62), "LI": (0.78, 0.55)}   # label block, fraction of max -log10 p
COUNT_Y = {"SI": 0.40, "LI": 0.22}                   # cell counts, axes fraction
TMAP = {"SI": {"V1": "T1", "V2": "T2", "V4": "T3"}, "LI": {"V1": "T1", "V3": "T2", "V4": "T3"}}
fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.8), gridspec_kw={"width_ratios": [1, 1.5], "wspace": 0.3, "hspace": 0.55})
for r, arm in enumerate(["SI", "LI"]):
    ax = axes[r, 0]
    t = mt[mt.arm == arm].sort_values("p").reset_index(drop=True)
    lp = -np.log10(t.p); sig = t.p_BH_matrix < 0.05
    ax.scatter(np.arange(len(t))[~sig] + 1, lp[~sig], s=4, color=fs.GREY, lw=0)
    ax.scatter(np.arange(len(t))[sig] + 1, lp[sig], s=9, color=fs.ACCENT, lw=0)
    for i, x in t[sig].iterrows():
        v, d_, sx = re.match(r"(V\d)(D\d)_(\w+)", x.original_trait).groups()
        ax.annotate(f"{x.module}, {TMAP[arm][v]} {d_} {sx.lower()}", (i + 1, lp[i]), xytext=(6, 0),
                    textcoords="offset points", fontsize=5.5, va="center")
    # vaccine-response module of the paper (SI M19, LI M23): cells with p < 0.05 ringed and labelled with FDR
    mod = HIGHLIGHT[arm]
    cells = t[(t.module == mod) & (t.p < 0.05)]
    ax.scatter(cells.index + 1, lp[cells.index], s=16, facecolor="none", edgecolor="black", lw=0.6, zorder=4)
    ys = np.linspace(*[f * lp.max() for f in LABEL_Y[arm]], len(cells))
    for k, (i, x) in enumerate(cells.iterrows()):
        v, d_, sx = re.match(r"(V\d)(D\d)_(\w+)", x.original_trait).groups()
        ax.annotate(f"{mod}, {TMAP[arm][v]} {d_} {sx.lower()}, FDR {x.p_BH_matrix:.2f}", (i + 1, lp[i]),
                    xytext=(len(t) * 0.3, ys[k]), textcoords="data", fontsize=5, va="center",
                    arrowprops=dict(arrowstyle="-", color=fs.GREY, lw=0.3, relpos=(0, 0.5)))
    ax.axhline(-np.log10(0.05), color="black", lw=0.5, ls=":")
    ax.text(len(t), -np.log10(0.05) - 0.03 * lp.max(), "p = 0.05", ha="right", va="top", fontsize=5.5)
    ax.set_xlabel("Module-trait cells, ranked by p"); ax.set_ylabel("-log10 p")
    ax.set_title(f"{'ac'[r]}  {arm}, module-trait correlations", loc="left")
    ax.text(0.97, COUNT_Y[arm], f"{len(t)} cells\n{int((t.p < 0.05).sum())} with p < 0.05\n{int(sig.sum())} with FDR < 0.05",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=5.5)
    ax = axes[r, 1]
    per = pd.read_csv(os.path.join(WGCNA, f"loso_jaccard_per_subject_{arm}.csv")).set_index("module")
    per = per.loc[sorted(per.index, key=num)]
    rng4 = np.random.default_rng(4)
    for i, (mod, vals) in enumerate(per.iterrows()):
        col = fs.ACCENT if mod == HIGHLIGHT[arm] else fs.POINT
        ax.scatter(i + rng4.uniform(-0.15, 0.15, len(vals)), vals.values, s=3, color=col, lw=0, alpha=0.8)
        ax.hlines(np.median(vals.values), i - 0.35, i + 0.35, color="black", lw=0.8)
    ax.axhline(0.5, color="black", lw=0.5, ls="--")
    ax.set_xticks(range(len(per))); ax.set_xticklabels(per.index, rotation=90, fontsize=5)
    ax.set_ylim(0, 1.02); ax.set_ylabel("Jaccard index with\nbest-matching module")
    ax.set_title(f"{'bd'[r]}  {arm}, networks rebuilt without one participant", loc="left")
fs.save(fig, FIG, "fig6_wgcna_robustness")

# ------------------------------------------------------------------ 7 variance partitioning and cell proportions
var = pd.read_csv(os.path.join(TAB, "variance_partitioning.csv"))
cells = pd.read_csv(os.path.join(PART, "cell_proportion_change.csv"))
CT = ["Neutrophils", "Monocytes", "Lymphocytes"]
fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.6), gridspec_kw={"width_ratios": [1.2, 1], "wspace": 0.45})
ax = axes[0]
y = np.arange(len(var))[::-1]
ax.barh(y + 0.19, var.all_genes_pct, height=0.38, color=fs.POINT, label="All expressed genes")
ax.barh(y - 0.19, var.PC1_pct, height=0.38, color=fs.LIGHT, label="First principal component")
ax.set_yticks(y); ax.set_yticklabels(var.factor, fontsize=6)
ax.set_xlabel("Variance explained (%)")
ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.55))
ax.set_title("a  Sources of variation", loc="left")
ax = axes[1]
chs = cells.groupby(["Subject", "Sex"])[CT].mean().reset_index()
rng2 = np.random.default_rng(2)
for i, ct in enumerate(CT):
    for k, sex in enumerate(fs.SEXES):
        v = chs.loc[chs.Sex == sex, ct].values
        xx = i + (k - 0.5) * 0.34
        ax.scatter(xx + rng2.uniform(-0.06, 0.06, len(v)), v, s=9, color=fs.SEX_COL[sex], edgecolor="black", lw=0.3, zorder=3)
        ax.hlines(np.median(v), xx - 0.13, xx + 0.13, color="black", lw=0.9, zorder=4)
ax.axhline(0, color=fs.GREY, lw=0.5, ls=":")
ax.set_xticks(range(len(CT))); ax.set_xticklabels(CT)
ax.set_ylabel("Change D1 vs D0\n(percentage points)")
ax.legend(handles=fs.sex_handles(markeredgecolor="black", markeredgewidth=0.3), loc="upper right")
ax.set_title("b  Estimated cell proportions", loc="left")
fs.save(fig, FIG, "fig7_variance_and_cell_proportions")
