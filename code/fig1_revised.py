"""
Revised Figure 1 (A-E). Numbers from fig1_revised_numbers.py. Panels C and D show only the ratios by age (lung cancer death rate ratio; BRCA odds ratio);
the median ages are reported in panel E.
Output: figures/Figure1_revised.png and .pdf
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from scipy import stats

from paths import OUT, FIG
N = json.loads((OUT / "fig1_revised_numbers.json").read_text())
st = pd.read_csv(OUT / "nmsc" / "state_level_results.csv")

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
plt.rcParams.update({"font.size": 7.5, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2})


def letter(ax, s):
    ax.text(-0.16, 1.06, s, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")


fig = plt.figure(figsize=(7.2, 9.0))
gs = GridSpec(3, 2, figure=fig, height_ratios=[1, 1.05, 0.62], hspace=0.42, wspace=0.28)

# ---------------- A and B: NMSC vs UV ----------------
for k, (pos, col, ylab, key) in enumerate(((gs[0, 0], "asr", "Age-standardised incidence rate\n(cases per 100,000)", "A"),
                                           (gs[0, 1], "meda_cases", "Median age of diagnosed cases (years)", "B"))):
    ax = fig.add_subplot(pos); letter(ax, key)
    sc = ax.scatter(st.uv, st[col], c=st.lat, cmap="inferno_r", vmin=25, vmax=48, s=16, edgecolor="white", linewidth=0.4)
    r = stats.linregress(st.uv, st[col]); xs = np.linspace(st.uv.min(), st.uv.max(), 50)
    ax.plot(xs, r.intercept + r.slope * xs, color=INK, lw=1)
    a = N[key]
    unit, nd = ("cases/100,000", 0) if key == "A" else ("y", 2)
    ax.text(0.03, 0.97, f"R² = {a['R2']:.3f}, p = {a['p']:.2g}, n = {a['n']}\n"
            f"slope per 1,000 Wh/m²: {a['slope_per_1000']:.{nd}f} {unit}\n(95% CI {a['slope_lo']:.{nd}f} to {a['slope_hi']:.{nd}f})",
            transform=ax.transAxes, va="top", fontsize=6.8, color=INK)
    ax.set_xlabel("UV exposure (Wh/m², mean of counties 2020–2024)"); ax.set_ylabel(ylab)
    if key == "A":
        ax.set_ylim(700, 1250)
    if key == "B":
        ax.set_ylim(69.5, 75.5)
        cb = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02); cb.set_label("latitude (°N)", fontsize=6.5); cb.ax.tick_params(labelsize=6)
    ax.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)

# ---------------- C: lung ----------------
axC = fig.add_subplot(gs[1, 0]); letter(axC, "C")
c, el = N["C"], N["E_lung"]
rr = np.array([r[1:] for r in c["rr_by_age"]], float); lab = [r[0] for r in c["rr_by_age"]]
xi = np.arange(len(rr))
axC.axhspan(el["RR_lo"], el["RR_hi"], color=GRID, lw=0)
axC.axhline(el["RR"], color=INK2, lw=0.9, ls="--")
axC.errorbar(xi, rr[:, 0], yerr=[rr[:, 0] - rr[:, 1], rr[:, 2] - rr[:, 0]], fmt="o", color=ORANGE, ms=4, elinewidth=1, capsize=0)
axC.set_yscale("log"); axC.set_ylim(4, 60); axC.set_yticks([5, 10, 20, 50]); axC.set_yticklabels(["5", "10", "20", "50"])
axC.minorticks_off()
axC.set_xticks(xi); axC.set_xticklabels(lab, fontsize=6.3); axC.set_xlim(-0.5, len(rr) - 0.5)
axC.set_xlabel("Age at lung cancer death (years)")
axC.set_ylabel("Lung cancer death rate ratio\n(current / never smokers)")
axC.text(0.03, 0.97, f"age-adjusted rate ratio {el['RR']:.1f} ({el['RR_lo']:.1f}–{el['RR_hi']:.1f}), dashed line\n"
         f"change per decade {el['trend']:.2f} ({el['trend_lo']:.2f}–{el['trend_hi']:.2f})",
         transform=axC.transAxes, va="top", fontsize=6.8, color=INK)
axC.yaxis.grid(True, color=GRID, lw=0.5); axC.set_axisbelow(True)

# ---------------- D: BRCA ----------------
axD = fig.add_subplot(gs[1, 1]); letter(axD, "D")
d = N["D"]
e1, e2 = N["E_BRCA1"], N["E_BRCA2"]
for g, col, off in (("BRCA1", BLUE, -0.08), ("BRCA2", ORANGE, 0.08)):
    o, ci = d["orplus"][g]; o = np.array(o); ci = np.array(ci)
    axD.errorbar(np.arange(4) + off, o, yerr=[o - ci[:, 0], ci[:, 1] - o], fmt="o-", color=col, ms=4, lw=1.2,
                 elinewidth=1, label=f"{g} carriers")
axD.axhline(1, color=INK2, lw=0.8, ls=":")
axD.set_yscale("log"); axD.set_ylim(0.8, 150); axD.set_yticks([1, 2, 5, 10, 20, 50, 100])
axD.set_yticklabels(["1", "2", "5", "10", "20", "50", "100"]); axD.minorticks_off()
axD.set_xticks(range(4)); axD.set_xticklabels(["<40", "40–49", "50–59", "≥60"]); axD.set_xlim(-0.5, 3.5)
axD.set_xlabel("Age at breast cancer diagnosis (years)")
axD.set_ylabel("Breast cancer odds ratio\n(carriers / non-carriers)")
axD.legend(loc="lower left", frameon=False, fontsize=6.8, bbox_to_anchor=(0.0, 0.04))
axD.text(0.03, 0.97, f"case-only OR per year of age\nBRCA1 {e1['case_only_or']:.3f} ({e1['case_only_lo']:.2f}–{e1['case_only_hi']:.2f})\n"
         f"BRCA2 {e2['case_only_or']:.3f} ({e2['case_only_lo']:.2f}–{e2['case_only_hi']:.2f})",
         transform=axD.transAxes, va="top", fontsize=6.8, color=INK)
axD.yaxis.grid(True, color=GRID, lw=0.5); axD.set_axisbelow(True)

# ---------------- E: summary table ----------------
axE = fig.add_subplot(gs[2, :]); axE.axis("off"); letter(axE, "E")
f1 = lambda v, lo, hi, nd=1: f"{v:+.{nd}f} ({lo:+.{nd}f} to {hi:+.{nd}f})"
en, el = N["E_NMSC"], N["E_lung"]
rows = [
    ["NMSC\n(high UV / low UV)", f"rate ratio {en['RR']:.2f}\n({en['RR_lo']:.2f}–{en['RR_hi']:.2f})",
     f"{en['meda_low']:.1f}", f"{en['meda_high']:.1f}", f1(en['dmeda'], en['dmeda_lo'], en['dmeda_hi']), f"{en['trend']:.2f}ᵃ"],
    ["Lung cancer death†\n(current / never smoker)", f"rate ratio {el['RR']:.1f}\n({el['RR_lo']:.1f}–{el['RR_hi']:.1f})",
     f"{el['meda_low']:.1f}", f"{el['meda_high']:.1f}", f1(el['dmeda'], el['dmeda_lo'], el['dmeda_hi']),
     f"{el['trend']:.2f} ({el['trend_lo']:.2f}–{el['trend_hi']:.2f})"],
    ["Breast cancer\n(BRCA1 / non-carriers)", f"odds ratio {e1['OR_young']:.0f} (<40 y)\n→ {e1['OR_old']:.1f} (≥60 y)",
     "n.p.ᶜ", f"{e1['meda_high']:.1f}\n({e1['meda_high_lo']:.1f}–{e1['meda_high_hi']:.1f})", "earlierᶜ",
     f"{e1['trend']:.2f} ({e1['trend_lo']:.2f}–{e1['trend_hi']:.2f})ᵇ"],
    ["Breast cancer\n(BRCA2 / non-carriers)", f"odds ratio {e2['OR_young']:.0f} (<40 y)\n→ {e2['OR_old']:.1f} (≥60 y)",
     "n.p.ᶜ", f"{e2['meda_high']:.1f}\n({e2['meda_high_lo']:.1f}–{e2['meda_high_hi']:.1f})", "earlierᶜ",
     f"{e2['trend']:.2f} ({e2['trend_lo']:.2f}–{e2['trend_hi']:.2f})ᵇ"],
]
cols = ["Comparison\n(higher / lower risk)", "Increase in risk\n(95% CI)", "MEDA (y)\nlower risk",
        "MEDA (y)\nhigher risk", "Difference in\nMEDA, y (95% CI)", "Relative risk change\nper decade (95% CI)"]
tab = axE.table(cellText=rows, colLabels=cols, loc="upper center", cellLoc="center",
                colWidths=[0.22, 0.19, 0.10, 0.12, 0.17, 0.20], bbox=[-0.06, 0.0, 1.10, 1])
tab.auto_set_font_size(False); tab.set_fontsize(6.2)
for (r_, c_), cell in tab.get_celld().items():
    cell.set_edgecolor("#bdbcb6"); cell.set_linewidth(0.5)
    if r_ == 0:
        cell.set_facecolor("#5c5b57"); cell.set_text_props(color="white", fontweight="bold")
    elif c_ == 0:
        cell.set_facecolor("#efeee9"); cell.set_text_props(fontweight="bold")

fig.savefig(FIG / "Figure1_revised.png", dpi=220, bbox_inches="tight")
fig.savefig(FIG / "Figure1_revised.pdf", bbox_inches="tight")
print("saved", FIG / "Figure1_revised.png")
