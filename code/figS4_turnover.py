"""
Supplementary Figure S4: matched-input turnover (Comment 9). Constant input, no aging; mutation input per cell per
generation (CTOR x events per division / genome factor) held at 0.03 or 0.10 across CTOR 0.01, 0.02 and 0.2.
A, B: cumulative occurrence (95% Wilson bands; dots, MEDA). C: fate of every run. D: observed - PH-predicted MEDA,
CTOR 0.02 vs 0.01 (output/turnover/exp2_matched_comparisons.csv; run phase2_analysis.py exp2 first).
Output: figures/FigureS4_turnover.png and .pdf
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

from paths import SIM, TURNOVER, OUT, FIG
from ph_benchmark import wilson
from phase2_analysis import load

T_END = 92.0
REF = SIM / "revision" / "constant_events_escc_no_aging"
RUNS = {  # matched input -> [(CTOR, label, file)]
    0.03: [(0.01, "CTOR 0.01", REF / "constant_3_automatics_ancestralResults.csv"),
           (0.02, "CTOR 0.02, genome ×2", TURNOVER / "ctor002_e3_g2_automatics_ancestralResults.csv"),
           (0.2, "CTOR 0.2, genome ×20", TURNOVER / "ctor_02_e3_g20_automatics_ancestralResults.csv")],
    0.10: [(0.01, "CTOR 0.01", REF / "constant_10_automatics_ancestralResults.csv"),
           (0.02, "CTOR 0.02, genome ×2", TURNOVER / "ctor002_e10_g2_automatics_ancestralResults.csv"),
           (0.2, "CTOR 0.2, genome ×20", TURNOVER / "ctor_02_e10_g20_automatics_ancestralResults.csv")],
}
COL = {0.01: "#3b0f70", 0.02: "#de4968", 0.2: "#fe9f6d"}
FATE_COL = {"tumour": "#3b0f70", "extinct": "#b9b8b3", "none": "#e6e5e0"}
INK2, GRID = "#52514e", "#e6e5e0"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2})


def letter(ax, s):
    ax.text(-0.17, 1.04, s, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")


def curve_panel(ax, inp, data):
    for c, lab, _ in RUNS[inp]:
        d = data[(inp, c)]
        last = min(d.age.max(), T_END); grid = np.linspace(0, last, 500)
        tum = np.sort(d.loc[(d.fate == "tumour") & (d.age <= T_END), "age"].values)
        k = np.searchsorted(tum, grid, side="right"); n = len(d); lo, hi = wilson(k, n)
        ax.fill_between(grid, np.nan_to_num(lo), np.nan_to_num(hi), color=COL[c], alpha=0.22, lw=0)
        ax.plot(grid, k / n, color=COL[c], lw=1.6, label=lab)
        ax.text(last + 1, k[-1] / n, f"{100 * k[-1] / n:.1f}%", color=COL[c], fontsize=7, va="center")
        if len(tum):
            m = np.median(tum)
            ax.plot(m, np.searchsorted(tum, m, side="right") / n, "o", color=COL[c], ms=5, mec="white", mew=0.8, zorder=5)
    ax.set_xlim(0, 105); ax.set_ylim(0, 1); ax.set_xticks(range(0, 100, 20))
    ax.set_xlabel("Age (years)"); ax.set_ylabel("Cumulative occurrence\n(proportion of runs)")
    ax.set_title(f"{3 if inp < 0.05 else 10} events per division (input {inp:g})", fontsize=8, color=INK2, loc="left")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)


data = {}
for inp, runs in RUNS.items():
    for c, _, f in runs:
        d = load(f)
        d.loc[(d.fate == "tumour") & (d.age > T_END), "fate"] = "none"
        data[(inp, c)] = d

fig = plt.figure(figsize=(9.2, 7.4))
gs = GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.38)
axA, axB, axC, axD = (fig.add_subplot(gs[i // 2, i % 2]) for i in range(4))
curve_panel(axA, 0.03, data); curve_panel(axB, 0.10, data)

# C: fate of every run
labels, x = [], 0
for inp in RUNS:
    for c, _, _ in RUNS[inp]:
        d = data[(inp, c)]; bottom = 0
        for fate in ("tumour", "extinct", "none"):
            v = (d.fate == fate).mean()
            axC.bar(x, v, bottom=bottom, color=FATE_COL[fate], edgecolor="white", lw=0.5,
                    label=fate if x == 0 else None)
            bottom += v
        labels.append(f"{c:g}"); x += 1
    x += 0.6
xs = [0, 1, 2, 3.6, 4.6, 5.6]
axC.set_xticks(xs); axC.set_xticklabels(labels)
axC.text(1, -0.24, "input 0.03", ha="center", transform=axC.get_xaxis_transform(), color=INK2)
axC.text(4.6, -0.24, "input 0.10", ha="center", transform=axC.get_xaxis_transform(), color=INK2)
axC.set_xlabel("CTOR", labelpad=2); axC.set_ylabel("Proportion of runs"); axC.set_ylim(0, 1)
axC.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=FATE_COL[f]) for f in ("tumour", "extinct", "none")],
           labels=["tumour", "extinction", "neither (censored)"], frameon=False, fontsize=7,
           loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3)

# D: residual vs PH expectation, CTOR 0.02 vs 0.01
cmp_ = pd.read_csv(OUT / "turnover" / "exp2_matched_comparisons.csv")
cmp_ = cmp_[cmp_.comparison.str.startswith("CTOR 0.02")].reset_index(drop=True)
for i, r in cmp_.iterrows():
    axD.errorbar(i, r.residual, yerr=[[r.residual - r.res_lo], [r.res_hi - r.residual]], fmt="o", color=COL[0.02],
                 ms=6, elinewidth=1.4, capsize=3)
axD.axhline(0, color=INK2, lw=0.9)
axD.set_xticks(range(len(cmp_))); axD.set_xticklabels(["input 0.03", "input 0.10"]); axD.set_xlim(-0.6, len(cmp_) - 0.4)
axD.set_ylabel("Observed − PH-predicted\nMEDA (years)")
axD.set_title("CTOR 0.02 vs 0.01 at matched input", fontsize=8, color=INK2, loc="left")
axD.grid(True, axis="y", color=GRID, lw=0.6); axD.set_axisbelow(True)

for ax, s in zip((axA, axB, axC, axD), "ABCD"):
    letter(ax, s)
fig.savefig(FIG / "FigureS4_turnover.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "FigureS4_turnover.pdf", bbox_inches="tight")
print("written", FIG / "FigureS4_turnover.png")
