"""
Revised Figure 3, same layout as the submitted figure (A-I) plus a fourth row (J, K).
A-E, G, H: cumulative occurrence of esiTumor (proportion of runs) with 95% Wilson bands; dot = MEDA.
F, I: observed minus PH-predicted MEDA for adjacent input steps (1, 2, 3 as in panel A), bootstrap 95% CI
      (output/pairs_ph_benchmark.csv, T = 92 y).
J: inherited first hit, carriers vs non-carriers (germline/, 3,000 runs per arm).
K: hazard ratio by 10-year age band, carrier vs non-carrier and the matched-input comparator (3 -> 4 events).
Output: figures/Figure3_revised.png and .pdf
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

from ph_benchmark import SERIES, load, wilson, OUT  # noqa: E402
from germline_analysis import load_groups, hr_trend  # noqa: E402

from paths import FIG, GERMLINE  # noqa: E402
T_END = 92.0
LEVEL_COL = ["#3b0f70", "#8c2981", "#de4968", "#fe9f6d"]     # magma, as in the submitted figure
STEP_COL = ["#440154", "#2a8f84", "#e3d33b"]                 # steps 1, 2, 3 (viridis), as in old F/I
NONCARRIER, CARRIER, GREY = "#2a78d6", "#eb6834", "#8a8984"
INK2, GRID = "#52514e", "#e6e5e0"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2})

PANELS = {  # panel -> series, legend unit, corner note
    "A": ("A_constant", "constant", None),
    "B": ("B_increasing", "increasing", None),
    "C": ("C_multistage", "multistage", None),
    "D": ("D_constant_aging", "constant + aging", None),
    "E": ("E_incr_aging_ctor001", "increasing + aging", "CTOR = 0.01"),
    "G": ("G_incr_aging_ctor002", "increasing + aging", "CTOR = 0.02"),
    "H": ("H_incr_aging_ctor02", "increasing + aging", "CTOR = 0.2"),
}


def letter(ax, s):
    ax.text(-0.17, 1.04, s, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")


def cum_curve(d, grid):
    tum = np.sort(d.loc[(d.fate == "tumour") & (d.age <= T_END), "age"].values)
    k = np.searchsorted(tum, grid, side="right"); n = len(d)
    lo, hi = wilson(k, n)
    return k / n, np.nan_to_num(lo), np.nan_to_num(hi, nan=0), tum


def occurrence_panel(ax, series, title, note):
    folder, files = SERIES[series]
    meda_pts = []
    for i, (prefix, level) in enumerate(files.items()):
        d = load(folder, prefix)
        last = min(d.age.max(), T_END)            # draw until the last run ended, as in the original
        grid = np.linspace(0, last, 500)
        p, lo, hi, tum = cum_curve(d, grid)
        col = LEVEL_COL[i]
        ax.fill_between(grid, lo, hi, color=col, alpha=0.22, lw=0)
        ax.plot(grid, p, color=col, lw=1.6, label=str(level))
        ax.text(last + 1, p[-1], f"{100 * p[-1]:.1f}%", color=col, fontsize=7, va="center")
        if len(tum):
            m = np.median(tum); pm = np.searchsorted(tum, m, side="right") / len(d)
            ax.plot(m, pm, "o", color=col, ms=5, mec="white", mew=0.8, zorder=5)
            meda_pts.append((m, pm))
    ax.set_xlim(0, 105); ax.set_ylim(0, 1.0); ax.set_xticks(range(0, 100, 20))
    ax.set_xlabel("Age (years)"); ax.set_ylabel("Cumulative occurrence\n(proportion of runs)")
    ax.legend(title="Events per division", ncol=2, frameon=False, fontsize=7, title_fontsize=7,
              loc="lower right" if series.startswith("H") else "upper left",
              bbox_to_anchor=(1, 0.1) if series.startswith("H") else None, handlelength=1.2, columnspacing=0.8)
    ax.set_title(title, fontsize=8, color=INK2, loc="right")
    if note:
        ax.text(0.97, 0.05, note, transform=ax.transAxes, fontsize=7.5, color=INK2, ha="right")
    ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    return meda_pts


def residual_panel(ax, series_list, names):
    pr = pd.read_csv(OUT / "pairs_ph_benchmark.csv")
    pr = pr[pr["T"] == 92]
    w = 0.22
    for g, s in enumerate(series_list):
        levels = list(SERIES[s][1].values())
        for j, (a, b) in enumerate(zip(levels[:-1], levels[1:])):
            r = pr[(pr.series == s) & (pr.input0 == a) & (pr.input1 == b)].iloc[0]
            x = g + (j - 1) * w
            ax.errorbar(x, r.residual_y, yerr=[[r.residual_y - r.residual_y_lo], [r.residual_y_hi - r.residual_y]],
                        fmt="o", color=STEP_COL[j], ms=5.5, elinewidth=1.4, capsize=2,
                        label=f"step {j + 1}" if g == 0 else None)
    ax.axhline(0, color=INK2, lw=0.9)
    ax.set_xticks(range(len(series_list))); ax.set_xticklabels(names)
    ax.set_xlim(-0.6, len(series_list) - 0.4)
    ax.set_ylabel("Observed − PH-predicted\nMEDA (years)")
    ax.legend(frameon=False, fontsize=7, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3,
              handletextpad=0.2, columnspacing=0.8)
    ax.grid(True, axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)


fig = plt.figure(figsize=(10.5, 12.6))
gs = GridSpec(4, 3, figure=fig, hspace=0.48, wspace=0.42)
axes = {k: fig.add_subplot(gs[i // 3, i % 3]) for i, k in enumerate("ABCDEFGHI")}
axes["J"] = fig.add_subplot(gs[3, 0]); axes["K"] = fig.add_subplot(gs[3, 1])

for k, (series, title, note) in PANELS.items():
    pts = occurrence_panel(axes[k], series, title, note)
    if k == "A":   # the three adjacent steps compared in F and I
        for j, ((x0, y0), (x1, y1)) in enumerate(zip(pts[:-1], pts[1:])):
            axes[k].annotate("", xy=(x1, y1), xytext=(x0, y0),
                             arrowprops=dict(arrowstyle="->", color=STEP_COL[j], lw=1.8))
            axes[k].text((x0 + x1) / 2 + 1.5, (y0 + y1) / 2 + 0.02, str(j + 1), color=STEP_COL[j],
                         fontsize=8, fontweight="bold")

residual_panel(axes["F"], ["A_constant", "B_increasing", "C_multistage", "D_constant_aging"], ["A", "B", "C", "D"])
residual_panel(axes["I"], ["E_incr_aging_ctor001", "G_incr_aging_ctor002", "H_incr_aging_ctor02"],
               ["E\n0.01", "G\n0.02", "H\n0.2"])

# J, K: inherited first hit
gr = load_groups(GERMLINE)
d0, d1, d4 = gr[("003", "noncarrier")], gr[("003", "carrier")], gr[("004", "noncarrier")]
ax = axes["J"]
for d, col, lab in ((d0, NONCARRIER, "non-carrier (3 events)"), (d1, CARRIER, "carrier (3 events)"),
                    (d4, GREY, "non-carrier (4 events)")):
    grid = np.linspace(0, T_END, 500)
    p, lo, hi, tum = cum_curve(d, grid)
    ls = "--" if col == GREY else "-"
    ax.fill_between(grid, lo, hi, color=col, alpha=0.22, lw=0)
    ax.plot(grid, p, color=col, lw=1.6, ls=ls, label=lab)
    ax.text(T_END + 1, p[-1], f"{100 * p[-1]:.1f}%", color=col, fontsize=7, va="center")
    m = np.median(tum); pm = np.searchsorted(tum, m, side="right") / len(d)
    ax.plot(m, pm, "o", color=col, ms=5, mec="white", mew=0.8, zorder=5)
ax.set_xlim(0, 105); ax.set_ylim(0, 1.0); ax.set_xticks(range(0, 100, 20))
ax.set_xlabel("Age (years)"); ax.set_ylabel("Cumulative occurrence\n(proportion of runs)")
ax.legend(title="Inherited first hit (multistage)", frameon=False, fontsize=7, title_fontsize=7, loc="upper left")
ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)

ax = axes["K"]
for dd, col, lab, off in ((d1, CARRIER, "carrier / non-carrier", -1.2), (d4, GREY, "4 / 3 events per division", 1.2)):
    tr, curve = hr_trend(d0, dd)
    x, y, lo, hi = map(np.array, zip(*curve))
    ax.errorbar(x + off, y, yerr=[y - lo, hi - y], fmt="o-", color=col, ms=4.5, lw=1.4, elinewidth=1, capsize=2,
                label=f"{lab}: ×{tr['trend']:.2f} per decade ({tr['trend_lo']:.2f}–{tr['trend_hi']:.2f})")
ax.set_yscale("log"); ax.axhline(1, color=INK2, lw=0.9); ax.set_ylim(0.8, 400)
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}"))
ax.set_xlabel("Age (years)"); ax.set_ylabel("Hazard ratio")
ax.legend(frameon=False, fontsize=6.8, loc="upper left", bbox_to_anchor=(0, 1.0))
ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)

for k, ax in axes.items():
    letter(ax, k)

fig.savefig(FIG / "Figure3_revised.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "Figure3_revised.pdf", bbox_inches="tight")
print("written", FIG / "Figure3_revised.png")
