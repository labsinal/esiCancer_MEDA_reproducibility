"""
Figure S2: validation of the PH benchmark and of the hazard-ratio trend on artificial data with a known answer.

Detection ages are drawn from the classical multistage model (cumulative hazard c t^k), as in
germline_analysis.py --synthetic, at the occurrences of the inherited first-hit experiment
(reference about 20%, comparison about 60% by 92 y) and 3,000 runs per group.
  A, B  proportional-hazards control: reference hazard (k = 5) multiplied by a constant.
  C, D  inherited first hit: one of the five steps present from birth (k = 4).
  E     share of 150 repetitions in which the change per decade was detected (95% CI excludes 1),
        by runs per group (output/germline/germline_power_check.csv, from power_check.py).
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

import germline_analysis as ga
from germline_analysis import synthetic, hr_trend, T_END
from paths import FIG
from ph_benchmark import at_endpoint, pair_stats, boot_pair

P0, P1, N = 0.20, 0.60, 3000
REF, ALT, GREY = "#2a78d6", "#eb6834", "#8a8984"
INK2, GRID = "#52514e", "#e6e5e0"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2})


def letter(ax, s):
    ax.text(-0.2, 1.04, s, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")


def occurrence(ax, d0, d1, lab1, title):
    grid = np.linspace(0, T_END, 400)
    a0, t0, _ = at_endpoint(d0, T_END); a1, t1, _ = at_endpoint(d1, T_END)
    s = pair_stats(a0, t0, a1, t1)
    ci = boot_pair(a0, t0, a1, t1, ["residual_y"])["residual_y"]
    for d, col, lab, med in ((d0, REF, "reference (k = 5)", s["meda0"]), (d1, ALT, lab1, s["meda1"])):
        tum = np.sort(d.loc[d.fate == "tumour", "age"].values)
        p = np.searchsorted(tum, grid, side="right") / len(d)
        ax.plot(grid, p, color=col, lw=2, label=lab)
        ax.plot(med, np.searchsorted(tum, med, side="right") / len(d), "o", color=col, mec="white", ms=6, zorder=3)
    ax.set_xlim(0, T_END); ax.set_ylim(0, 0.7)
    ax.set_xlabel("Age (years)"); ax.set_ylabel("Cumulative occurrence\n(proportion of runs)")
    ax.legend(frameon=False, loc="upper left", handlelength=1.2)
    ax.text(0.02, 0.62, f"observed − PH-predicted MEDA\n{s['residual_y']:+.1f} y (95% CI {ci[0]:+.1f} to {ci[1]:+.1f})",
            transform=ax.transAxes, fontsize=7, color=INK2, va="top")
    ax.set_title(title, loc="right", fontsize=8, color=INK2)


def hazard_ratio(ax, d0, d1):
    tr, curve = hr_trend(d0, d1)
    x, y, lo, hi = map(np.array, zip(*curve))
    ax.errorbar(x, y, yerr=[y - lo, hi - y], fmt="o-", color=ALT, ms=4, lw=1.5, elinewidth=1, capsize=2)
    ax.axhline(1, color=INK2, lw=0.8)
    ax.set_yscale("log"); ax.set_xlabel("Age (years)"); ax.set_ylabel("Hazard ratio")
    ax.text(0.98, 0.95, f"×{tr['trend']:.2f} per decade\n(95% CI {tr['trend_lo']:.2f}–{tr['trend_hi']:.2f})",
            transform=ax.transAxes, ha="right", va="top", fontsize=7, color=INK2)


def power(ax):
    pw = pd.read_csv(ga.OUT / "germline_power_check.csv")
    groups = [("PH", 0.05, "PH\n5→30%"), ("PH", 0.2, "PH\n20→60%"),
              ("multistage", 0.05, "hit\n5→30%"), ("multistage", 0.2, "hit\n20→60%")]
    w = 0.38
    for i, (scen, p0, _) in enumerate(groups):
        for j, (n, col) in enumerate(((1000, GREY), (3000, ALT))):
            v = pw[(pw.scenario == scen) & (pw.p0 == p0) & (pw.n_per_arm == n)].share_trend_CI_excludes_1.iloc[0]
            ax.bar(i + (j - 0.5) * w, v, w * 0.92, color=col, label=f"{n:,} runs per group" if i == 0 else None)
            ax.text(i + (j - 0.5) * w, v + 0.015, f"{v:.0%}", ha="center", fontsize=6, color=INK2)
    ax.axhline(0.05, color=INK2, lw=0.8, ls="--")
    ax.set_xticks(range(len(groups))); ax.set_xticklabels([g[2] for g in groups], fontsize=7)
    ax.set_ylim(0, 1.3); ax.set_yticks(np.arange(0, 1.01, 0.2)); ax.set_ylabel("Change with age detected\n(share of 150 repetitions)")
    ax.set_xlabel("PH control or first hit; occurrence")
    ax.legend(frameon=False, loc="upper left")


def main():
    scale = np.log(1 - P1) / np.log(1 - P0)
    ref = synthetic(5, P0, N)
    ph = synthetic(5, P0, N, scale=scale)
    hit = synthetic(4, P1, N)

    fig = plt.figure(figsize=(7.2, 5.6))
    gs = GridSpec(2, 3, figure=fig, wspace=0.55, hspace=0.5, width_ratios=[1, 1, 1.15])
    ax = {k: fig.add_subplot(gs[r, c]) for k, (r, c) in
          dict(A=(0, 0), B=(0, 1), C=(1, 0), D=(1, 1)).items()}
    ax["E"] = fig.add_subplot(gs[:, 2])
    occurrence(ax["A"], ref, ph, "risk × constant", "proportional hazards")
    hazard_ratio(ax["B"], ref, ph)
    occurrence(ax["C"], ref, hit, "first hit (k = 4)", "inherited first hit")
    hazard_ratio(ax["D"], ref, hit)
    power(ax["E"])
    for k, a in ax.items():
        letter(a, k); a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
    for ext, kw in (("png", dict(dpi=300)), ("pdf", {})):
        fig.savefig(FIG / f"FigureS2_synthetic.{ext}", bbox_inches="tight", **kw)


if __name__ == "__main__":
    main()
