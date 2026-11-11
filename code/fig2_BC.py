"""
Figure 2B, C (panel A is a hand-drawn schematic).
B  Single runs with aging events (seeds 1-14; ESCC esiTable, 20 events per division, CTOR 0.01): cell
   population against age, coloured by fate (esiTumor, extinction, neither by 92 y); threshold at 10,000 cells;
   aging events at 60, 70, 80 and 90 years. Data: data/simulations/revision/single_runs_aging/
   (Seed, Generation, Population extracted from the esiCancer ancestral output of each run).
C  Events per division: constant (1, 3, 5, 10) and increasing, 1 + floor(age / step) with steps of 20, 10, 7 and
   5 years (levels 5, 10, 14, 20 of Figure 3B, E, G, H; level = 100 / step, a nominal label).
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from paths import FIG, ROOT

SINGLES = ROOT / "data" / "simulations" / "revision" / "single_runs_aging" / "single_runs_aging_population.csv"
WEEKS_PER_YEAR, THRESHOLD, T_END = 52.143, 10_000, 92.0
FATE_COL = {"esiTumor": "#2a78d6", "extinction": "#8a8984", "neither by 92 y": "#eb6834"}
LEVEL_COL = ["#3b0f70", "#8c2981", "#de4968", "#fe9f6d"]
INK2, GRID = "#52514e", "#e6e5e0"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2})


def letter(ax, s):
    ax.text(-0.14, 1.04, s, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")


def runs_panel(ax):
    for _, d in pd.read_csv(SINGLES).groupby("Seed"):
        age, pop = d.Generation / WEEKS_PER_YEAR, d.Population.clip(lower=1)
        last = d.Population.iloc[-1]
        fate = "esiTumor" if last >= THRESHOLD else "extinction" if last == 0 else "neither by 92 y"
        ax.plot(age, pop, color=FATE_COL[fate], lw=1)
    for a in (60, 70, 80, 90):
        ax.axvline(a, color=INK2, lw=0.8, ls=":")
    ax.axhline(THRESHOLD, color=INK2, lw=0.8, ls="--")
    ax.text(1, THRESHOLD * 1.15, "esiTumor threshold (10× initial population)", fontsize=7, color=INK2)
    ax.text(60.5, 1.6, "aging events", fontsize=7, color=INK2)
    ax.set_yscale("log"); ax.set_ylim(1, 3e4); ax.set_xlim(0, T_END)
    ax.set_xlabel("Age (years)"); ax.set_ylabel("Number of cells")
    for fate, col in FATE_COL.items():
        ax.plot([], [], color=col, lw=2, label=fate)
    ax.legend(frameon=False, loc="lower left", title="Fate of the run", title_fontsize=7)


def events_panel(ax):
    age = np.linspace(0, T_END, 1000)
    for c, col in zip((1, 3, 5, 10), LEVEL_COL):
        ax.plot(age, np.full_like(age, c), color=col, lw=1.5, ls="--")
    for (level, step), col in zip(((5, 20), (10, 10), (14, 7), (20, 5)), LEVEL_COL):
        ax.step(age, 1 + age // step, where="post", color=col, lw=1.8, label=f"{level} (+1 every {step} y)")
    ax.set_xlim(0, T_END); ax.set_ylim(0, 21)
    ax.set_xlabel("Age (years)"); ax.set_ylabel("Events per division")
    ax.plot([], [], color=INK2, ls="--", label="constant: 1, 3, 5, 10")
    ax.legend(frameon=False, loc="upper left", title="Increasing, level", title_fontsize=7)


def main():
    fig, (b, c) = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw=dict(wspace=0.3))
    runs_panel(b); events_panel(c)
    for ax, k in ((b, "B"), (c, "C")):
        letter(ax, k); ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    for ext, kw in (("png", dict(dpi=300)), ("pdf", {})):
        fig.savefig(FIG / f"Figure2_BC.{ext}", bbox_inches="tight", **kw)


if __name__ == "__main__":
    main()
