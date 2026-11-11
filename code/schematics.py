"""
Schematics replacing the screenshots of Figure 2A and Figure S1 (reviewer: legible schematics).
2A  How esiCancer works: inputs (esiTable, configuration), one generation of one cell, and the three ends of a run.
S1  esiCancer guide: the esiTable columns with example records of each type (repo/esitables/esiTable_escc.csv),
    and the simulation configuration used in this study (repo/esicancer/CONFIGURATIONS.md).
Colours match Figure 2B (esiTumor blue, extinction grey, neither orange).
"""
import matplotlib.pyplot as plt
import textwrap

from matplotlib.patches import Ellipse, FancyArrowPatch, FancyBboxPatch

from paths import FIG

INK, INK2, LINE, PANEL = "#1f1e1c", "#52514e", "#bdbbb4", "#f6f5f1"
TUMOR, EXTINCT, NEITHER = "#2a78d6", "#8a8984", "#eb6834"
DRIVER, ESSENTIAL, PASSENGER, AGING = "#c8442b", "#7b4fa8", "#8a8984", "#2f8f6b"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans"})


def box(ax, x, y, w, h, fc=PANEL, ec=LINE, lw=0.8, r=0.012):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw))


def arrow(ax, a, b, color=INK2, rad=0.0, lw=1.0):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=9, color=color, lw=lw,
                                 connectionstyle=f"arc3,rad={rad}", shrinkA=0, shrinkB=0))


def dot(ax, x, y, r, **kw):
    """A circle in axes-fraction coordinates on a non-square figure (r in units of figure width)."""
    w, h = ax.figure.get_size_inches()
    ax.add_patch(Ellipse((x, y), 2 * r, 2 * r * w / h, **kw))


def cell(ax, x, y, r=0.022, dots=()):
    dot(ax, x, y, r, fc="#ffffff", ec=INK2, lw=0.9)
    w, h = ax.figure.get_size_inches()
    for dx, dy, c in dots:
        dot(ax, x + dx * r, y + dy * r * w / h, r * 0.28, fc=c, ec="none")


def title(ax, x, y, s, size=9):
    ax.text(x, y, s, fontsize=size, fontweight="bold", color=INK, va="top")


def figure2a():
    fig = plt.figure(figsize=(7.2, 3.6))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.005, 0.985, "A", fontsize=14, fontweight="bold", va="top")
    S = 6.3

    # ---- inputs
    box(ax, 0.03, 0.04, 0.265, 0.88)
    title(ax, 0.045, 0.9, "Inputs")
    ax.text(0.045, 0.835, "esiTable: one record\nper gene event", fontsize=7.2, color=INK, va="top", linespacing=1.1)
    types = [(DRIVER, "Driver (58)", "raises division, lowers death"),
             (ESSENTIAL, "Essential (47)", "lowers division, raises death"),
             (PASSENGER, "Passenger (74)", "no effect on fitness"),
             (AGING, "Aging (4)", "all cells, at 60, 70, 80, 90 y")]
    for i, (c, name, eff) in enumerate(types):
        y = 0.735 - i * 0.085
        dot(ax, 0.055, y - 0.013, 0.007, fc=c, ec="none")
        ax.text(0.07, y, name, fontsize=7.2, fontweight="bold", color=c, va="top")
        ax.text(0.07, y - 0.037, eff, fontsize=S, color=INK2, va="top")
    ax.text(0.045, 0.38, "Each record has an event weight\n(PROBEVENT), a dominance and\neffects on division and death",
            fontsize=S, color=INK2, va="top", linespacing=1.15)
    ax.plot([0.045, 0.28], [0.265, 0.265], color=LINE, lw=0.6)
    ax.text(0.045, 0.245, "Configuration", fontsize=7.2, color=INK, va="top")
    ax.text(0.045, 0.2, "1,000 initial cells\nP(division) = P(death) = 0.01\nevents per division, constant\n"
            "or increasing with age\n1 generation = 1 week", fontsize=S, color=INK2, va="top", linespacing=1.15)

    arrow(ax, (0.3, 0.5), (0.325, 0.5))

    # ---- one generation of one cell
    box(ax, 0.33, 0.04, 0.40, 0.88)
    title(ax, 0.345, 0.9, "Each generation, for each cell")
    cell(ax, 0.405, 0.63, r=0.032, dots=[(-0.35, 0.25, DRIVER)])
    ax.text(0.405, 0.555, "cell with its own\nP(division), P(death)", fontsize=S, color=INK2, ha="center", va="top")
    arrow(ax, (0.442, 0.67), (0.515, 0.77))
    arrow(ax, (0.442, 0.63), (0.515, 0.63))
    arrow(ax, (0.442, 0.59), (0.515, 0.49))
    ax.text(0.468, 0.74, "divides", fontsize=S, color=INK2, ha="right")
    ax.text(0.479, 0.642, "neither", fontsize=S, color=INK2, ha="center")
    ax.text(0.488, 0.565, "dies", fontsize=S, color=INK2, ha="left")
    cell(ax, 0.537, 0.785, r=0.018, dots=[(-0.35, 0.25, DRIVER)])
    cell(ax, 0.582, 0.785, r=0.018, dots=[(-0.35, 0.25, DRIVER), (0.35, -0.25, ESSENTIAL)])
    cell(ax, 0.537, 0.63, r=0.018, dots=[(-0.35, 0.25, DRIVER)])
    ax.text(0.537, 0.49, "×", fontsize=12, color=EXTINCT, ha="center", va="center")
    ax.text(0.608, 0.785, "k events per\ndivision", fontsize=S, color=INK2, va="center", linespacing=1.1)
    box(ax, 0.35, 0.1, 0.36, 0.28, fc="#ffffff")
    ax.text(0.362, 0.36, "At each division", fontsize=7.2, fontweight="bold", color=INK, va="top")
    ax.text(0.362, 0.31,
            "1. each of the k events hits an esiTable\n    record, in proportion to its PROBEVENT\n"
            "2. P(division) and P(death) of the cell are\n    updated by the record's effects, scaled\n"
            "    by its dominance when one allele is hit\n3. both daughter cells carry the new values",
            fontsize=S, color=INK2, va="top", linespacing=1.2)
    arrow(ax, (0.582, 0.745), (0.582, 0.39))

    arrow(ax, (0.735, 0.5), (0.76, 0.5))

    # ---- end of a run
    box(ax, 0.765, 0.04, 0.225, 0.88)
    title(ax, 0.78, 0.9, "End of a run")
    ax.text(0.78, 0.84, "one run = one individual", fontsize=S, color=INK2, va="top")
    ends = [(TUMOR, "esiTumor", "reaches 10,000 cells\n(10× the start);\nage = detection age"),
            (EXTINCT, "Extinction", "no cells left"),
            (NEITHER, "Neither", "92 years reached\n(4,800 generations);\ncensored at 92 y")]
    for i, (c, name, txt) in enumerate(ends):
        y = 0.75 - i * 0.235
        if i == 0:
            for dx, dy in [(0, 0), (0.014, 0.022), (0.014, -0.022), (0.028, 0), (0.028, 0.044), (0.0, 0.044)]:
                dot(ax, 0.795 + dx, y - 0.05 + dy, 0.0075, fc=c, ec="white", lw=0.4)
        elif i == 1:
            dot(ax, 0.81, y - 0.03, 0.016, fc="none", ec=c, lw=1.0, ls=(0, (2, 1.5)))
        else:
            for dx in (0, 0.02):
                dot(ax, 0.8 + dx, y - 0.03, 0.0075, fc=c, ec="white", lw=0.4)
        ax.text(0.845, y, name, fontsize=7.2, fontweight="bold", color=c, va="top")
        ax.text(0.845, y - 0.042, txt, fontsize=S, color=INK2, va="top", linespacing=1.15)
    return fig


ROWS = [  # GENE, EVENT, PROBEVENT, DOMINANCE, DIVRATE, DEATHRATE, type
    ("TP53", "mutation", "1", "0.342", "ADD +0.0752", "ADD −0.0752", DRIVER, "driver"),
    ("CCND1", "amplification", "1", "0.863", "ADD +0.1292", "ADD −0.1292", DRIVER, "driver"),
    ("CDKN2A", "deletion", "1", "0.829", "ADD +0.0854", "ADD −0.0854", DRIVER, "driver"),
    ("MDN1", "mutation", "1", "0.68", "ADD −0.005", "ADD +0.005", ESSENTIAL, "essential"),
    ("TTN", "mutation", "45", "1", "ADD 0", "ADD 0", PASSENGER, "passenger"),
    ("aging_1", "mutation", "0", "1", "MULT ×0.9", "MULT ×1.1", AGING, "aging"),
]
COLS = [("GENE", 0.0), ("EVENT", 0.12), ("PROBEVENT", 0.27), ("DOMINANCE", 0.39), ("DIVRATE", 0.52),
        ("DEATHRATE", 0.665), ("TYPE", 0.83)]


def figure_s1():
    fig = plt.figure(figsize=(7.2, 6.0))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    S, LH = 6.5, 6.5 * 1.25 / 72 / 6.0   # font size, line height in axes fraction

    # ---- esiTable
    box(ax, 0.01, 0.6, 0.98, 0.385)
    title(ax, 0.025, 0.972, "esiTable (Table S2): one record per gene event")
    x0, ytop = 0.035, 0.905
    for name, dx in COLS:
        ax.text(x0 + dx, ytop, name, fontsize=7, fontweight="bold", color=INK, va="center")
    ax.plot([x0, 0.965], [ytop - 0.02, ytop - 0.02], color=INK2, lw=0.7)
    for i, row in enumerate(ROWS):
        y = ytop - 0.05 - i * 0.042
        c = row[6]
        ax.add_patch(FancyBboxPatch((x0 - 0.01, y - 0.018), 0.94, 0.036, boxstyle="round,pad=0,rounding_size=0.005",
                                    fc=c, ec="none", alpha=0.10))
        for (_, dx), v in zip(COLS, row[:6]):
            ax.text(x0 + dx, y, v, fontsize=7, color=INK, va="center")
        ax.text(x0 + COLS[6][1], y, row[7], fontsize=7, color=c, va="center", fontweight="bold")
    ax.text(x0, ytop - 0.05 - 6 * 0.042 + 0.01,
            "…  183 records in total: 58 driver, 47 essential, 74 passenger and 4 aging records, in 173 genes",
            fontsize=S, color=INK2, va="center")

    # ---- column definitions
    box(ax, 0.01, 0.015, 0.53, 0.57)
    title(ax, 0.025, 0.57, "What the columns mean")
    notes = [
        ("PROBEVENT", "relative event weight: the chance that an event hits a record is proportional to it; 1 for "
                      "drivers and essential genes, 5–45 for passengers (mutational target size)"),
        ("DOMINANCE", "fraction of the effect applied when one allele is hit, from the 20/20 rule"),
        ("DIVRATE, DEATHRATE", "effect on the cell's probability of division and of death, added (ADD) or "
                               "multiplied (MULT); set from the frequency of the event in ESCC tumors"),
        ("Genome size", "16,647,721, the sum of all gene lengths; oncogenome size 58, the sum of driver weights"),
        ("Not used", "mutation-rate and microenvironment columns (inactive or neutral)"),
    ]
    y = 0.515
    for k, v in notes:
        ax.text(0.03, y, k, fontsize=7, fontweight="bold", color=INK, va="top")
        lines = textwrap.wrap(v, 62)
        ax.text(0.03, y - LH * 1.05, "\n".join(lines), fontsize=S, color=INK2, va="top", linespacing=1.15)
        y -= LH * (1.25 + 1.15 * len(lines)) + 0.012

    # ---- configuration
    box(ax, 0.55, 0.015, 0.44, 0.57)
    title(ax, 0.565, 0.57, "Simulation configuration (Table S3)")
    cfg = [("Initial cells", "1,000"),
           ("P(division) = P(death)", "0.01 (CTOR 0.01); 0.02 or 0.2 in Figure 3G, H"),
           ("Maximum divisions", "1,000 per cell"),
           ("Events per division", "constant 1, 3, 5 or 10; or increasing, 1 + ⌊age / step⌋ with steps of 20, 10, 7, 5 y"),
           ("Manual events", "aging_1–4 in all cells at 60, 70, 80, 90 y"),
           ("Six-gene model", "6 genes, genome 750; if both genes of a pair are hit, division effect ×3, death effect ÷3"),
           ("Stop", "10,000 cells (esiTumor), 0 cells, or 4,800 generations"),
           ("Time", "1 generation = 1 week; 4,800 generations = 92 y"),
           ("Runs", "seeds 1–1,000 per condition; 1–3,000 for the inherited first hit")]
    y = 0.515
    for n, (k, v) in enumerate(cfg):
        lines = textwrap.wrap(v, 34)
        ax.text(0.565, y, k, fontsize=7, fontweight="bold", color=INK, va="top")
        ax.text(0.74, y, "\n".join(lines), fontsize=S, color=INK2, va="top", linespacing=1.15)
        y -= LH * 1.15 * max(1, len(lines)) + 0.012
        if n < len(cfg) - 1:
            ax.plot([0.565, 0.975], [y + 0.007, y + 0.007], color=LINE, lw=0.4)
    return fig


if __name__ == "__main__":
    for fn, name in ((figure2a, "Figure2A_schematic"), (figure_s1, "FigureS1_guide")):
        f = fn()
        for ext in ("png", "pdf"):
            f.savefig(FIG / f"{name}.{ext}", dpi=300)
        plt.close(f)
        print("saved", FIG / name)
