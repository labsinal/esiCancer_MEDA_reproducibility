"""
Figure S3: sensitivity of the PH benchmark to the end of follow-up.
Observed minus PH-predicted MEDA for the adjacent input steps of Figure 3F, I, with follow-up ending at
T = 60, 75 and 92 years (output/pairs_ph_benchmark.csv). Comparisons that could not be estimated at a given T
(too few tumours) are marked n.e.
"""
import matplotlib.pyplot as plt
import pandas as pd

from paths import FIG, OUT

SERIES = [("A_constant", "A"), ("B_increasing", "B"), ("C_multistage", "C"), ("D_constant_aging", "D"),
          ("E_incr_aging_ctor001", "E"), ("G_incr_aging_ctor002", "G"), ("H_incr_aging_ctor02", "H")]
STEP_COL = ["#440154", "#2a8f84", "#e3d33b"]
INK2, GRID = "#52514e", "#e6e5e0"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2})


def main():
    pr = pd.read_csv(OUT / "pairs_ph_benchmark.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.8), sharey=True, gridspec_kw=dict(wspace=0.08))
    for ax, T in zip(axes, (60, 75, 92)):
        t = pr[pr["T"] == T]
        for g, (s, name) in enumerate(SERIES):
            levels = sorted(set(pr[(pr.series == s)].input0) | set(pr[(pr.series == s)].input1))
            for j, (a, b) in enumerate(zip(levels[:-1], levels[1:])):
                x = g + (j - 1) * 0.24
                r = t[(t.series == s) & (t.input0 == a) & (t.input1 == b)]
                if r.empty or r.residual_y.isna().iloc[0]:
                    ax.plot(x, 0, "x", color=STEP_COL[j], ms=4, mew=1.2)
                    continue
                r = r.iloc[0]
                ax.errorbar(x, r.residual_y, yerr=[[r.residual_y - r.residual_y_lo], [r.residual_y_hi - r.residual_y]],
                            fmt="o", color=STEP_COL[j], ms=4, elinewidth=1.2, capsize=1.5,
                            label=f"step {j + 1}" if (g == 0 and T == 60) else None)
        ax.axhline(0, color=INK2, lw=0.8)
        ax.set_xticks(range(len(SERIES))); ax.set_xticklabels([n for _, n in SERIES])
        ax.set_xlim(-0.6, len(SERIES) - 0.4)
        ax.set_title(f"follow-up to {T} years", fontsize=8, color=INK2)
        ax.set_xlabel("Series (Figure 3)")
        ax.grid(True, axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
    axes[0].set_ylabel("Observed − PH-predicted\nMEDA (years)")
    axes[0].plot([], [], "x", color=INK2, label="not estimable (too few tumours)")
    fig.legend(*axes[0].get_legend_handles_labels(), frameon=False, loc="upper center", ncol=4,
               bbox_to_anchor=(0.5, 1.08), handletextpad=0.3, columnspacing=1.2)
    for ext, kw in (("png", dict(dpi=300)), ("pdf", {})):
        fig.savefig(FIG / f"FigureS3_endpoint.{ext}", bbox_inches="tight", **kw)


if __name__ == "__main__":
    main()
