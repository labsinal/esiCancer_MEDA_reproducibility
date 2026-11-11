"""
Germline first-hit simulation: carriers vs non-carriers (see RUN_PROTOCOL.md).

Usage
  python germline_analysis.py --dir ../data/simulations/germline
      reads germ_<level>_carrier_automatics_ancestralResults.csv and
            germ_<level>_noncarrier_ (or _non_carrier_) automatics_ancestralResults.csv
      a non-carrier file without a carrier partner is used only as a "more mutation input" comparator
  python germline_analysis.py --synthetic
      validation with a known answer: classic multistage model (non-carriers 5 stages, carriers 4)
      and a proportional-hazards negative control. Expected: multistage -> HR trend ~0.86/decade and
      negative MEDA residual; PH control -> trend ~1 and residual ~0.

For each level reports (T = 92 y):
  occurrence (Wilson 95% CI), MEDA (bootstrap 95% CI), observed MEDA shift, PH-expected shift,
  residual (bootstrap CI), HR at 60 y and HR change per decade (Poisson GLM, cause-specific hazard).
Also, for consecutive input levels, the comparator "more mutation input" (non-carriers at level L vs
non-carriers at the next level): the germline hit must be judged against this, within the same model,
not only against the flat line, because the 6-gene model shows some age dependence for input alone.
Writes output/germline_summary.csv and figures/germline_<label>.png.
"""
import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

HERE = Path(__file__).resolve().parent
from paths import OUT as _OUT, FIG as _FIG, GERMLINE
from ph_benchmark import WEEKS_PER_YEAR, THRESHOLD, at_endpoint, pair_stats, boot_pair, wilson  # noqa: E402

OUT = _OUT / "germline"; FIG = _FIG / "germline"
T_END = 92.0
EDGES = np.arange(0, 100, 10.0)
NONCARRIER, CARRIER = "#2a78d6", "#eb6834"
RNG = np.random.default_rng(1971)


def load_file(path):
    d = pd.read_csv(path, sep=";", usecols=[0, 1, 3])
    d.columns = ["seed", "gen", "pop"]
    d["age"] = d["gen"] / WEEKS_PER_YEAR
    d["fate"] = np.where(d["pop"] >= THRESHOLD, "tumour", np.where(d["pop"] == 0, "extinct", "none"))
    return d


def load_groups(folder):
    """Return {(level, 'carrier'|'noncarrier'): DataFrame}. Several files per group (batches of seeds, e.g.
    germ_003_carrier_automatics_... and germ_003_carrier_s1001_automatics_...) are concatenated; duplicate
    seeds are an error. Accepts 'noncarrier' and 'non_carrier'."""
    parts = {}
    for f in sorted(Path(folder).glob("germ_*automatics_ancestralResults.csv")):
        m = re.match(r"germ_([^_]+)_(non_?carrier|carrier)(?:_[A-Za-z0-9]+)?_automatics", f.name)
        if m:
            key = (m.group(1), "carrier" if m.group(2) == "carrier" else "noncarrier")
            parts.setdefault(key, []).append(load_file(f))
    groups = {}
    for key, ds in parts.items():
        d = pd.concat(ds, ignore_index=True)
        if d.seed.duplicated().any():
            sys.exit(f"duplicate seeds in group {key}: check the batch files")
        groups[key] = d
    return groups


def hr_trend(d0, d1):
    """Poisson GLM on 10-y bands: events ~ C(band) + carrier + carrier:dec, offset log person-time."""
    rows = []
    for arm, d in ((0, d0), (1, d1)):
        age, tum = d["age"].values, d["fate"].values == "tumour"
        for lo, hi in zip(EDGES[:-1], EDGES[1:]):
            pt = np.clip(np.minimum(age, hi) - lo, 0, None).sum()
            if pt > 0:
                rows.append(dict(arm=arm, band=lo, mid=lo + 5, pt=pt, ev=((age >= lo) & (age < hi) & tum).sum()))
    t = pd.DataFrame(rows)
    t = t[t.groupby("band").ev.transform("sum") > 0].copy()
    t["dec"] = (t.mid - 60) / 10
    m = smf.glm("ev ~ C(band) + arm + arm:dec", t, family=sm.families.Poisson(), offset=np.log(t.pt)).fit()
    ci = m.conf_int()
    piv = t.pivot(index="mid", columns="arm", values=["ev", "pt"])
    curve = []
    for mid, r in piv.iterrows():
        e0, e1 = r[("ev", 0)], r[("ev", 1)]
        if e0 > 0 and e1 > 0:
            hr = (e1 / r[("pt", 1)]) / (e0 / r[("pt", 0)]); se = np.sqrt(1 / e0 + 1 / e1)
            curve.append((mid, hr, hr * np.exp(-1.96 * se), hr * np.exp(1.96 * se)))
    return dict(HR_at_60=np.exp(m.params["arm"]), HR_lo=np.exp(ci.loc["arm", 0]), HR_hi=np.exp(ci.loc["arm", 1]),
                trend=np.exp(m.params["arm:dec"]), trend_lo=np.exp(ci.loc["arm:dec", 0]),
                trend_hi=np.exp(ci.loc["arm:dec", 1])), curve


def analyse(label, d0, d1):
    a0, t0, _ = at_endpoint(d0, T_END); a1, t1, _ = at_endpoint(d1, T_END)
    s = pair_stats(a0, t0, a1, t1)
    ci = boot_pair(a0, t0, a1, t1, ["dmeda_y", "residual_y"])
    tr, curve = hr_trend(d0, d1)
    row = dict(level=label, n0=len(d0), n1=len(d1), events0=int(t0.sum()), events1=int(t1.sum()))
    for arm, t in (("0", t0), ("1", t1)):
        k, n = t.sum(), len(t); lo, hi = wilson(k, n)
        row[f"p{arm}"], row[f"p{arm}_lo"], row[f"p{arm}_hi"] = k / n, lo, hi
    row.update(meda0=s["meda0"], meda1=s["meda1"], obs_shift=s["dmeda_y"],
               obs_lo=ci["dmeda_y"][0], obs_hi=ci["dmeda_y"][1],
               ph_expected_shift=s["ph_pred_shift"], residual=s["residual_y"],
               res_lo=ci["residual_y"][0], res_hi=ci["residual_y"][1], **tr)
    plot(label, d0, d1, curve, row)
    return row


def plot(label, d0, d1, curve, row):
    FIG.mkdir(exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw=dict(wspace=0.35))
    grid = np.linspace(0, T_END, 400)
    for d, col, lab in ((d0, NONCARRIER, "non-carrier"), (d1, CARRIER, "carrier")):
        tum_age = np.sort(d.loc[d.fate == "tumour", "age"].values)
        ax1.plot(grid, np.searchsorted(tum_age, grid, side="right") / len(d), color=col, lw=2, label=lab)
    for med, col in ((row["meda0"], NONCARRIER), (row["meda1"], CARRIER)):
        ax1.axvline(med, color=col, lw=1, ls=":")
    ax1.set_xlabel("Age (years)"); ax1.set_ylabel("Cumulative occurrence"); ax1.set_xlim(0, T_END)
    ax1.legend(frameon=False, fontsize=7, loc="upper left")
    ax1.set_title(f"{label}: MEDA {row['meda0']:.1f} → {row['meda1']:.1f} y", fontsize=8, loc="left")
    if curve:
        x, y, lo, hi = map(np.array, zip(*curve))
        ax2.errorbar(x, y, yerr=[y - lo, hi - y], fmt="o-", color=CARRIER, ms=4, lw=1.5, elinewidth=1)
    ax2.set_yscale("log"); ax2.set_xlabel("Age (years)"); ax2.set_ylabel("Hazard ratio, carrier / non-carrier")
    ax2.set_title(f"change per decade {row['trend']:.2f} ({row['trend_lo']:.2f}–{row['trend_hi']:.2f})", fontsize=8, loc="left")
    for ax in (ax1, ax2):
        ax.spines[["top", "right"]].set_visible(False); ax.grid(True, color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)
    fig.savefig(FIG / f"germline_{label}.png", dpi=180, bbox_inches="tight"); plt.close(fig)


def synthetic(k, p_target, n=1000, scale=1.0):
    """Event times with cumulative hazard H(t) = c t^k, c set so P(event by 92) = p_target (times scale)."""
    c = -np.log(1 - p_target) / T_END ** k * scale
    t = (RNG.exponential(size=n) / c) ** (1 / k)
    age = np.minimum(t, T_END)
    return pd.DataFrame(dict(seed=np.arange(1, n + 1), gen=np.round(age * WEEKS_PER_YEAR), pop=0, age=age,
                             fate=np.where(t <= T_END, "tumour", "none")))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=str(GERMLINE)); ap.add_argument("--synthetic", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    rows = []
    if a.synthetic:
        base = synthetic(5, 0.05)
        rows.append(analyse("SYNTH_multistage_k5_vs_k4", base, synthetic(4, 0.30)))
        # PH control: same shape, hazard x 7 (gives ~30% occurrence)
        rows.append(analyse("SYNTH_PH_control_x7", base, synthetic(5, 0.05, scale=7.0)))
        out = OUT / "germline_SYNTHETIC_validation.csv"
    else:
        groups = load_groups(a.dir)
        carr = {lv: d for (lv, arm), d in groups.items() if arm == "carrier"}
        nonc = {lv: d for (lv, arm), d in groups.items() if arm == "noncarrier"}
        print("runs per group:", {f"{lv}_{arm}": len(d) for (lv, arm), d in sorted(groups.items())})
        if not carr:
            sys.exit(f"no germ_*_carrier files in {a.dir}")
        for lvl in sorted(carr):
            r = analyse(f"germline_{lvl}", nonc[lvl], carr[lvl]); r["contrast"] = "germline first hit"; rows.append(r)
        # Comparator in the same model: the same baseline, but occurrence raised by more mutation input
        # instead of a germline hit (non-carriers at the next input level). The question for the paper is
        # whether the germline hit gives a lower trend / more negative residual than this.
        lv = sorted(nonc)
        for l0, l1 in zip(lv[:-1], lv[1:]):
            r = analyse(f"input_{l0}_to_{l1}", nonc[l0], nonc[l1]); r["contrast"] = "more mutation input"; rows.append(r)
        out = OUT / "germline_summary.csv"
    r = pd.DataFrame(rows); r.to_csv(out, index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    if "contrast" not in r:
        r["contrast"] = "synthetic"
    cols = ["contrast", "level", "events0", "events1", "p0", "p1", "meda0", "meda1", "obs_shift", "ph_expected_shift",
            "residual", "res_lo", "res_hi", "HR_at_60", "trend", "trend_lo", "trend_hi"]
    print(r[cols].round(3).to_string(index=False)); print("saved", out)


if __name__ == "__main__":
    main()
