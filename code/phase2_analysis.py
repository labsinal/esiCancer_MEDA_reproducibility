"""
Phase 2 controls (see RUN_PROTOCOL.md): verification that the esiTable mutation-rate field is inert (Part A),
the pilot check of fractional events per division, and matched-input turnover (Part B).

Usage (from this folder):
  python phase2_analysis.py pilot      # genome-size scaling: mutated cells at 10 ev/div (genome x1, x2) and 5 ev/div
  python phase2_analysis.py verify     # verify_mutator10 vs germline/germ_003_carrier, seeds 1-50
  python phase2_analysis.py exp2       # matched-input turnover (CTOR 0.01 arms = Fig 3A runs)
  options: --runs DIR (default runs/)

File names (in --runs): verify_mutator10, ctor<code>_e<ev>[_g<k>], pilot_e<ev>[_g<k>], each followed by
_automatics_ancestralResults.csv (and _genesMutationResults.csv for the pilot).
<code> 001 = 0.01, 002 = 0.02, 02 = 0.2; _g<k> = esiTable with genome size multiplied by k, which divides the
per-gene mutation probability by k (effective events per division = ev / k). esiCancer accepts only integer events
per division, so matched input at high turnover is obtained with a larger genome instead of fractional events.
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
from ph_benchmark import WEEKS_PER_YEAR, THRESHOLD, at_endpoint, pair_stats, boot_pair, wilson  # noqa: E402
from germline_analysis import hr_trend  # noqa: E402

from paths import SIM, GERMLINE, TURNOVER, OUT as _OUT  # noqa: E402
DATA = SIM / "revision" / "constant_events_escc_no_aging"
FIG3A = {1: "constant_1", 3: "constant_3", 5: "constant_5", 10: "constant_10"}
OUT = _OUT / "turnover"
T_END = 92.0
B = 2000
RNG = np.random.default_rng(2024)


def ev_value(tag):
    return float(tag.replace("p", "."))


def ctor_value(code):
    return int(code) / 10 ** (len(code) - 1)


def load(path):
    d = pd.read_csv(path, sep=";", usecols=[0, 1, 2, 3])
    d.columns = ["seed", "gen", "ndiv", "pop"]
    d["age"] = d["gen"] / WEEKS_PER_YEAR
    d["fate"] = np.where(d["pop"] >= THRESHOLD, "tumour", np.where(d["pop"] == 0, "extinct", "none"))
    return d.sort_values("seed").reset_index(drop=True)


def find(runs, pattern):
    out = {}
    for f in sorted(Path(runs).glob("*_automatics_ancestralResults.csv")):
        m = re.fullmatch(pattern + r"_automatics_ancestralResults\.csv", f.name)
        if m:
            out[m.groups()] = f
    return out


def describe(d):
    a, t, e = at_endpoint(d, T_END)
    k, n = t.sum(), len(t); lo, hi = wilson(k, n)
    bm = [np.median(a[t][RNG.integers(0, k, k)]) for _ in range(B)] if k > 1 else [np.nan]
    return dict(n=n, tumours=int(k), extinct=int(e.sum()), p_occ=k / n, p_lo=lo, p_hi=hi,
                meda=np.median(a[t]) if k else np.nan, meda_lo=np.nanpercentile(bm, 2.5), meda_hi=np.nanpercentile(bm, 97.5))


def compare(d0, d1):
    """d1 vs d0: occurrence and MEDA differences, PH residual, HR at 60 and its trend per decade."""
    a0, t0, _ = at_endpoint(d0, T_END); a1, t1, _ = at_endpoint(d1, T_END)
    s = pair_stats(a0, t0, a1, t1)
    ci = boot_pair(a0, t0, a1, t1, ["dp_pp", "dmeda_y", "residual_y"])
    tr, _ = hr_trend(d0, d1)
    return dict(dp_pp=s["dp_pp"], dp_lo=ci["dp_pp"][0], dp_hi=ci["dp_pp"][1],
                dmeda=s["dmeda_y"], dmeda_lo=ci["dmeda_y"][0], dmeda_hi=ci["dmeda_y"][1],
                ph_expected=s["ph_pred_shift"], residual=s["residual_y"], res_lo=ci["residual_y"][0],
                res_hi=ci["residual_y"][1], HR_at_60=tr["HR_at_60"], trend=tr["trend"],
                trend_lo=tr["trend_lo"], trend_hi=tr["trend_hi"])


def show(df, title):
    pd.set_option("display.width", 260); pd.set_option("display.max_columns", 40)
    print(f"\n=== {title} ===")
    print(df.round(3).to_string(index=False))


# ------------------------------------------------------------------ pilot
def pilot(runs):
    res = {}
    for tag in ("e10", "e10_g2", "e5"):
        f = Path(runs) / f"pilot_{tag}_automatics_genesMutationResults.csv"
        if not f.exists():
            sys.exit(f"missing {f.name}")
        g = pd.read_csv(f, sep=";"); g.columns = [c.strip() for c in g.columns]
        cols = [c for c in g.columns if c.endswith("monoAl") or c.endswith("biAl")]
        alive = g["Total Population"] > 0
        res[tag] = (g.loc[alive, cols].sum(axis=1) / g.loc[alive, "Total Population"]).mean()
    rg, r5 = res["e10_g2"] / res["e10"], res["e5"] / res["e10"]
    print("mean mutated gene-events per cell at end: "
          f"10 ev/div = {res['e10']:.4f}; 10 ev/div genome x2 = {res['e10_g2']:.4f}; 5 ev/div = {res['e5']:.4f}")
    print(f"ratios to 10 ev/div: genome x2 = {rg:.2f}, 5 ev/div = {r5:.2f} (both expected about 0.5)")
    if 0.3 <= rg <= 0.7 and 0.3 <= r5 <= 0.7:
        print("-> genome size scales the mutation input as expected. Proceed with the matched design (_g2, _g20).")
    elif rg > 0.85:
        print("-> genome size does not change the input: use Plan B (integer events, higher input).")
    else:
        print("-> unclear: send me the files before deciding.")


# ------------------------------------------------------------------ part A: mutation-rate field
def verify(runs):
    f = Path(runs) / "verify_mutator10_automatics_ancestralResults.csv"
    if not f.exists():
        sys.exit(f"missing {f.name}")
    test = load(f)
    ref = load(GERMLINE / "germ_003_carrier_automatics_ancestralResults.csv")
    m = test.merge(ref, on="seed", suffixes=("_test", "_ref"))
    same = (m.gen_test == m.gen_ref) & (m.pop_test == m.pop_ref) & (m.ndiv_test == m.ndiv_ref)
    print(f"seeds compared: {len(m)} (expected 50)")
    print(f"identical generation, population and number of divisions: {same.sum()} of {len(m)}")
    if same.all():
        print("-> the mutation-rate field (A: +10 events/division) has NO effect: outputs are identical seed by seed.")
    else:
        print("-> outputs DIFFER: the field is active. Tumours by 92 y: test "
              f"{(m.pop_test >= THRESHOLD).sum()} vs reference {(m.pop_ref >= THRESHOLD).sum()}; "
              f"median divisions {m.ndiv_test.median():.0f} vs {m.ndiv_ref.median():.0f}. Stop and check before Part B.")


# ------------------------------------------------------------------ experiment 2
def exp2(runs):
    found = find(runs, r"ctor_?(\d+)_e([0-9p]+)(?:_g(\d+))?")
    if not found:
        sys.exit("no ctor*_e* runs found")
    conds = {}
    for (code, tag, g), f in found.items():
        c, ev = ctor_value(code), ev_value(tag) / int(g or 1)     # effective events per division
        conds[(c, round(ev, 6))] = load(f)
    # CTOR 0.01 arms: the existing Fig 3A runs (same table, no aging; max divisions never binds at CTOR 0.01),
    # unless new ctor001_* runs are present
    for ev, name in FIG3A.items():
        if (0.01, float(ev)) not in conds and (DATA / f"{name}_automatics_ancestralResults.csv").exists():
            conds[(0.01, float(ev))] = load(DATA / f"{name}_automatics_ancestralResults.csv")
    rows = []
    for (c, ev), d in sorted(conds.items()):
        # realised input = total events / total years at risk (ratio of sums, all runs), and the same in
        # tumour-free runs only, where clonal expansion does not inflate the number of divisions
        expected = c * ev * 1000 * WEEKS_PER_YEAR            # 1000 cells x divisions/cell/generation x events x generations/year
        realised = (d["ndiv"] * ev).sum() / d["age"].sum()
        nt = d[d.fate != "tumour"]
        realised_nt = (nt["ndiv"] * ev).sum() / nt["age"].sum() if len(nt) else np.nan
        rows.append(dict(ctor=c, events_per_div=ev, input_per_cell_gen=round(c * ev, 6), **describe(d),
                         expected_events_per_year=expected, realised_all_over_expected=realised / expected,
                         realised_tumourfree_over_expected=realised_nt / expected))
    cond = pd.DataFrame(rows)
    # matched comparisons: each CTOR vs CTOR 0.01 at the same input per cell per generation
    cmp_rows = []
    ref = {round(c * ev, 6): (c, ev) for (c, ev) in conds if c == 0.01}
    for (c, ev), d in sorted(conds.items()):
        key = round(c * ev, 6)
        if c != 0.01 and key in ref:
            cmp_rows.append(dict(comparison=f"CTOR {c:g} vs 0.01, matched input {key:g}", **compare(conds[ref[key]], d)))
    # supply / dynamics split for unmatched arms: CTOR c at the same events/division as the reference
    split = []
    for (c, ev), d in sorted(conds.items()):
        if c == 0.01 or (0.01, ev) not in conds:
            continue
        ref_d = conds[(0.01, ev)]
        matched = conds.get((c, round(0.01 * ev / c, 6)))
        if matched is None:
            continue
        p = lambda x: (x.fate == "tumour").mean()
        med = lambda x: x.loc[x.fate == "tumour", "age"].median()
        split.append(dict(ctor=c, events_per_div=ev,
                          total_logRR_occ=np.log(p(d) / p(ref_d)), dynamics_logRR=np.log(p(matched) / p(ref_d)),
                          supply_logRR=np.log(p(d) / p(matched)),
                          total_dMEDA=med(d) - med(ref_d), dynamics_dMEDA=med(matched) - med(ref_d),
                          supply_dMEDA=med(d) - med(matched)))
    OUT.mkdir(exist_ok=True)
    cond.to_csv(OUT / "exp2_conditions.csv", index=False); show(cond, "Experiment 2: conditions")
    if cmp_rows:
        c2 = pd.DataFrame(cmp_rows); c2.to_csv(OUT / "exp2_matched_comparisons.csv", index=False)
        show(c2, "Experiment 2: matched input, each CTOR vs 0.01")
    if split:
        s2 = pd.DataFrame(split); s2.to_csv(OUT / "exp2_supply_vs_dynamics.csv", index=False)
        show(s2, "Experiment 2: total turnover effect = population dynamics (matched) + mutation supply")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["verify", "pilot", "exp2"])
    ap.add_argument("--runs", default=str(TURNOVER))
    a = ap.parse_args()
    {"verify": verify, "pilot": pilot, "exp2": exp2}[a.what](a.runs)
