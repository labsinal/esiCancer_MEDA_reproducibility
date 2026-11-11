"""
Reanalysis of the Figure 3 simulations for the response to reviewers (Manuscript 1790300).

Addresses Reviewer comments 2, 3 and 10:
  - occurrence probability with Wilson 95% CI (binomial),
  - MEDA = median(tau | tau <= T) with bootstrap 95% CI (resampling runs),
  - absolute differences and direction-aligned log changes, log(rInc) and -log(rMEDA),
  - proportional-hazards (PH) null benchmark: the MEDA expected for the higher-input
    condition if its tumour hazard were the baseline hazard multiplied by a constant,
  - fate classification (tumour / extinction / no event by T),
  - sensitivity of all of the above to the observation endpoint T.

Input: esiCancer *_automatics_ancestralResults.csv files (one row per seed; the run stops at
the generation the population reaches the esiTumor threshold, or at generation 4799).

PH benchmark (with extinction as a competing event)
  F0(t)  = empirical cumulative incidence of esiTumor in the baseline condition
           (all censoring happens at T, so the empirical proportion is the Aalen-Johansen estimate).
  Under PH on the subdistribution hazard, F_a(t) = 1 - (1 - F0(t))**a.
  a is chosen so that F_a(T) equals the observed occurrence p1 of the higher-input condition:
           a = log(1 - p1) / log(1 - p0).
  The predicted conditional median is the age where F_a(t) = p1 / 2, i.e. the age at which
  F0 reaches q* = 1 - (1 - p1/2)**(1/a)  (closed form, always <= p0).
  Residual = observed MEDA1 - PH-predicted MEDA1. Residual ~ 0 means the MEDA shift is what
  the outcome definition and hazard scaling alone would produce.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from paths import SIM as DATA, OUT
WEEKS_PER_YEAR = 52.143          # same conversion as derivadas.ipynb (4800 gen ~ 92 y)
THRESHOLD = 10_000               # esiTumor = 10x the initial 1,000 cells
N_BOOT = 2000
RNG = np.random.default_rng(20261015)

# Figure 3 panel -> folder, file prefix. The numeric input level is parsed from the file name.
# Meaning of the input level per series must be confirmed by the authors (see follow-up doc).
SERIES = {
    "A_constant":          ("revision/constant_events_escc_no_aging", {"constant_1": 1, "constant_3": 3, "constant_5": 5, "constant_10": 10}),
    "B_increasing":        ("revision/linear_events_escc_no_aging",   {"linear_5": 5, "linear_10": 10, "linear_14": 14, "linear_20": 20}),
    "C_multistage":        ("revision/6_genes_interaction_mut_rate",  {"6_genes_mut_002": 2, "6_genes_mut_003": 3, "6_genes_mut_004": 4, "6_genes_mut_006": 6}),
    "D_constant_aging":    ("revision/constant_events_escc_aging",               {"mut_rate_const_001": 1, "mut_rate_const_003": 3, "mut_rate_const_005": 5, "mut_rate_const_010": 10}),
    "E_incr_aging_ctor001": ("revision/mut_rate_PA_mut",              {"mut_rate_PA5_mut_2": 5, "mut_rate_PA10_mut": 10, "mut_rate_PA14_mut": 14, "mut_rate_PA20_mut": 20}),
    "G_incr_aging_ctor002": ("revision/ctor_002",                     {"ctor_002_5": 5, "ctor_002_10": 10, "ctor_002_14": 14, "ctor_002_20": 20}),
    "H_incr_aging_ctor02":  ("revision/ctor_02",                      {"ctor_02_5": 5, "ctor_02_10": 10, "ctor_02_14": 14, "ctor_02_20": 20}),
}


def load(folder, prefix):
    f = DATA / folder / f"{prefix}_automatics_ancestralResults.csv"
    d = pd.read_csv(f, sep=";", usecols=[0, 1, 3])
    d.columns = ["seed", "gen", "pop"]
    d["age"] = d["gen"] / WEEKS_PER_YEAR
    d["fate"] = np.where(d["pop"] >= THRESHOLD, "tumour", np.where(d["pop"] == 0, "extinct", "none"))
    return d


def at_endpoint(d, T):
    """Re-define outcomes as if the simulation had stopped at age T."""
    tum = (d["fate"].values == "tumour") & (d["age"].values <= T)
    ext = (d["fate"].values == "extinct") & (d["age"].values <= T)
    return d["age"].values, tum, ext


def wilson(k, n, z=1.96):
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return c - h, c + h


def ph_predicted_median(age0, tum0, p1):
    """Age at which baseline cumulative incidence reaches q* (see module docstring)."""
    n0 = len(age0)
    p0 = tum0.sum() / n0
    if not (0 < p0 < 1) or not (0 < p1 < 1):
        return np.nan
    a = np.log(1 - p1) / np.log(1 - p0)
    q = 1 - (1 - p1 / 2) ** (1 / a)
    k = q * n0                          # number of baseline events needed
    t = np.sort(age0[tum0])
    # linear interpolation on the empirical CIF step function
    idx = np.clip(k - 1, 0, len(t) - 1)
    lo = int(np.floor(idx)); hi = min(lo + 1, len(t) - 1)
    return t[lo] + (idx - lo) * (t[hi] - t[lo])


def pair_stats(age0, tum0, age1, tum1):
    n0, n1 = len(age0), len(age1)
    p0, p1 = tum0.mean(), tum1.mean()
    m0 = np.median(age0[tum0]) if tum0.any() else np.nan
    m1 = np.median(age1[tum1]) if tum1.any() else np.nan
    pred = ph_predicted_median(age0, tum0, p1)
    return dict(
        p0=p0, p1=p1, dp_pp=100 * (p1 - p0),
        meda0=m0, meda1=m1, dmeda_y=m1 - m0,
        log_rInc=np.log(p1 / p0) if p0 > 0 and p1 > 0 else np.nan,
        neg_log_rMEDA=-np.log(m1 / m0),
        Q_original=(p1 / p0) / (m1 / m0) if p0 > 0 else np.nan,
        ph_pred_meda1=pred,
        ph_pred_shift=pred - m0,
        residual_y=m1 - pred,
    )


def boot_pair(age0, tum0, age1, tum1, keys):
    n0, n1 = len(age0), len(age1)
    out = {k: [] for k in keys}
    for _ in range(N_BOOT):
        i0 = RNG.integers(0, n0, n0); i1 = RNG.integers(0, n1, n1)
        s = pair_stats(age0[i0], tum0[i0], age1[i1], tum1[i1])
        for k in keys:
            out[k].append(s[k])
    return {k: np.nanpercentile(v, [2.5, 97.5]) for k, v in out.items()}


def main():
    OUT.mkdir(exist_ok=True)
    cond_rows, pair_rows = [], []
    for T in (60, 75, 92):
        for series, (folder, files) in SERIES.items():
            loaded = {lvl: load(folder, pfx) for pfx, lvl in files.items()}
            levels = sorted(loaded)
            for lvl in levels:
                age, tum, ext = at_endpoint(loaded[lvl], T)
                n, k = len(age), tum.sum()
                lo, hi = wilson(k, n)
                bm = [np.median(age[tum][RNG.integers(0, k, k)]) for _ in range(N_BOOT)] if k > 1 else [np.nan]
                cond_rows.append(dict(T=T, series=series, input=lvl, n_runs=n, n_tumour=int(k),
                                      n_extinct=int(ext.sum()), n_none=int(n - k - ext.sum()),
                                      p_occ=k / n, p_lo=lo, p_hi=hi,
                                      meda=np.median(age[tum]) if k else np.nan,
                                      meda_lo=np.nanpercentile(bm, 2.5), meda_hi=np.nanpercentile(bm, 97.5)))
            # consecutive pairs (as in Fig 3F/3I) plus lowest -> highest
            pairs = list(zip(levels[:-1], levels[1:])) + [(levels[0], levels[-1])]
            for l0, l1 in pairs:
                a0, t0, _ = at_endpoint(loaded[l0], T)
                a1, t1, _ = at_endpoint(loaded[l1], T)
                if t0.sum() < 5 or t1.sum() < 5:
                    continue
                s = pair_stats(a0, t0, a1, t1)
                ci = boot_pair(a0, t0, a1, t1, ["dp_pp", "dmeda_y", "log_rInc", "neg_log_rMEDA", "residual_y", "ph_pred_meda1"])
                row = dict(T=T, series=series, input0=l0, input1=l1, log_input_ratio=np.log(l1 / l0), **s)
                for k, (lo, hi) in ci.items():
                    row[f"{k}_lo"], row[f"{k}_hi"] = lo, hi
                row["elast_inc"] = s["log_rInc"] / row["log_input_ratio"]
                row["elast_meda"] = s["neg_log_rMEDA"] / row["log_input_ratio"]
                row["share_of_meda_shift_explained_by_PH"] = s["ph_pred_shift"] / s["dmeda_y"] if s["dmeda_y"] != 0 else np.nan
                pair_rows.append(row)
    cond = pd.DataFrame(cond_rows)
    pairs = pd.DataFrame(pair_rows)
    cond.to_csv(OUT / "conditions_occurrence_meda.csv", index=False)
    pairs.to_csv(OUT / "pairs_ph_benchmark.csv", index=False)

    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40); pd.set_option("display.max_rows", 400)
    r = lambda x: x.round(3)
    print("=== Conditions at T = 92 y ===")
    print(r(cond[cond["T"] == 92].drop(columns="T")).to_string(index=False))
    cols = ["series", "input0", "input1", "p0", "p1", "meda0", "meda1", "dmeda_y", "dmeda_y_lo", "dmeda_y_hi",
            "ph_pred_meda1", "residual_y", "residual_y_lo", "residual_y_hi", "share_of_meda_shift_explained_by_PH",
            "log_rInc", "neg_log_rMEDA", "elast_inc", "elast_meda", "Q_original"]
    for T in (92, 75, 60):
        print(f"\n=== Pairs, PH benchmark, T = {T} y ===")
        print(r(pairs[pairs["T"] == T][cols]).to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
