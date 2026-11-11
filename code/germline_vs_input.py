"""
Direct test: does a germline first hit differ from a mutation-input increase that raises occurrence by a
similar amount (same model, same baseline)?

Groups (6-gene multistage model; all seed batches of each group are pooled, see load_groups):
  base  = non-carriers, 3 events/division   (germ_003_non_carrier)
  germ  = carriers,     3 events/division   (germ_003_carrier; gene A in 100% of cells at generation 0)
  input = non-carriers, 4 events/division   (Fig 3C run 6_genes_mut_004, reproduced exactly at 003)

1. Age trend of the hazard ratio: one Poisson GLM on 10-y bands with all three groups,
   events ~ C(band) + group + group:dec, offset log(person-time), base = reference.
   Difference of trends = exp(beta_germ:dec - beta_input:dec), Wald CI from the joint covariance.
2. MEDA residual (observed minus PH-expected), each vs base, and their difference, with a bootstrap that
   resamples all three groups jointly (so the shared baseline is handled correctly).
3. MEDA of germ vs input directly (similar occurrence, different mechanism).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

from paths import GERMLINE, OUT
from ph_benchmark import at_endpoint, pair_stats  # noqa: E402
from germline_analysis import load_groups, EDGES, T_END  # noqa: E402

B = 2000
RNG = np.random.default_rng(1971)
_g = load_groups(GERMLINE)
G = {"base": _g[("003", "noncarrier")], "germ": _g[("003", "carrier")], "input": _g[("004", "noncarrier")]}
print("runs per group:", {k: len(v) for k, v in G.items()})

# ---- 1. joint Poisson model ----
rows = []
for g, d in G.items():
    age, tum = d.age.values, d.fate.values == "tumour"
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        pt = np.clip(np.minimum(age, hi) - lo, 0, None).sum()
        if pt > 0:
            rows.append(dict(group=g, band=lo, dec=(lo + 5 - 60) / 10, pt=pt, ev=((age >= lo) & (age < hi) & tum).sum()))
t = pd.DataFrame(rows)
t = t[t.groupby("band").ev.transform("sum") > 0]
t["germ"] = (t.group == "germ").astype(int); t["inp"] = (t.group == "input").astype(int)
m = smf.glm("ev ~ C(band) + germ + inp + germ:dec + inp:dec", t, family=sm.families.Poisson(), offset=np.log(t.pt)).fit()
V = m.cov_params()
d_trend = m.params["germ:dec"] - m.params["inp:dec"]
se = np.sqrt(V.loc["germ:dec", "germ:dec"] + V.loc["inp:dec", "inp:dec"] - 2 * V.loc["germ:dec", "inp:dec"])
z = d_trend / se
from scipy import stats  # noqa: E402
print("=== 1. Age trend of the hazard ratio (per decade), joint model, reference = non-carriers 3 ev/div ===")
for k, lab in (("germ:dec", "germline first hit"), ("inp:dec", "more input (4 ev/div)")):
    lo, hi = m.conf_int().loc[k]
    print(f"  {lab:24s} {np.exp(m.params[k]):.3f} ({np.exp(lo):.3f}-{np.exp(hi):.3f})")
print(f"  ratio germline / input   {np.exp(d_trend):.3f} ({np.exp(d_trend - 1.96 * se):.3f}-{np.exp(d_trend + 1.96 * se):.3f}), "
      f"z = {z:.2f}, p = {2 * stats.norm.sf(abs(z)):.3f}")


# ---- 2/3. residuals and direct MEDA comparison, joint bootstrap ----
def stats_of(base, germ, inp):
    ab, tb, _ = at_endpoint(base, T_END); ag, tg, _ = at_endpoint(germ, T_END); ai, ti, _ = at_endpoint(inp, T_END)
    sg, si = pair_stats(ab, tb, ag, tg), pair_stats(ab, tb, ai, ti)
    return dict(res_germ=sg["residual_y"], res_input=si["residual_y"], res_diff=sg["residual_y"] - si["residual_y"],
                meda_germ=sg["meda1"], meda_input=si["meda1"], meda_diff=sg["meda1"] - si["meda1"],
                p_germ=sg["p1"], p_input=si["p1"])


point = stats_of(*G.values())
bs = []
for _ in range(B):
    bs.append(stats_of(*[d.iloc[RNG.integers(0, len(d), len(d))] for d in G.values()]))
bs = pd.DataFrame(bs)
print("\n=== 2. MEDA residual (observed - PH-expected), years; joint bootstrap ===")
for k, lab in (("res_germ", "germline first hit"), ("res_input", "more input"), ("res_diff", "difference germ - input")):
    lo, hi = np.nanpercentile(bs[k], [2.5, 97.5])
    print(f"  {lab:24s} {point[k]:+.2f} ({lo:+.2f} to {hi:+.2f})   P(boot >= 0) = {(bs[k] >= 0).mean():.3f}")
print("\n=== 3. Similar occurrence, different mechanism ===")
lo, hi = np.nanpercentile(bs["meda_diff"], [2.5, 97.5])
print(f"  occurrence: germline {point['p_germ']:.3f} vs input {point['p_input']:.3f}")
print(f"  MEDA: germline {point['meda_germ']:.2f} y vs input {point['meda_input']:.2f} y; "
      f"difference {point['meda_diff']:+.2f} ({lo:+.2f} to {hi:+.2f})")
pd.DataFrame([point]).to_csv(OUT / "germline" / "germline_vs_input.csv", index=False)
