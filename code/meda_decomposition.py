"""
Links the MEDA analysis to the age-trend analysis.

For each Figure 3 series, lowest input vs every higher input (T = 92 y):
  observed MEDA shift      = MEDA1 - MEDA0
  PH-expected shift        = PH-predicted MEDA1 - MEDA0   (what a constant multiplier, i.e. a flat
                             relative risk, produces through depletion of the at-risk pool)
  residual                 = observed - PH-predicted       (the part due to the RR changing with age)
and the change in hazard ratio per decade (sim_age_trend_rr.py).

If the two approaches describe the same thing, the residual should be ~0 when the trend is ~1,
negative when the trend is < 1 (risk concentrated young) and positive when > 1.
"""
import numpy as np
import pandas as pd
from scipy import stats

from ph_benchmark import SERIES, load, at_endpoint, pair_stats, boot_pair, OUT

trend = pd.read_csv(OUT / "fig3_HR_trend_summary.csv")
rows = []
for series, (folder, files) in SERIES.items():
    data = {lvl: load(folder, pfx) for pfx, lvl in files.items()}
    levels = sorted(data)
    a0, t0, _ = at_endpoint(data[levels[0]], 92)
    for lvl in levels[1:]:
        a1, t1, _ = at_endpoint(data[lvl], 92)
        s = pair_stats(a0, t0, a1, t1)
        ci = boot_pair(a0, t0, a1, t1, ["residual_y"])["residual_y"]
        tr = trend[(trend.series == series) & (trend.input1 == lvl)].iloc[0]
        rows.append(dict(series=series, input0=levels[0], input1=lvl, p0=s["p0"], p1=s["p1"],
                         obs_shift=s["dmeda_y"], ph_expected_shift=s["ph_pred_shift"],
                         residual=s["residual_y"], res_lo=ci[0], res_hi=ci[1],
                         trend_per_decade=tr.trend_per_decade, trend_lo=tr.trend_lo, trend_hi=tr.trend_hi))
r = pd.DataFrame(rows)
r.to_csv(OUT / "meda_decomposition.csv", index=False)
pd.set_option("display.width", 220)
print(r.round(2).to_string(index=False))
ok = r.dropna()
rho, p = stats.spearmanr(np.log(ok.trend_per_decade), ok.residual)
print(f"\nSpearman(log trend, residual) = {rho:.2f}, p = {p:.3g}, n = {len(ok)} pairs (pairs share baselines, so p is optimistic)")
print(f"Share of observed MEDA shift accounted for by PH expectation, all pairs pooled: "
      f"{ok.ph_expected_shift.sum() / ok.obs_shift.sum():.2f}")
