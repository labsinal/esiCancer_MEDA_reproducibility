"""
How many runs per arm does the germline experiment need?

Replicates the synthetic validation (germline_analysis.py --synthetic) many times:
  multistage: non-carriers k = 5 stages, carriers k = 4 (theory: HR falls as 1/t, ~0.85 per decade at 60)
  PH control: carriers = non-carrier hazard x constant (truth: trend 1, residual 0)
for baseline occurrence p0 (non-carriers, by 92 y) and runs per arm n.
Carriers are calibrated to the same occurrence in both scenarios, so the only difference is the age pattern.
Reports: mean trend, share of replicates whose 95% CI for trend excludes 1 (power, or false-positive rate
for the PH control), and mean / SD of the MEDA residual.
"""
import numpy as np
import pandas as pd

from germline_analysis import synthetic, hr_trend, T_END
import germline_analysis as ga
import sys
from ph_benchmark import at_endpoint, pair_stats  # noqa: E402

REPS = 150
rows = []
for p0, p1 in ((0.05, 0.30), (0.20, 0.60)):
    # carriers' occurrence p1 is the same in both scenarios; PH scale chosen to reach p1
    ph_scale = np.log(1 - p1) / np.log(1 - p0)
    for n in (1000, 3000):
        for scen in ("multistage", "PH"):
            tr, sig, res = [], 0, []
            for _ in range(REPS):
                d0 = synthetic(5, p0, n)
                d1 = synthetic(4, p1, n) if scen == "multistage" else synthetic(5, p0, n, scale=ph_scale)
                t, _ = hr_trend(d0, d1)
                tr.append(t["trend"]); sig += (t["trend_hi"] < 1) or (t["trend_lo"] > 1)
                a0, e0, _ = at_endpoint(d0, T_END); a1, e1, _ = at_endpoint(d1, T_END)
                res.append(pair_stats(a0, e0, a1, e1)["residual_y"])
            rows.append(dict(p0=p0, p1=p1, n_per_arm=n, scenario=scen, mean_trend=np.mean(tr),
                             share_trend_CI_excludes_1=sig / REPS, mean_residual=np.mean(res), sd_residual=np.std(res),
                             residual_power_approx=np.mean(np.abs(np.mean(res)) / np.std(res) > 1.96) if scen else None))
r = pd.DataFrame(rows).drop(columns="residual_power_approx")
r["residual_z"] = r.mean_residual / r.sd_residual
r.to_csv(ga.OUT / "germline_power_check.csv", index=False)
pd.set_option("display.width", 200)
print(r.round(3).to_string(index=False))
