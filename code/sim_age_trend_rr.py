"""
Same measure as age_trend_rr.py, applied to the Figure 3 simulations, so that simulated and
epidemiological examples sit on the same axes:

  hazard ratio (HR) of esiTumor formation, higher vs lowest input of each series, by age, and
  its trend per decade of age.

Person-time: each run contributes time at risk from age 0 until esiTumor, extinction or the end of
simulation (T = 92 y). Events and person-time are tallied in 10-year age bands.
Model: Poisson GLM, events ~ C(band) + high + high:dec, offset log(person-time),
where dec = (band midpoint - 60) / 10. exp(high) = HR at 60 y; exp(high:dec) = change per decade.
Extinction is treated as leaving the risk set (cause-specific hazard).
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

from ph_benchmark import SERIES, load, OUT

EDGES = np.arange(0, 100, 10.0)


def tally(d, label):
    age = d["age"].values
    tum = d["fate"].values == "tumour"
    out = []
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        pt = np.clip(np.minimum(age, hi) - lo, 0, None).sum()
        ev = ((age >= lo) & (age < hi) & tum).sum()
        if pt > 0:
            out.append(dict(cond=label, band=lo, mid=(lo + min(hi, age.max())) / 2, pt=pt, ev=ev))
    return out


rows, curves = [], []
for series, (folder, files) in SERIES.items():
    data = {lvl: load(folder, pfx) for pfx, lvl in files.items()}
    levels = sorted(data)
    base = levels[0]
    for lvl in levels[1:]:
        t = pd.DataFrame(tally(data[base], 0) + tally(data[lvl], 1))
        t = t[t.groupby("band").ev.transform("sum") > 0]           # drop bands with no events in either
        t["dec"] = (t.mid - 60) / 10
        t["high"] = t.cond
        m = smf.glm("ev ~ C(band) + high + high:dec", t, family=sm.families.Poisson(), offset=np.log(t.pt)).fit()
        ci = m.conf_int()
        rows.append(dict(series=series, input0=base, input1=lvl,
                         HR_at_60=np.exp(m.params["high"]), HR_lo=np.exp(ci.loc["high", 0]), HR_hi=np.exp(ci.loc["high", 1]),
                         trend_per_decade=np.exp(m.params["high:dec"]),
                         trend_lo=np.exp(ci.loc["high:dec", 0]), trend_hi=np.exp(ci.loc["high:dec", 1])))
        # age-specific HR for plotting
        p = t.pivot(index=["band", "mid"], columns="cond", values=["ev", "pt"]).reset_index()
        for _, r in p.iterrows():
            e0, e1 = r[("ev", 0)], r[("ev", 1)]
            if e0 > 0 and e1 > 0:
                hr = (e1 / r[("pt", 1)]) / (e0 / r[("pt", 0)])
                se = np.sqrt(1 / e0 + 1 / e1)
                curves.append(dict(series=series, input0=base, input1=lvl, age=r[("mid", "")], HR=hr,
                                   lo=hr * np.exp(-1.96 * se), hi=hr * np.exp(1.96 * se), ev0=e0, ev1=e1))

res = pd.DataFrame(rows)
pd.DataFrame(curves).to_csv(OUT / "fig3_age_specific_HR.csv", index=False)
res.to_csv(OUT / "fig3_HR_trend_summary.csv", index=False)
pd.set_option("display.width", 200)
print(res.round(3).to_string(index=False))
