"""
One common measure for every example in Figure 1 (and for the simulations in Figure 3):

  - age-specific relative risk (RR) of the higher-risk vs lower-risk group, by age band;
  - its overall level (age-adjusted RR): how much the risk factor raises incidence;
  - its trend with age, as the multiplicative change in RR per decade of age (exp of the slope of
    log RR on age). 1.0 = the same multiplier at all ages (proportional hazards): more cases, same
    age distribution, MEDA unchanged. < 1 = extra risk concentrated at young ages: MEDA earlier.

Data
  Lung   : Halpern 1993 Table 4 (lung_halpern_table4.csv), deaths and person-years.
           Poisson GLM, deaths ~ smoking * age(continuous, decades), offset log(person-years).
  BRCA   : BCAC 2021 appendix. Table S11 (case-only OR per year, population-based studies) is the
           primary estimate. Check: weighted regression of log OR+ (Table S12) on age-band midpoints.
  NMSC   : GBD 2023 age-specific Number and Rate, states pooled into high-UV (>5000 Wh/m2) and
           low-UV (<4000) groups; RR = pooled rate ratio by age band, ages 30+. GBD values are
           modelled estimates, not observed counts, so no sampling CI is given (descriptive only).
"""
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

from paths import OUT, EPI
rows = []          # tidy age-specific RRs for plotting
summary = []       # one line per example

# ---------------- Lung ----------------
d = pd.read_csv(EPI / "lung" / "lung_halpern_table4.csv")
lo = d.age_range.str.split("-").str[0].astype(float)
hi = d.age_range.str.split("-").str[1].astype(float) + 1
d["mid"] = (lo + hi) / 2
d["current"] = (d.Status == "current").astype(int)
d["dec"] = (d.mid - 60) / 10
m = smf.glm("deaths ~ current * dec", d, family=sm.families.Poisson(), offset=np.log(d.n)).fit()
ci = m.conf_int()
summary.append(dict(example="Lung cancer death, current vs never smokers", source="Halpern 1993, Table 4",
                    RR_at_60=np.exp(m.params["current"]), RR_lo=np.exp(ci.loc["current", 0]), RR_hi=np.exp(ci.loc["current", 1]),
                    trend_per_decade=np.exp(m.params["current:dec"]),
                    trend_lo=np.exp(ci.loc["current:dec", 0]), trend_hi=np.exp(ci.loc["current:dec", 1])))
nv, cu = d[d.Status == "never"].reset_index(), d[d.Status == "current"].reset_index()
for i in range(len(nv)):
    if nv.deaths[i] == 0:
        continue
    rr = (cu.deaths[i] / cu.n[i]) / (nv.deaths[i] / nv.n[i])
    se = np.sqrt(1 / cu.deaths[i] + 1 / nv.deaths[i])
    rows.append(dict(example="lung", age=nv.mid[i], RR=rr, lo=rr * np.exp(-1.96 * se), hi=rr * np.exp(1.96 * se)))

# ---------------- BRCA ----------------
# Table S12 OR+ (cases in age group vs all controls) with 95% CI; midpoints assume <40 ~ 35, 60+ ~ 67
mids = np.array([35, 45, 55, 67.])
S12 = {"BRCA1": ([46.3, 14.2, 6.52, 2.33], [(33.4, 64.1), (10.4, 19.5), (4.66, 9.14), (1.57, 3.48)]),
       "BRCA2": ([18.7, 7.85, 5.14, 2.81], [(14.4, 24.1), (6.23, 9.88), (4.09, 6.47), (2.20, 3.60)])}
S11 = {"BRCA1": (0.941, 0.94, 0.95), "BRCA2": (0.968, 0.96, 0.97)}   # case-only OR per year
for g, (orp, cis) in S12.items():
    y = np.log(orp)
    se = np.array([(np.log(h) - np.log(l)) / 3.92 for l, h in cis])
    w = sm.WLS(y, sm.add_constant((mids - 60) / 10), weights=1 / se**2).fit()
    slope, s_se = w.params[1], w.bse[1]
    for a, o, (l, h) in zip(mids, orp, cis):
        rows.append(dict(example=g, age=a, RR=o, lo=l, hi=h))
    p, l, h = S11[g]
    summary.append(dict(example=f"Breast cancer, {g} PTV carriers vs non-carriers", source="BCAC 2021, Tables S11/S12",
                        RR_at_60=np.exp(w.params[0]), RR_lo=np.nan, RR_hi=np.nan,
                        trend_per_decade=p**10, trend_lo=l**10, trend_hi=h**10,
                        check_trend_from_S12=np.exp(slope),
                        check_lo=np.exp(slope - 1.96 * s_se), check_hi=np.exp(slope + 1.96 * s_se)))

# ---------------- NMSC ----------------
raw = pd.read_csv(EPI / "nmsc" / "IHME_data_USA_states.csv")
raw = raw[(raw.measure == "Incidence") & raw.metric.isin(["Number", "Rate"]) & (raw.year == 2023)]
w_ = raw.pivot_table(index=["location", "sex", "age"], columns="metric", values="val").reset_index()
w_ = w_[w_.Rate > 0]
w_["pop"] = w_.Number / w_.Rate * 1e5
uv = pd.read_csv(OUT / "nmsc" / "state_level_results.csv")[["location", "uv"]]
w_ = w_.merge(uv, on="location")
w_["grp"] = np.where(w_.uv > 5000, "high", np.where(w_.uv < 4000, "low", None))
w_ = w_[w_.grp.notna()]


def age_bounds(s):
    s = s.replace(" years", "").replace("+", "-100")
    a, b = s.split("-")
    return float(a), float(b) + 1


w_["lo"] = w_.age.map(lambda s: age_bounds(s)[0]); w_["hi"] = w_.age.map(lambda s: age_bounds(s)[1])
w_ = w_[w_.lo >= 30]
g = w_.groupby(["grp", "lo", "hi"], as_index=False)[["Number", "pop"]].sum()
g["rate"] = g.Number / g["pop"]
t = g.pivot(index=["lo", "hi"], columns="grp", values="rate").reset_index()
t["RR"] = t.high / t.low
t["mid"] = (t.lo + t.hi) / 2
for _, r in t.iterrows():
    rows.append(dict(example="nmsc", age=r.mid, RR=r.RR, lo=np.nan, hi=np.nan))
cases = g.pivot(index=["lo", "hi"], columns="grp", values="Number").sum(axis=1).values
fit = sm.WLS(np.log(t.RR.values), sm.add_constant((t.mid.values - 60) / 10), weights=cases).fit()
summary.append(dict(example="NMSC incidence, high- vs low-UV states (ecological)", source="GBD 2023 + NCI UV",
                    RR_at_60=np.exp(fit.params[0]), RR_lo=np.nan, RR_hi=np.nan,
                    trend_per_decade=np.exp(fit.params[1]), trend_lo=np.nan, trend_hi=np.nan))

rr = pd.DataFrame(rows); sm_ = pd.DataFrame(summary)
rr.to_csv(OUT / "fig1_age_specific_RR.csv", index=False)
sm_.to_csv(OUT / "fig1_RR_trend_summary.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 20)
print(sm_.round(3).to_string(index=False))
print("\nNMSC age-specific RR (high/low UV):"); print(t[["lo", "hi", "RR"]].round(3).to_string(index=False))
