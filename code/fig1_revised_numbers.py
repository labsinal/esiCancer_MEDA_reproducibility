"""
Numbers behind the revised Figure 1 (same layout as the submitted figure: A-E).
Writes output/fig1_revised_numbers.json and prints a readable summary.

A  NMSC age-standardised incidence vs UV (unchanged analysis; values from the executed NMSC notebook)
B  NMSC median age of diagnosed cases vs UV (corrected grouped-median formula)
C  Lung cancer death, current vs never smokers (Halpern 1993 Table 4, verified):
     bars  = share of deaths in each age band if both groups had the same person-time (pooled)
     inset = age-specific death rate ratio (current / never) with 95% CI
D  Breast cancer, BCAC 2021 (population-based studies, Table S12):
     bars  = share of cases in each age group: BRCA1 and BRCA2 PTV carriers (published counts only;
             no non-carrier distribution is derived)
     inset = age-specific OR+ with 95% CI (published)
     age shift = case-only OR per year of age (published, Table S11)
E  Summary table: risk increase, MEDA in each group, MEDA difference, change in relative risk per decade
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from paths import OUT, EPI
RNG = np.random.default_rng(1)
B = 10000
res = {}


def gmedian(lo, hi, w):
    c = np.cumsum(w) / np.sum(w); i = int(np.searchsorted(c, 0.5)); p = c[i - 1] if i else 0.0
    return lo[i] + (0.5 - p) / (c[i] - p) * (hi[i] - lo[i])


# ---------------- A, B: NMSC ----------------
st = pd.read_csv(OUT / "nmsc" / "state_level_results.csv")
for key, col in (("A", "asr"), ("B", "meda_cases")):
    r = stats.linregress(st.uv, st[col]); t = stats.t.ppf(0.975, len(st) - 2)
    res[key] = dict(n=len(st), R2=r.rvalue ** 2, p=r.pvalue, slope_per_1000=1000 * r.slope,
                    slope_lo=1000 * (r.slope - t * r.stderr), slope_hi=1000 * (r.slope + t * r.stderr),
                    min=st[col].min(), max=st[col].max())
hi_uv, lo_uv = st[st.uv > 5000], st[st.uv < 4000]
bd, br = [], []
for _ in range(B):
    h = hi_uv.sample(len(hi_uv), replace=True, random_state=RNG.integers(1e9))
    l = lo_uv.sample(len(lo_uv), replace=True, random_state=RNG.integers(1e9))
    br.append(h.asr.mean() / l.asr.mean()); bd.append(h.meda_cases.mean() - l.meda_cases.mean())
res["E_NMSC"] = dict(n_high=len(hi_uv), n_low=len(lo_uv), RR=hi_uv.asr.mean() / lo_uv.asr.mean(),
                     RR_lo=np.percentile(br, 2.5), RR_hi=np.percentile(br, 97.5),
                     meda_low=lo_uv.meda_cases.mean(), meda_high=hi_uv.meda_cases.mean(),
                     dmeda=hi_uv.meda_cases.mean() - lo_uv.meda_cases.mean(),
                     dmeda_lo=np.percentile(bd, 2.5), dmeda_hi=np.percentile(bd, 97.5))
tr = pd.read_csv(OUT / "fig1_RR_trend_summary.csv")
res["E_NMSC"]["trend"] = float(tr.loc[tr.example.str.startswith("NMSC"), "trend_per_decade"].iloc[0])

# ---------------- C: lung ----------------
L = pd.read_csv(EPI / "lung" / "lung_halpern_table4.csv")
nv, cu = L[L.Status == "never"].reset_index(drop=True), L[L.Status == "current"].reset_index(drop=True)
lo = nv.age_range.str.split("-").str[0].astype(float).values
hi = nv.age_range.str.split("-").str[1].astype(float).values + 1
pooled = (nv.n + cu.n).values
w_n, w_c = nv.deaths.values / nv.n.values * pooled, cu.deaths.values / cu.n.values * pooled
std_rate = lambda d, py: 1e5 * (d / py * pooled).sum() / pooled.sum()     # directly age-standardised, pooled PY
med_n, med_c = gmedian(lo, hi, w_n), gmedian(lo, hi, w_c)
bm = []
for _ in range(B):
    dn, dc = RNG.poisson(nv.deaths.values), RNG.poisson(cu.deaths.values)
    bm.append(gmedian(lo, hi, dc / cu.n.values * pooled) - gmedian(lo, hi, dn / nv.n.values * pooled))
rr_age = []
for i in range(len(nv)):
    dn, dc = nv.deaths[i], cu.deaths[i]
    if dn and dc:
        rr = (dc / cu.n[i]) / (dn / nv.n[i]); se = np.sqrt(1 / dn + 1 / dc)
        rr_age.append([nv.age_range[i], rr, rr * np.exp(-1.96 * se), rr * np.exp(1.96 * se)])
summ = pd.read_csv(OUT / "fig1_RR_trend_summary.csv")
lung_tr = summ[summ.example.str.startswith("Lung")].iloc[0]
res["C"] = dict(bands=list(nv.age_range), share_never=list(w_n / w_n.sum()), share_current=list(w_c / w_c.sum()),
                std_rate_never=std_rate(nv.deaths.values, nv.n.values), std_rate_current=std_rate(cu.deaths.values, cu.n.values),
                rr_by_age=rr_age, median_never=med_n, median_current=med_c)
res["E_lung"] = dict(RR=21.14, RR_lo=18.98, RR_hi=23.56,       # Poisson, offset log PY (lung_halpern_results.txt)
                     meda_low=med_n, meda_high=med_c, dmeda=med_c - med_n,
                     dmeda_lo=np.percentile(bm, 2.5), dmeda_hi=np.percentile(bm, 97.5),
                     trend=lung_tr.trend_per_decade, trend_lo=lung_tr.trend_lo, trend_hi=lung_tr.trend_hi,
                     interaction_p=0.796)

# ---------------- D: BRCA ----------------
# published values: BCAC 2021 Supplementary Tables S12 (population-based studies) and S11 (case-only)
s12 = pd.read_csv(EPI / "brca" / "bcac2021_tableS12_population_based.csv")
s11 = pd.read_csv(EPI / "brca" / "bcac2021_tableS11_case_only.csv").set_index("gene")
b1 = s12[s12.gene == "BRCA1"]
blo, bhi = b1.age_lo.to_numpy(float), b1.age_hi.to_numpy(float)   # open groups closed at 20 and 90 y
counts = {g: t.carrier_cases.to_numpy() for g, t in s12.groupby("gene")}
orp = {g: (list(t.OR_vs_all_controls), list(zip(t.OR_lo, t.OR_hi))) for g, t in s12.groupby("gene")}
res["D"] = dict(groups=["<40", "40-49", "50-59", "60+"], share={g: list(c / c.sum()) for g, c in counts.items()},
                orplus=orp)
case_only = {g: tuple(r) for g, r in s11[["case_only_OR_per_year_of_age", "lo", "hi"]].iterrows()}
for g, c in counts.items():
    med = gmedian(blo, bhi, c)
    boot = [gmedian(blo, bhi, RNG.multinomial(c.sum(), c / c.sum())) for _ in range(B)]
    p, l, h = case_only[g]
    res[f"E_{g}"] = dict(OR_young=orp[g][0][0], OR_old=orp[g][0][-1], meda_high=med,
                         meda_high_lo=np.percentile(boot, 2.5), meda_high_hi=np.percentile(boot, 97.5),
                         case_only_or=p, case_only_lo=l, case_only_hi=h,
                         trend=p ** 10, trend_lo=l ** 10, trend_hi=h ** 10)

OUT.mkdir(exist_ok=True)
(OUT / "fig1_revised_numbers.json").write_text(json.dumps(res, indent=2, default=float))
for k, v in res.items():
    print(k, {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in v.items() if not isinstance(vv, (list, dict))})
