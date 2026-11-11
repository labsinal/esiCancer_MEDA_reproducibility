"""
Reanalysis of Figure 1C (lung cancer mortality, current vs never smokers) for Reviewer
comments 1, 4 and 6.

Source: Halpern, Gillespie & Warner 1993, JNCI 85:457-464, Appendix Table 4 (ACS CPS-II),
transcribed in lung_halpern_table4.csv (checked cell by cell against the printed table; the older
lung_halpern.csv had 248279 for 248278 and 611 for 622). `n` = person-years, `deaths` = lung cancer
deaths (ICD-9 trachea, bronchus, lung), sexes pooled in the source table.

Outputs
  1. Crude and age-adjusted rate ratios (Poisson GLM with log person-year offset; Mantel-Haenszel).
  2. Age-specific rate ratios with 95% CI.
  3. LRT for smoking x age interaction (a test of proportional rates across age, NOT of MEDA equality).
  4. Median age at death three ways, so the choice is explicit:
       a. as in the original notebook: median of the normalised rate profile (each band weighted 1);
       b. observed deaths in the cohort (depends on each group's person-time age structure);
       c. age-standardised: both groups' rates applied to the same person-time distribution
          (pooled person-years), i.e. the distribution of ages at death the two groups would
          show with identical age structure.
     Medians are interpolated linearly within bands; CIs by parametric bootstrap (Poisson deaths).
"""
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

from paths import EPI
d = pd.read_csv(EPI / "lung" / "lung_halpern_table4.csv")   # verified against printed Table 4 (2026-10-15)
d["lo"] = d.age_range.str.split("-").str[0].astype(float)
d["hi"] = d.age_range.str.split("-").str[1].astype(float) + 1     # band 40-43 covers [40, 44)
d["rate"] = d.deaths / d.n
RNG = np.random.default_rng(1993)
B = 5000


def interp_median(lo, hi, w):
    w = np.asarray(w, float)
    c = np.cumsum(w) / w.sum()
    i = int(np.searchsorted(c, 0.5))
    prev = c[i - 1] if i else 0.0
    return lo[i] + (0.5 - prev) / (c[i] - prev) * (hi[i] - lo[i])


def medians(never_deaths, curr_deaths, never_py, curr_py, lo, hi):
    rn, rc = never_deaths / never_py, curr_deaths / curr_py
    pooled = never_py + curr_py
    return {
        "a_normalised_rate": (interp_median(lo, hi, rn), interp_median(lo, hi, rc)),
        "b_observed_deaths": (interp_median(lo, hi, never_deaths), interp_median(lo, hi, curr_deaths)),
        "c_age_standardised": (interp_median(lo, hi, rn * pooled), interp_median(lo, hi, rc * pooled)),
    }


nv, cu = d[d.Status == "never"].reset_index(drop=True), d[d.Status == "current"].reset_index(drop=True)
lo, hi = nv.lo.values, nv.hi.values

print("=== Totals ===")
for name, g in (("never", nv), ("current", cu)):
    print(f"{name:8s} deaths={g.deaths.sum():5d}  person-years={g.n.sum():,}  rate/100k PY={1e5 * g.deaths.sum() / g.n.sum():.2f}")
crude = (cu.deaths.sum() / cu.n.sum()) / (nv.deaths.sum() / nv.n.sum())
print(f"Crude rate ratio current/never = {crude:.3f}   (manuscript Fig 1E rInc = 16.164)")

# Poisson models with person-time offset
d["current"] = (d.Status == "current").astype(int)
d["age_band"] = pd.Categorical(d.age_range)
m0 = smf.glm("deaths ~ current + age_band", d, family=sm.families.Poisson(), offset=np.log(d.n)).fit()
m1 = smf.glm("deaths ~ current * age_band", d, family=sm.families.Poisson(), offset=np.log(d.n)).fit()
lr = 2 * (m1.llf - m0.llf); dfd = m1.df_model - m0.df_model
print(f"\nAge-adjusted RR (Poisson, common RR across age) = {np.exp(m0.params['current']):.2f} "
      f"[{np.exp(m0.conf_int().loc['current', 0]):.2f}, {np.exp(m0.conf_int().loc['current', 1]):.2f}]")
print(f"LRT smoking x age interaction: chi2 = {lr:.2f}, df = {dfd:.0f}, p = {stats.chi2.sf(lr, dfd):.3g}")
print(f"Main-effects model deviance = {m0.deviance:.2f} on {m0.df_resid:.0f} df (goodness of fit p = {stats.chi2.sf(m0.deviance, m0.df_resid):.3g})")
# Mantel-Haenszel rate ratio
w = nv.n * cu.n / (nv.n + cu.n)
mh = (cu.deaths * nv.n / (nv.n + cu.n)).sum() / (nv.deaths * cu.n / (nv.n + cu.n)).sum()
print(f"Mantel-Haenszel RR = {mh:.2f}")

print("\n=== Age-specific rates and rate ratios (exact Poisson CI for ratio via conditional binomial) ===")
rows = []
for i in range(len(nv)):
    dn, dc, pn, pc = nv.deaths[i], cu.deaths[i], nv.n[i], cu.n[i]
    rr = (dc / pc) / (dn / pn) if dn else np.inf
    # conditional binomial: dc | dn+dc ~ Bin(n, pc*RR/(pn + pc*RR))
    t = dn + dc
    if dn and dc:
        plo, phi = stats.beta.ppf(0.025, dc, t - dc + 1), stats.beta.ppf(0.975, dc + 1, t - dc)
        rlo, rhi = plo / (1 - plo) * pn / pc, phi / (1 - phi) * pn / pc
    else:
        rlo, rhi = np.nan, np.nan
    rows.append(dict(age=nv.age_range[i], never_rate_100k=1e5 * dn / pn, current_rate_100k=1e5 * dc / pc,
                     RR=rr, RR_lo=rlo, RR_hi=rhi))
print(pd.DataFrame(rows).round(2).to_string(index=False))

print("\n=== Median age at lung cancer death (years), with parametric-bootstrap 95% CI ===")
point = medians(nv.deaths.values, cu.deaths.values, nv.n.values, cu.n.values, lo, hi)
boots = {k: [] for k in point}
for _ in range(B):
    bn, bc = RNG.poisson(nv.deaths.values), RNG.poisson(cu.deaths.values)
    for k, (mn, mc) in medians(bn, bc, nv.n.values, cu.n.values, lo, hi).items():
        boots[k].append((mn, mc, mc - mn, mc / mn))
for k, (mn, mc) in point.items():
    b = np.array(boots[k])
    q = lambda j: np.percentile(b[:, j], [2.5, 97.5]).round(2)
    print(f"{k:20s} never={mn:.2f} {q(0)}  current={mc:.2f} {q(1)}  diff={mc - mn:+.2f} {q(2)}  rMEDA={mc / mn:.3f} {q(3).round(3)}")
print("\nNote: the open top band (74-80) is closed at 81; the cohort has no follow-up beyond 80,")
print("so every median here is the median among deaths occurring at ages 40-80 within 6 years of follow-up.")
