# FRB Likelihood Verification (Referee Comment M2)

## Context

The original manuscript equation for the marginal FRB probability density,
`P(DM_ext | z)`, contained a spurious external factor `<DM_cos>^-1` and a
rescaled argument of `P_cos`. This was corrected in the manuscript (current
`Eq. (6.7)` in `main_arxiv.tex`): `P_cos` is evaluated directly in physical
DM units, and the convolution with the host-galaxy log-normal PDF integrates
to unity without any additional normalization factor.

The response letter (`response_to_referee.tex`, item **M2**) states:

> "We confirm that this inconsistency was confined to the manuscript
> equation: the likelihood code already evaluated the truncated Gaussian
> directly in physical DM units, with no external factor and no rescaled
> argument, identically across all six models. This was verified directly
> against the FRB-likelihood implementation in each model's MCMC notebook,
> so the reported constraints are unaffected."

That claim was previously only qualitative ("we verified the code"), with no
numbers behind it. This folder numerically substantiates it.

**Nothing outside this folder was modified.** All scripts only *read* the
production notebooks and the real FRB catalog (`Data/FRB/frb_catalog.txt`);
they reimplement Eqs. (3.7) and (6.7) from scratch for cross-validation.

## Files

- `extract_and_diff_functions.py` — extracts the `P_DM_ext_vectorized`
  function (the FRB likelihood convolution) from the six main-analysis
  notebooks and checks they are identical up to whitespace/comments.
- `extracted_P_DM_ext_functions.txt` — the six extracted function bodies,
  for manual inspection.
- `verify_frb_likelihood.py` — the three numerical tests described below.
- `results_summary.json` — machine-readable results of the three tests.

## Method and results

### Step 1 — Code identity across the six models

`P_DM_ext_vectorized` (the truncated-Gaussian/log-normal convolution that
implements Eq. 6.7) was extracted from the six `*_CC_SNe_BAO_FRB.ipynb`
Baseline+FRB notebooks (LCDM, wCDM, CPL, AB, HS, ST — the ones that produced
Tables 2–7 and Table 10 of the manuscript) and compared after stripping
whitespace, comments, and line-wrapping.

**Result: identical in all six models.** Differences only appear afterward,
in each model's `lnlike_*` wrapper, where the model-specific `H(z)` is used
to compute the Macquart-relation mean `DM_IGM_mean` — exactly as expected,
since the FRB statistical treatment (Eqs. 3.7–3.9, 6.7) does not depend on
the underlying cosmological model.

### Step 2 — Normalization checks (adaptive quadrature)

Using the real FRB sample (104 bursts, `0.0085 ≤ z ≤ 1.354`) and the
ΛCDM Baseline+FRB best-fit point (`H0=70.46`, `Ωm=0.305`, `Ωb=0.0501`,
`σ_host=0.833`, `e^μ=104.2`, from Table 2 of the manuscript):

| Check | Redshifts tested | Max deviation from unity |
|---|---|---|
| `∫ P_cos(DM_cos\|z) dDM_cos = 1` (Eq. 3.7) | 20, spanning full sample range | **5.6 × 10⁻¹⁶** (machine precision) |
| `∫ P(DM_ext\|z) dDM_ext = 1` (Eq. 6.7) | 6, spanning full sample range | **3.0 × 10⁻⁴** (quadrature-limited, at the lowest z) |

Both normalizations hold to numerical precision across the full redshift
range of the sample, confirming the truncated Gaussian and its convolution
with the host PDF are correctly normalized as claimed in the text
surrounding Eq. (6.7).

### Step 3 — Production code vs. independent quadrature, on the real data

The production `P_DM_ext_vectorized` (fixed trapezoidal grid, `npts=300`,
the default used in every MCMC run) was compared point-by-point, for all
104 real FRBs at the same best-fit parameters, against an independent
from-scratch implementation of Eq. (6.7) using adaptive quadrature
(`scipy.integrate.quad`).

| | `npts=300` (production default) | `npts=20000` |
|---|---|---|
| Max relative difference | 6.75 % | 0.10 % |
| Median relative difference | 1.3 × 10⁻⁴ | 1.2 × 10⁻⁶ |

The large *maximum* relative difference at the production resolution is
**not a formula error** — it is a grid-resolution artifact confined to a
handful of the *lowest-redshift* FRBs in the sample, where the host-DM
integration grid (fixed at `npts=300` over `[0.01, max(DM_ext_obs)]`) is
coarse relative to the narrow region that actually contributes
(`DM_host ∈ [0, DM_ext_obs·(1+z)]`, which is small for nearby, low-DM_ext
bursts). Only **4 of 104 FRBs** (all at `z ≲ 0.08`) show a relative
difference above 1%; the worst offenders are:

| FRB | z | DM_ext,obs | rel. diff (npts=300) |
|---|---|---|---|
| FRB 20231230A | 0.0298 | 49.2 | 6.8 % |
| FRB 20181223C | 0.0302 | 41.7 | 3.1 % |
| FRB 20200723B | 0.0085 | 161.0 | 2.1 % |
| FRB 20231011A | 0.0783 | 77.7 | 1.3 % |

**Net impact on the reported likelihood is negligible.** Summed over all
104 FRBs:

- `ln L_FRB` (production, npts=300) = −649.728102
- `ln L_FRB` (independent quadrature) = −649.699917
- **Difference: −0.028** (compare to `Δln L` between models in Table 10,
  which range from ≈ 2 to ≈ 10 — two to three orders of magnitude larger).

This confirms that, while the fixed-grid production code has a known,
quantified, and small per-point numerical resolution limitation for a
handful of the closest FRBs, it does not materially affect `ln L_max`,
the reported cosmological constraints, or the model-comparison conclusions
of the paper.

## Bottom line

1. The FRB likelihood formula actually used in the MCMC code matches the
   corrected Eq. (6.7) — no spurious external factor, no rescaled argument
   — identically across all six models (confirmed by direct code
   comparison, not just visual inspection).
2. Both the truncated-Gaussian `P_cos` and its convolution `P(DM_ext|z)`
   are correctly normalized to unity, to within numerical-quadrature
   precision, across the full redshift range of the FRB sample.
3. The production code's fixed-grid trapezoidal integration introduces a
   small numerical resolution effect for a few of the lowest-redshift FRBs
   (max 6.8% per-point at npts=300), but this has a negligible effect
   (Δln L ≈ 0.03) on the total FRB log-likelihood and therefore on the
   reported constraints.

## How to reproduce

```bash
source /home/brunowesley/projetos/venv/bin/activate
python3 Codes/FRB-Likelihood-Verification/extract_and_diff_functions.py
python3 Codes/FRB-Likelihood-Verification/verify_frb_likelihood.py
```
