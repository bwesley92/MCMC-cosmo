# Robustness Tests

Two additional Baseline+FRB (`CC_SNe_BAO_FRB`) MCMC runs per model, each
relaxing one assumption that is held fixed in the main analysis
(`Codes/MCMC-Models/`), to check whether the reported constraints depend on
it. Run for three representative models spanning both families — `LCDM`,
`CPL` (GR-based) and `Hu-Sawicki` n=1 (MG-based) — and compared against
each model's fiducial Baseline+FRB chain.

## Tests

### `DMbudget` — Galactic DM budget treatment

The main analysis fixes the Galactic contribution to the FRB dispersion
measure: `DM_MW,ISM` to the NE2001 point estimate, and `DM_MW,halo` to a
single fiducial value. This test instead propagates the uncertainty on
both:

- `DM_MW,ISM`: a systematic uncertainty `sigma_MW = 30 pc/cm^3` (ATNF
  Pulsar Catalogue calibration, Manchester et al. 2005) is added in
  quadrature to the cosmological scatter.
- `DM_MW,halo`: marginalized via a Gaussian convolution,
  `DM_MW,halo ~ N(mu=50, sigma)` truncated to the plausible range
  `[25, 110] pc/cm^3` (`sigma` = range/√12), instead of being held fixed.

This follows the treatment of Kalita et al. (2024).

### `DMbudget_CCsubset` — `DMbudget` + correlated CC removal

Same Galactic-DM treatment as `DMbudget`, plus: 2 of the 35 Cosmic
Chronometer measurements used in the main analysis (Jimenez et al. 2023 and
Jiao et al. 2023) are recalibrations of the same underlying LEGA-C sample
as a third (Borghi et al. 2022, kept). This test removes the two derived
measurements, leaving 33 independent CC points, to check sensitivity to
this hidden correlation (per Moresco 2024).

## Layout

Per model (`LCDM/`, `CPL/`, `Hu-Sawicki/`):

- `CC_SNe_BAO_FRB/` — `<Model>_CC_SNe_BAO_FRB_DMbudget.ipynb` and
  `..._DMbudget_CCsubset.ipynb`, the two MCMC runs; each saves its
  `flat_samples_*.npy` and `log_prob_*.npy` here.
- `Corner/` — `<model>_corner_DMbudget.ipynb` and
  `..._DMbudget_CCsubset.ipynb` load the corresponding robustness chain
  together with the model's fiducial Baseline+FRB chain from
  `../../../MCMC-Models/.../CC_SNe_BAO_FRB/`, overlay them in a single
  corner plot (`Corner_*.pdf`/`.png`), and save marginalized constraints
  to `stats_*.txt`.

## How to reproduce

1. The model's fiducial Baseline+FRB chain must already exist under
   `Codes/MCMC-Models/` (both corner notebooks load it for the overlay).
2. Run `<Model>_CC_SNe_BAO_FRB_DMbudget.ipynb` and/or
   `..._DMbudget_CCsubset.ipynb`.
3. Run the matching `<model>_corner_....ipynb` in `Corner/`.

## References

1. Kalita, S., Bhatporia, S., & Weltman, A. (2025), *Fast radio bursts as
   probes of the late-time universe: A new insight on the Hubble tension*,
   [arXiv:2410.01974](https://arxiv.org/abs/2410.01974) — the Galactic-DM
   marginalization treatment adopted in `DMbudget`.
2. Manchester, R. N., Hobbs, G. B., Teoh, A., & Hobbs, M. (2005), *The
   Australia Telescope National Facility Pulsar Catalogue*,
   [arXiv:astro-ph/0412641](https://arxiv.org/abs/astro-ph/0412641) — source
   of the `sigma_MW` calibration.
3. Moresco, M. (2024), *Measuring the expansion history of the Universe
   with Cosmic Chronometers*,
   [arXiv:2412.01994](https://arxiv.org/abs/2412.01994) — identifies the
   LEGA-C correlation addressed in `DMbudget_CCsubset`.
