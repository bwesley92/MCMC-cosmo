# MCMC-cosmo

Bayesian parameter estimation for six expansion-history models, using MCMC
(`emcee`) to confront them against four observational probes — Cosmic
Chronometers (CC), Type Ia Supernovae (Pantheon+SH0ES), DESI BAO DR2, and
Fast Radio Burst (FRB) dispersion measures — individually and in every
combination. Each fit is followed by a Figure-of-Merit analysis quantifying
how much the FRB data tightens the constraints, and a report of the
`Ωb·H0·f_d` combination that FRB dispersion measures directly probe. The
resulting posterior chains are then compared across models (information
criteria and direct Bayesian evidence), stress-tested against specific
modeling assumptions, and the FRB likelihood implementation is
independently cross-checked.

## Models

**GR-based** (`Codes/MCMC-Models/GR-based/`) — standard dark-energy
parametrizations on top of GR:

| Model | Expansion history | Free parameters |
|-------|--------------------|------------------|
| `LCDM` | Flat ΛCDM: `E(z) = sqrt(Ωm0(1+z)³ + (1-Ωm0))` | `H0`, `Ωm0` |
| `wCDM` | Constant dark-energy EoS `w` | `H0`, `Ωm0`, `w` |
| `CPL`  | Chevallier–Polarski–Linder: `w(z) = w0 + wa·z/(1+z)` | `H0`, `Ωm0`, `w0`, `wa` |

**MG-based** (`Codes/MCMC-Models/MG-based/`) — `f(R)` modified-gravity
models. `H(a)`/`yH(ln a)` is obtained by numerically solving the modified
Friedmann ODE system (`scipy.integrate.solve_ivp`) rather than a
closed-form `E(z)`, with a scalaron regularization scale `delta_s = 1e-7`:

| Model | `f(R)` form | Extra parameter |
|-------|-------------|------------------|
| `Appleby-Battye` (AB) | tanh-based broken `f(R)` | `b` (also fit with a Jeffreys prior on `log b` as an alternative run) |
| `Hu-Sawicki` (HS) | broken power law, `n = 1` | `mu_HS` (also saved reparametrized as `mu_HS/100` and `ln(mu_HS)` for prior/numerical-stability checks) |
| `Starobinsky` (ST) | `n = 1` | `lbd` (Λ-scale parameter) |

All six models additionally fit `Ωb`, `MB` (SNe absolute magnitude), and,
whenever FRB data is included, `sigma_host` and `exp_mu` (host-galaxy DM
log-normal nuisance parameters).

## Data (`Data/`)

All four datasets below are third-party observational data, redistributed
only for reproducibility; see `Data/README.md` for full citations and
provenance of each.

- `CC/CC_Hz_data.txt` — cosmic chronometer `H(z)` measurements with errors.
- `SNe/Pantheon+SH0ES.dat` + `Pantheon+SH0ES_STAT+SYS.cov` — 1701 SNe Ia
  apparent magnitudes with full stat+sys covariance (a `z > 0.01` cut is
  applied before fitting).
- `BAO/desi_gaussian_bao_ALL_GCcomb_mean.txt` + `..._cov.txt` — DESI DR2
  Gaussianized BAO distance measurements and covariance.
- `FRB/frb_catalog.txt` — 104 localized FRBs: observed dispersion measure
  (`DM_obs`), Milky Way ISM contribution (NE2001 model), redshift, host
  type, and reference. The Galactic halo and host-galaxy DM contributions
  are treated as model parameters rather than read from the catalog (see
  `Codes/MCMC-Models/README.md`).

## Analysis pipeline (`Codes/`)

- **`MCMC-Models/`** — the primary analysis: one `emcee` fit per model per
  dataset combination, plus, per model, a Figure-of-Merit notebook (area of
  the 2D `H0`–`Ωb` credible region, with vs. without FRB data) and a report
  of the marginalized `Ωb·H0·f_d` combination for the Baseline+FRB fit. See
  `Codes/MCMC-Models/README.md` for the full per-model layout, output
  conventions, and the interpolation-grid strategy used by the `f(R)`
  models.
- **`Model-comparison/`** — compares the six models' full-joint-fit
  posteriors two ways: `AIC-BIC-LRT/` computes the information criteria
  and a likelihood-ratio test; `Bayesian-Evidence/` computes a genuine
  k-nearest-neighbour estimate of the Bayesian evidence (Heavens et al.
  2017) directly from the chains, for a Bayes-factor cross-check of the
  BIC-based approximation.
- **`Robustness-test/`** — re-runs the Baseline+FRB fit for three
  representative models under two relaxed assumptions (the fixed Galactic
  dispersion-measure budget, and a possible hidden correlation between
  Cosmic Chronometer measurements), to check whether the reported
  constraints depend on them. See `Codes/Robustness-test/README.md`.
- **`FRB-likelihood-verification/`** — an independent, from-scratch
  reimplementation of the FRB likelihood, used to numerically confirm that
  the production code is implemented identically across all six models
  and correctly normalized. See `Codes/FRB-likelihood-verification/README.md`.

## Requirements

```bash
pip install numpy pandas scipy matplotlib emcee getdist scikit-learn
```

(`scikit-learn` is only needed for the Bayesian-evidence k-NN estimator in
`Codes/Model-comparison/Bayesian-Evidence/`.)

## Usage

Each analysis is a self-contained Jupyter notebook. Open the relevant
notebook (e.g. `Codes/MCMC-Models/GR-based/LCDM/CC_SNe_BAO/LCDM_CC_SNe_BAO.ipynb`)
and run all cells — it loads the data from `Data/`, runs the `emcee`
sampler, and writes the chain, log-probability, and corner plot into its
own directory.

**Note**: data paths in the notebooks are currently absolute (e.g.
`/home/brunowesley/projetos/MCMC-cosmo/Data/CC/CC_Hz_data.txt`); update
them if running from a different location.
