# MCMC-cosmo

Bayesian parameter estimation for cosmological models, using MCMC (`emcee`) to confront six expansion-history models against four observational probes — Cosmic Chronometers (CC), Type Ia Supernovae (Pantheon+SH0ES), DESI BAO DR2, and Fast Radio Burst (FRB) dispersion measures — individually and in combination, followed by a unified model comparison (AIC/BIC/LRT) and Figure-of-Merit analysis.

## Models

**GR-based** (`Codes/GR-based/`) — standard dark-energy parametrizations on top of GR:

| Model | Expansion history | Free parameters |
|-------|--------------------|------------------|
| `LCDM` | Flat ΛCDM: `E(z) = sqrt(Ωm0(1+z)³ + (1-Ωm0))` | `H0`, `Ωm0` |
| `wCDM` | Constant dark-energy EoS `w` | `H0`, `Ωm0`, `w` |
| `CPL`  | Chevallier–Polarski–Linder: `w(z) = w0 + wa·z/(1+z)` | `H0`, `Ωm0`, `w0`, `wa` |

**MG-based** (`Codes/MG-based/`) — `f(R)` modified-gravity models. `H(a)`/`yH(ln a)` is obtained by numerically solving the modified Friedmann ODE system (`scipy.integrate.solve_ivp`) rather than a closed-form `E(z)`, with a scalaron regularization scale `delta_s = 1e-7`:

| Model | `f(R)` form | Extra parameter |
|-------|-------------|------------------|
| `Appleby-Battye` (AB) | tanh-based broken `f(R)` | `b` (also fit with a Jeffreys prior on `log b` as an alternative run) |
| `Hu-Sawicki` (HS) | broken power law, `n = 1` | `mu_HS` (also saved reparametrized as `mu_HS/100` and `ln(mu_HS)` for prior/numerical-stability checks) |
| `Starobinsky` (ST) | `n = 1` | `lbd` (Λ-scale parameter) |

## Data (`Data/`)

- `CC/CC_Hz_data.txt` — cosmic chronometer `H(z)` measurements with errors.
- `SNe/Pantheon+SH0ES.dat` + `Pantheon+SH0ES_STAT+SYS.cov` — 1701 SNe Ia apparent magnitudes with full stat+sys covariance (a `z > 0.01` cut is applied before fitting).
- `BAO/desi_gaussian_bao_ALL_GCcomb_mean.txt` + `..._cov.txt` — DESI DR2 Gaussianized BAO distance measurements and covariance.
- `FRB/frb_catalog.txt` — FRB dispersion measures (`DM_obs`), Milky Way ISM + halo contributions, redshifts and host metadata, used to extract the extragalactic DM.

## Code layout (`Codes/`)

Each model directory (`LCDM`, `wCDM`, `CPL`, `Appleby-Battye`, `Hu-Sawicki`, `Starobinsky`) follows the same pattern:

- `CC/`, `SNe/`, `BAO/`, `FRB/` — single-probe MCMC fits.
- `CC_FRB/`, `BAO_FRB/`, `SNe_FRB/` — two-probe combinations with FRB.
- `CC_SNe_BAO/` — the three "standard" probes combined.
- `CC_SNe_BAO_FRB/` — full joint fit; also stores `log_prob_*.npy` alongside the chain.
- `Corner/` (or `Corners/`) — overlaid corner plots comparing dataset combinations for that model, plus `stats_*_full.txt` with parameter constraints at 68% and 95% CL.
- `FoM/` — Figure-of-Merit notebook: computes the inverse-area of the 2D highest-posterior-density contour (via Gaussian KDE) to quantify how much adding FRB data shrinks the `H0`–`Ωb` posterior.

Each fit notebook saves its flattened chain as `flat_samples_*.npy` (and `log_prob_*.npy` for the full joint fits) plus a GetDist corner plot (`Corner_*.png`).

`Codes/Model-comparison/Model_comparison.ipynb` builds a unified comparison table (`model_comparison.txt`) across all six models using the full `CC+SNe+BAO+FRB` fit: χ², AIC, BIC, ΔAIC/ΔBIC, an approximate Bayes factor (`lnBF ≈ -0.5·ΔBIC`), and a likelihood-ratio-test p-value, with verdicts against LCDM as the reference model.

## Requirements

```bash
pip install numpy pandas scipy matplotlib emcee getdist
```

MG-based notebooks additionally use `multiprocessing.Pool` for parallel walker evaluation and `scipy.interpolate.RegularGridInterpolator` to speed up repeated `H(z)` evaluations from precomputed grids (`H_array_*.npy`, `Hgrid_*.npz`).

## Usage

Each analysis is a self-contained Jupyter notebook. Open the relevant notebook (e.g. `Codes/GR-based/LCDM/CC_SNe_BAO/LCDM_CC_SNe_BAO.ipynb`) and run all cells — it loads the data from `Data/`, runs the `emcee` sampler, and writes the chain, log-probability, and corner plot into its own directory.

**Note**: data paths in the notebooks are currently absolute (e.g. `/home/brunowesley/projetos/MCMC-cosmo/Data/CC/CC_Hz_data.txt`); update them if running from a different location.

## Model comparison summary (CC + SNe + BAO + FRB)

| Model | k | Δ AIC vs LCDM | Δ BIC vs LCDM | LRT vs LCDM |
|-------|---|----------------|----------------|-------------|
| LCDM | 6 | reference | reference | reference |
| wCDM | 7 | −8.63 (strong for wCDM) | −3.17 (positive for wCDM) | p = 0.0011 |
| CPL  | 8 | −7.06 (strong for CPL) | +3.86 (positive for LCDM) | not significant |
| AB   | 7 | +2.26 (positive for LCDM) | +7.72 (strong for LCDM) | — |
| HS   | 7 | −8.17 (strong for HS) | −2.70 (positive for HS) | p = 0.0014 |
| ST   | 7 | −8.19 (strong for ST) | −2.72 (positive for ST) | p = 0.0014 |

Full numbers are in `Codes/Model-comparison/model_comparison.txt`.
