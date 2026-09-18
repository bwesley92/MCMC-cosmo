# MCMC Models

MCMC parameter estimation (`emcee`) for six expansion-history models against
four observational probes — Cosmic Chronometers (CC), Type Ia Supernovae
(Pantheon+SH0ES), DESI BAO DR2, and Fast Radio Burst (FRB) dispersion
measures — individually and in combination. This is the primary analysis;
`Codes/Model-comparison/`, `Codes/Robustness-test/`, and
`Codes/FRB-likelihood-verification/` all post-process the chains produced
here.

## Models

**`GR-based/`** — standard dark-energy parametrizations on top of General
Relativity:

| Model | Expansion history | Free parameters |
|-------|--------------------|------------------|
| `LCDM` | Flat ΛCDM: `E(z) = sqrt(Ωm0(1+z)³ + (1-Ωm0))` | `H0`, `Ωm0` |
| `wCDM` | Constant dark-energy EoS `w` | `H0`, `Ωm0`, `w` |
| `CPL`  | Chevallier–Polarski–Linder: `w(z) = w0 + wa·z/(1+z)` | `H0`, `Ωm0`, `w0`, `wa` |

**`MG-based/`** — `f(R)` modified-gravity models. `H(a)`/`yH(ln a)` is
obtained by numerically solving the modified Friedmann ODE system
(`scipy.integrate.solve_ivp`) rather than a closed-form `E(z)`, with a
scalaron regularization scale `delta_s = 1e-7`:

| Model | `f(R)` form | Extra parameter |
|-------|-------------|------------------|
| `Appleby-Battye` (AB) | tanh-based broken `f(R)` | `b` (also fit with a Jeffreys prior on `log b` as an alternative run) |
| `Hu-Sawicki` (HS) | broken power law, `n = 1` | `mu_HS` (also saved reparametrized as `mu_HS/100` and `ln(mu_HS)` for prior/numerical-stability checks) |
| `Starobinsky` (ST) | `n = 1` | `lbd` (Λ-scale parameter) |

All six models additionally fit `Ωb`, `MB` (SNe absolute magnitude), and,
when FRB data is included, `sigma_host` and `exp_mu` (host-galaxy DM
log-normal nuisance parameters).

## Layout

Each model directory (`LCDM`, `wCDM`, `CPL`, `Appleby-Battye`, `Hu-Sawicki`,
`Starobinsky`) follows the same pattern:

- `CC/`, `SNe/`, `BAO/`, `FRB/` — single-probe MCMC fits.
- `CC_FRB/`, `BAO_FRB/`, `SNe_FRB/` — two-probe combinations with FRB.
- `CC_SNe_BAO/` — the three "standard" probes combined (**Baseline**).
- `CC_SNe_BAO_FRB/` — full joint fit (**Baseline + FRB**); also stores
  `log_prob_*.npy` alongside the chain.
- `Corner/` (or `Corners/`) — notebook that overlays the corner plots for
  every dataset combination of that model, plus `stats_*_full.txt` with
  parameter constraints at 68% and 95% CL.
- `FoM/` — Figure-of-Merit notebook: computes the inverse-area of the 2D
  highest-posterior-density contour (via Gaussian KDE) to quantify how much
  adding FRB data shrinks the `H0`–`Ωb` posterior.
- `ObH0fd/` — reports the marginalized `Ωb·H0·f_d` combination (the
  parameter combination FRB dispersion measures actually constrain,
  via the Macquart relation) for the Baseline+FRB fit.

Each fit notebook is self-contained: it loads the data from `Data/`, runs
the `emcee` sampler, and saves its flattened chain as `flat_samples_*.npy`
(plus `log_prob_*.npy` for the full joint fits) and a GetDist corner plot
(`Corner_*.png`) into its own directory.

### MG-based grids

`MG-based/` notebooks that use the full joint or Baseline dataset combine
`multiprocessing.Pool` for parallel walker evaluation with
`scipy.interpolate.RegularGridInterpolator`, evaluating `H(z)` from a
precomputed 4D grid (`H0`, `Ωm0`, extra parameter, `z`) instead of solving
the ODE system at every likelihood call. Most of these notebooks contain
**two ways to run**, both present in the same cell block:

1. A pre-built grid already exists (`H_array_*.npy`) — the default: just
   run the "Load grid" cell(s).
2. No saved grid, or a different resolution/parameter range is needed —
   uncomment and run the grid-building block first (parallel construction
   + interpolation-error validation against direct ODE solutions), then
   load it.

## Requirements

```bash
pip install numpy pandas scipy matplotlib emcee getdist
```

## Usage

Open the relevant notebook (e.g.
`GR-based/LCDM/CC_SNe_BAO/LCDM_CC_SNe_BAO.ipynb`) and run all cells.

**Note**: data paths in the notebooks are currently absolute (e.g.
`/home/brunowesley/projetos/MCMC-cosmo/Data/CC/CC_Hz_data.txt`); update
them if running from a different location.
