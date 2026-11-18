# LRT calibration by simulations under the ΛCDM null

`lrt_bootstrap.py` calibrates the likelihood-ratio test (LRT) used in the
model comparison (`../AIC-BIC-LRT/`) with a parametric bootstrap, instead of
relying on Wilks' theorem.

Why: Wilks' theorem (Λ_LR ~ χ²_ΔK under the null) requires the null model to
be nested in the alternative at an *interior* point of its parameter space,
with the regularity conditions satisfied. This holds for wCDM (w = −1) and
CPL ((w0, wa) = (−1, 0)), but not for

- **Hu–Sawicki**: ΛCDM is exactly recovered at μ_HS = 0, the lower *boundary*
  of the parameter space (and H(z) depends on μ_HS only through μ_HS²);
- **Starobinsky**: ΛCDM is only approached for λ_S → ∞, which is *outside* the
  parameter range λ_S ∈ [1.54, 10] used in the fits (the models are not
  nested within the fitted parameter space);
- **Appleby–Battye**: same as ST (b → ∞, outside b ∈ [1.6, 12]); its observed
  statistic is negative, so it is not calibrated here.

## Method

1. `observed`: maximum-likelihood fits of ΛCDM, wCDM, CPL, HS and ST to the
   real CC + SNe Ia + BAO + FRB data (bounded, multi-start Powell /
   Nelder–Mead), within the prior bounds of the MCMC analysis. The
   likelihoods are copied from the production notebooks
   (`Codes/MCMC-Models/*/CC_SNe_BAO_FRB/`), including the same f(R) H(z)
   interpolation grids; M_B is profiled analytically (exact).
2. `mocks`: mock data sets generated at the ΛCDM best fit: CC with Gaussian
   errors, SNe Ia and BAO with their full covariances (Cholesky), FRBs with
   DM_cos drawn from the truncated Gaussian and DM_host from the log-normal
   at the best-fit host parameters (redshifts fixed). Each mock is re-fitted
   by all five models and Λ_LR = 2[ln L_max(alt) − ln L_max(ΛCDM)] stored.
3. `summary`: empirical p-values p = (N_{≥obs} + 1)/(N + 1), compared with the
   Wilks (χ²_ΔK) and Chernoff (½χ²_0 + ½χ²_1) asymptotic references.

## Usage

```bash
python lrt_bootstrap.py observed
python lrt_bootstrap.py mocks --n 5000 --nproc 80 --seed 1
python lrt_bootstrap.py summary
```

Outputs: `observed_fits.npz`, `mocks_seed*_*.npz`, `lrt_calibration_summary.txt`.
Each mock takes ~50 s on one core (5 models × several starting points).
