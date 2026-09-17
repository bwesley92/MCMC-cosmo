# Bayesian Evidence

Direct numerical estimate of the Bayesian evidence `ln Z` (in the
manuscript's notation; the code stores it as `lnE`) and the resulting
Bayes factors `ln BF` (relative to ΛCDM) for all six models, as a genuine
alternative to the BIC-based asymptotic approximation
`ln BF ≈ -ΔBIC/2` used in `../AIC-BIC-LRT/`.

## Method

We use the **k-th-nearest-neighbour (kNN) evidence estimator** of
Heavens et al. (2017), implemented by the vendored [`MCEvidence/`](./MCEvidence)
package. Unlike nested sampling or thermodynamic integration, this estimator
computes `ln Z` directly from an existing MCMC chain — no re-fitting, no
additional likelihood evaluations — which makes it applicable as a
post-processing step on the chains already produced under
`Codes/GR-based-models/` and `Codes/MG-based-models/`.

For each model, for both the Baseline (CC+SNe+BAO) and Baseline+FRB chains:

1. **Load** the flattened chain (`flat_samples_*.npy`) and its
   log-probability (`log_prob_*.npy`).
2. **Prior volume**: computed as the product of the (flat) prior widths for
   that chain's parameters, from the bounds in Table 1 of the manuscript.
3. **Thinning**: the kNN estimator assumes independent samples, but raw
   `emcee` chains are autocorrelated. Each chain is reshaped back to
   `(step, walker, dim)`, its integrated autocorrelation time `tau` is
   recomputed (`emcee.autocorr.integrated_time`), and the chain is thinned
   by `tau_max / 2` — the same criterion already used for the convergence
   diagnostics in the original fitting notebooks — leaving a set of
   quasi-independent posterior samples.
4. **Evidence estimate**: the thinned `(weight, -log_prob, samples)` chain is
   passed to `MCEvidence(..., kmax=k+1)` with `k = 1` (first nearest
   neighbour). The returned `ln Z` (the code's `lnE` variable) carries an
   analytic uncertainty `sigma_lnE = 1 / sqrt(N*k + 1)`, with `N` the number
   of thinned samples.
5. **Bayes factor**: `ln BF(model) = ln Z(model) - ln Z(ΛCDM)`, with
   uncertainties combined in quadrature (the two evidence estimates are
   independent), interpreted on the **Jeffreys scale** (`|ln BF| < 1`: not
   worth a mention; `< 2.5`: substantial; `< 5`: strong; `≥ 5`: decisive).

This supplements the AIC/BIC/LRT comparison with a direct, genuine
evidence-based comparison, and explains why the two can disagree: the
BIC approximation assumes a large-sample Gaussian posterior, which is
less reliable for parameters only weakly constrained by the data, e.g.
`(w0, wa)` in CPL or `b` in AB.

## Files

- `<Model>/<model>_evidence.ipynb` (one per model: `AB`, `CPL`, `HS`, `LCDM`,
  `ST`, `wCDM`) — runs the five steps above for that model's Baseline and
  Baseline+FRB chains; saves `evidence_<model>.txt` (human-readable) and
  `evidence_<model>.npz` (arrays `lnE`, `sigma_lnE`, `ndim`, `N_used`,
  `priorvolume`, for both datasets).
- `BayesFactors_summary.ipynb` — loads all six `evidence_*.npz` files,
  computes the Bayes factors relative to ΛCDM for both datasets, compares
  the Baseline+FRB values against the BIC-based approximation from
  `../AIC-BIC-LRT/model_comparison.txt`, and saves `BayesFactors_summary.txt`
  / `.npz`.
- `MCEvidence/` — the third-party kNN evidence-estimator package (MIT
  license), vendored so the evidence notebooks can `import` it directly;
  not modified from upstream.

## How to reproduce

1. Run the relevant `*_CC_SNe_BAO.ipynb` and `*_CC_SNe_BAO_FRB.ipynb` chains
   under `Codes/GR-based-models/` / `Codes/MG-based-models/` first (the
   evidence notebooks only read `flat_samples_*.npy` / `log_prob_*.npy`,
   they don't regenerate them).
2. Run each `<Model>/<model>_evidence.ipynb`.
3. Run `BayesFactors_summary.ipynb`.

## References

1. Heavens, A., Fantaye, Y., Mootoovaloo, A., Eggers, H., Hosenie, Z.,
   Kroon, S., & Sellentin, E. (2017), *Marginal Likelihoods from Monte Carlo
   Markov Chains*, [arXiv:1704.03472](https://arxiv.org/abs/1704.03472) —
   the kNN evidence-estimator method used here.
2. `MCEvidence` package: [github.com/yabebalFantaye/MCEvidence](https://github.com/yabebalFantaye/MCEvidence)
   — the reference implementation vendored in [`MCEvidence/`](./MCEvidence).
3. Manuscript, Sec. "Model Comparison" (Table 10) — the reported
   constraints and the discussion of the disagreement with the
   BIC-based approximation.
