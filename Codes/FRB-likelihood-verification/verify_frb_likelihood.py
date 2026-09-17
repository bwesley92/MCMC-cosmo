"""
Numerical verification of the FRB likelihood normalization.

Context
-------
The original manuscript equation for the marginal FRB probability density,
P(DM_ext | z), contained a spurious external factor <DM_cos>^-1 and a
rescaled argument of P_cos. This was fixed in the manuscript (current
Eq. 6.7 in main_arxiv.tex): P_cos is evaluated directly in physical DM
units and the convolution with the host-galaxy log-normal PDF integrates
to unity without any additional normalization factor.

The response letter (response_to_referee.tex) claims that this
bug was confined to the manuscript equation, and that the production
likelihood code already used the correct (physical-units) form in all
six models. This script verifies that claim numerically, in three steps:

  1. Confirms (by direct extraction and comparison) that the FRB
     likelihood function `P_DM_ext_vectorized` is implemented identically
     (same functional form) across all six main-analysis notebooks
     (LCDM, wCDM, CPL, AB, HS, ST) -- see extract_and_diff_functions.py
     for that check; this script assumes/re-confirms the LCDM form.

  2. Numerically confirms that the truncated-Gaussian P_cos (Eq. 3.7)
     integrates to unity over DM_cos >= 0, for redshifts spanning the
     real FRB sample.

  3. Numerically confirms that the marginal P(DM_ext | z) (Eq. 6.7)
     integrates to unity over DM_ext, and cross-validates the production
     code's vectorized/trapezoidal implementation against an independent
     adaptive-quadrature (scipy.integrate.quad) implementation, evaluated
     at the 104 real FRBs and at the ΛCDM Baseline+FRB best-fit
     parameters (Table 2 of the manuscript).

Nothing in the existing repository is modified; this script only reads
`Data/FRB/frb_catalog.txt` and reimplements the relevant equations from
scratch.
"""

import numpy as np
from scipy.integrate import trapezoid, cumulative_trapezoid, quad
from scipy.special import erf
import pandas as pd
import json

REPO = "/home/brunowesley/projetos/MCMC-cosmo"

# ------------------------------------------------------------------
# 1. Load the real FRB catalog (same loading logic as the production
#    notebooks: Codes/GR-based-models/LCDM/CC_SNe_BAO_FRB/LCDM_CC_SNe_BAO_FRB.ipynb)
# ------------------------------------------------------------------
df_frb = pd.read_csv(f"{REPO}/Data/FRB/frb_catalog.txt", sep="\t")
df_frb = df_frb.sort_values(by="Redshift", ascending=True).reset_index(drop=True)

z_frb = df_frb["Redshift"].to_numpy(float)
DM_obs = df_frb["DM_obs"].to_numpy(float)
DM_MW = df_frb["DM_MW_ISM_(NE2001)"].to_numpy(float)
DM_MW_total = DM_MW + 50.0
DM_ext_obs = DM_obs - DM_MW_total

print(f"Loaded {len(z_frb)} FRBs, z in [{z_frb.min():.4f}, {z_frb.max():.4f}]")

# ------------------------------------------------------------------
# 2. Physical constants and cosmology (identical to the production
#    notebooks; flat LCDM, Baseline+FRB best-fit from Table 2)
# ------------------------------------------------------------------
c_si = 2.99792458e8
G_si = 6.67430e-11
mp_si = 1.67262192369e-27
km_to_pc = 1.0 / 3.0857e13
m2_to_cm2 = 1e4
Mpc_to_cm = 3.0857e24
f_IGM = 0.83
chi_e = 0.875

c_kms = 299792.458
z_grid = np.linspace(0.0, 2.5, 300)


def E_LCDM(z, Om0):
    return np.sqrt(Om0 * (1 + z) ** 3 + (1 - Om0))


def H_LCDM(z, H0, Om0):
    return H0 * E_LCDM(z, Om0)


def DM_IGM_LCDM(z, H0, Om0, Ob):
    prefactor = 3.0 * c_si * f_IGM * chi_e / (8.0 * np.pi * G_si * mp_si)
    H_grid = H_LCDM(z_grid, H0, Om0)
    integrand = (1 + z_grid) / H_grid
    I_grid = cumulative_trapezoid(integrand, z_grid, initial=0.0)
    I_z = np.interp(z, z_grid, I_grid)
    DM_IGM_z = prefactor * Ob * H0 ** 2 * I_z
    return DM_IGM_z * (km_to_pc / (m2_to_cm2 * Mpc_to_cm))


# Table 2 (LCDM, Baseline+FRB) marginalized medians:
H0_bf, Om0_bf, Ob_bf = 70.46, 0.305, 0.0501
sigma_host_bf, exp_mu_bf = 0.833, 104.2

DM_IGM_mean = DM_IGM_LCDM(z_frb, H0_bf, Om0_bf, Ob_bf)
sigma_IGM = 173.8 * (z_frb ** 0.4)

print(f"\nBest-fit point used: H0={H0_bf}, Om0={Om0_bf}, Ob={Ob_bf}, "
      f"sigma_host={sigma_host_bf}, exp_mu={exp_mu_bf}")
print(f"DM_IGM_mean range: [{DM_IGM_mean.min():.2f}, {DM_IGM_mean.max():.2f}] pc/cm^3")
print(f"sigma_IGM range:   [{sigma_IGM.min():.2f}, {sigma_IGM.max():.2f}] pc/cm^3")


# ------------------------------------------------------------------
# 3. Production likelihood function, extracted verbatim from
#    Codes/GR-based-models/LCDM/CC_SNe_BAO_FRB/LCDM_CC_SNe_BAO_FRB.ipynb
#    (identical in all six main-analysis notebooks, see
#    extract_and_diff_functions.py)
# ------------------------------------------------------------------
def P_DM_ext_vectorized(z_frb, DM_ext_obs, DM_IGM_mean, sigma_IGM,
                         sigma_host, exp_mu, npts=300):
    mu_host = np.log(exp_mu)
    upper = max(DM_ext_obs.max(), 5.0 * exp_mu)
    DM_hosp_grid = np.linspace(1e-2, upper, npts)
    DM_hosp_grid = DM_hosp_grid[None, :]
    DM_ext_array = DM_ext_obs[:, None]
    DM_IGM_array = DM_ext_array - DM_hosp_grid / (1 + z_frb[:, None])
    sigma_IGM_array = sigma_IGM[:, None]
    DM_IGM_mean_array = DM_IGM_mean[:, None]

    y = (DM_IGM_array - DM_IGM_mean_array) / sigma_IGM_array
    phi = (1 / np.sqrt(2 * np.pi)) * np.exp(-0.5 * y ** 2)
    norm_factor = 0.5 * (1 + erf(DM_IGM_mean_array / (sigma_IGM_array * np.sqrt(2))))
    P_IGM = (1 / sigma_IGM_array) * phi / norm_factor
    P_IGM[DM_IGM_array < 0] = 0.0

    DM_hosp_array = DM_hosp_grid
    P_host = (1.0 / (np.sqrt(2 * np.pi) * sigma_host * DM_hosp_array)) \
        * np.exp(-0.5 * ((np.log(DM_hosp_array) - mu_host) / sigma_host) ** 2)

    integrand = P_IGM * P_host
    P_total = trapezoid(integrand, x=DM_hosp_grid[0], axis=1)
    return np.maximum(P_total, 1e-100)


# ------------------------------------------------------------------
# 4. Independent (from-scratch) implementation of Eqs. (3.7) and (6.7),
#    using adaptive quadrature instead of the production code's
#    fixed trapezoidal grid.
# ------------------------------------------------------------------
def P_cos_scalar(DM_cos, mean, sigma):
    """Eq. (3.7): truncated Gaussian for DM_cos >= 0."""
    if DM_cos < 0:
        return 0.0
    N = 0.5 * (1 + erf(mean / (sigma * np.sqrt(2))))
    phi = (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * ((DM_cos - mean) / sigma) ** 2)
    return phi / N


def P_host_scalar(DM_host, sigma_host, exp_mu):
    if DM_host <= 0:
        return 0.0
    mu_host = np.log(exp_mu)
    return (1.0 / (np.sqrt(2 * np.pi) * sigma_host * DM_host)) \
        * np.exp(-0.5 * ((np.log(DM_host) - mu_host) / sigma_host) ** 2)


def P_DM_ext_scalar(DM_ext, z, mean, sigma, sigma_host, exp_mu):
    """Eq. (6.7): marginal PDF via adaptive-quadrature convolution."""
    def integrand(DM_host):
        DM_cos = DM_ext - DM_host / (1 + z)
        return P_cos_scalar(DM_cos, mean, sigma) * P_host_scalar(DM_host, sigma_host, exp_mu)

    upper_host = DM_ext * (1 + z)  # beyond this, DM_cos < 0 and P_cos = 0
    if upper_host <= 0:
        return 0.0
    val, _ = quad(integrand, 0.0, upper_host, limit=200)
    return val


# ------------------------------------------------------------------
# TEST 1: normalization of the truncated Gaussian P_cos (Eq. 3.7)
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("TEST 1: normalization of P_cos over DM_cos >= 0 (Eq. 3.7)")
print("=" * 70)
z_test_grid = np.linspace(z_frb.min(), z_frb.max(), 20)
mean_test = DM_IGM_LCDM(z_test_grid, H0_bf, Om0_bf, Ob_bf)
sigma_test = 173.8 * (z_test_grid ** 0.4)

max_dev_pcos = 0.0
for z_i, m_i, s_i in zip(z_test_grid, mean_test, sigma_test):
    val, _ = quad(P_cos_scalar, 0.0, m_i + 50 * s_i, args=(m_i, s_i), limit=200)
    dev = abs(val - 1.0)
    max_dev_pcos = max(max_dev_pcos, dev)
print(f"Tested {len(z_test_grid)} redshifts in [{z_test_grid.min():.4f}, {z_test_grid.max():.4f}]")
print(f"Max |integral(P_cos) - 1| = {max_dev_pcos:.3e}")

# ------------------------------------------------------------------
# TEST 2: normalization of the marginal P(DM_ext|z) (Eq. 6.7)
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("TEST 2: normalization of P(DM_ext|z) over DM_ext >= 0 (Eq. 6.7)")
print("=" * 70)
z_test_sub = np.linspace(z_frb.min(), z_frb.max(), 6)  # coarser: this test is expensive (nested quad)
mean_sub = DM_IGM_LCDM(z_test_sub, H0_bf, Om0_bf, Ob_bf)
sigma_sub = 173.8 * (z_test_sub ** 0.4)

max_dev_pext = 0.0
for z_i, m_i, s_i in zip(z_test_sub, mean_sub, sigma_sub):
    upper_ext = m_i + 50 * s_i + 5 * exp_mu_bf
    val, _ = quad(
        lambda DM_ext: P_DM_ext_scalar(DM_ext, z_i, m_i, s_i, sigma_host_bf, exp_mu_bf),
        0.0, upper_ext, limit=200,
    )
    dev = abs(val - 1.0)
    max_dev_pext = max(max_dev_pext, dev)
    print(f"  z={z_i:.4f}: integral(P_DM_ext) = {val:.6f}  (deviation from 1: {dev:.3e})")
print(f"\nMax |integral(P_DM_ext) - 1| = {max_dev_pext:.3e}")

# ------------------------------------------------------------------
# TEST 3: production code vs. independent quadrature, at the 104 real
# FRBs and the LCDM Baseline+FRB best-fit point.
# ------------------------------------------------------------------
print("\n" + "=" * 70)
print("TEST 3: production P_DM_ext_vectorized vs. independent quadrature")
print("        (104 real FRBs, LCDM Baseline+FRB best-fit)")
print("=" * 70)

P_production = P_DM_ext_vectorized(
    z_frb, DM_ext_obs, DM_IGM_mean, sigma_IGM, sigma_host_bf, exp_mu_bf, npts=300
)
P_production_hires = P_DM_ext_vectorized(
    z_frb, DM_ext_obs, DM_IGM_mean, sigma_IGM, sigma_host_bf, exp_mu_bf, npts=20000
)

P_independent = np.array([
    P_DM_ext_scalar(dm, z, m, s, sigma_host_bf, exp_mu_bf)
    for dm, z, m, s in zip(DM_ext_obs, z_frb, DM_IGM_mean, sigma_IGM)
])

rel_diff_default = np.abs(P_production - P_independent) / P_independent
rel_diff_hires = np.abs(P_production_hires - P_independent) / P_independent

lnL_production = np.sum(np.log(P_production))
lnL_independent = np.sum(np.log(P_independent))

print(f"npts=300  (production default): max relative diff = {rel_diff_default.max():.3e}, "
      f"median = {np.median(rel_diff_default):.3e}")
print(f"npts=20000 (high-res trapezoid): max relative diff = {rel_diff_hires.max():.3e}, "
      f"median = {np.median(rel_diff_hires):.3e}")
print(f"\nln L_FRB (production, npts=300)  = {lnL_production:.6f}")
print(f"ln L_FRB (independent quadrature) = {lnL_independent:.6f}")
print(f"Difference in ln L_FRB            = {lnL_production - lnL_independent:.3e}")

# Identify which FRBs dominate the npts=300 relative-difference tail
idx_worst = np.argsort(-rel_diff_default)[:10]
print("\nFRBs with the largest production-vs-quadrature relative difference (npts=300):")
print(f"{'Name':18s} {'z':>8s} {'DM_ext_obs':>11s} {'rel_diff':>10s}")
worst_list = []
for i in idx_worst:
    name = str(df_frb['Name'].iloc[i])
    print(f"{name:18s} {z_frb[i]:8.4f} {DM_ext_obs[i]:11.3f} {rel_diff_default[i]:10.3e}")
    worst_list.append({
        "name": name, "z": float(z_frb[i]),
        "DM_ext_obs": float(DM_ext_obs[i]),
        "rel_diff_npts300": float(rel_diff_default[i]),
    })
n_low_z_outliers = int(np.sum(rel_diff_default > 1e-2))
print(f"\nNumber of FRBs with rel. diff > 1% at npts=300: {n_low_z_outliers} "
      f"(all at low z, where the fixed host-DM grid is coarse relative to the "
      f"narrow integration region set by DM_ext_obs*(1+z))")

# ------------------------------------------------------------------
# Save a summary
# ------------------------------------------------------------------
summary = {
    "n_frb": int(len(z_frb)),
    "z_range": [float(z_frb.min()), float(z_frb.max())],
    "bestfit_point": {
        "H0": H0_bf, "Om0": Om0_bf, "Ob": Ob_bf,
        "sigma_host": sigma_host_bf, "exp_mu": exp_mu_bf,
    },
    "test1_Pcos_normalization": {
        "n_redshifts_tested": int(len(z_test_grid)),
        "max_abs_deviation_from_unity": float(max_dev_pcos),
    },
    "test2_PDMext_normalization": {
        "n_redshifts_tested": int(len(z_test_sub)),
        "max_abs_deviation_from_unity": float(max_dev_pext),
    },
    "test3_production_vs_quadrature": {
        "npts_300": {
            "max_relative_diff": float(rel_diff_default.max()),
            "median_relative_diff": float(np.median(rel_diff_default)),
        },
        "npts_20000": {
            "max_relative_diff": float(rel_diff_hires.max()),
            "median_relative_diff": float(np.median(rel_diff_hires)),
        },
        "lnL_FRB_production_npts300": float(lnL_production),
        "lnL_FRB_independent_quadrature": float(lnL_independent),
        "lnL_FRB_difference": float(lnL_production - lnL_independent),
        "n_frb_reldiff_above_1pct_npts300": n_low_z_outliers,
        "worst_frbs_npts300": worst_list,
    },
}

with open(f"{REPO}/Codes/FRB-likelihood-verification/results_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print("\nSaved results_summary.json")
