#!/usr/bin/env python3
"""
Parametric-bootstrap calibration of the likelihood-ratio test (LRT).

Null hypothesis: flat LCDM.  Mock CC + SNe Ia + BAO + FRB data sets are
generated at the LCDM maximum-likelihood point of the real data, and each
mock is re-fitted (maximum likelihood, within the same prior bounds used
in the MCMC analysis) by LCDM, wCDM, CPL, Hu-Sawicki (n=1) and
Starobinsky (n=1).  The distribution of

    Lambda_LR = 2 [ ln L_max(alt) - ln L_max(LCDM) ]

under the null is then compared with the observed value, which is
re-computed here with the very same optimizer.

The likelihoods are copied from the production notebooks
(Codes/MCMC-Models/*/CC_SNe_BAO_FRB/*.ipynb): same redshift grids,
same f(R) H(z) interpolation grids, same FRB convolution.  The only
difference is that M_B is profiled analytically (the SNe Ia chi^2 is
quadratic in M_B), which gives the exact maximum over M_B.

Usage
-----
    python lrt_bootstrap.py observed                # real-data fits
    python lrt_bootstrap.py mocks --n 1000 --nproc 48 --seed 1
    python lrt_bootstrap.py summary                 # p-values etc.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import sys
import time
import argparse
import numpy as np
import pandas as pd
from scipy.integrate import trapezoid, cumulative_trapezoid
from scipy.interpolate import RegularGridInterpolator
from scipy.special import erf
from scipy.optimize import minimize
from scipy import stats

ROOT = "/home/brunowesley/projetos/MCMC-cosmo"
OUT = os.path.join(ROOT, "Codes/Model-comparison/LRT-calibration")

# ====================================================================
# Data (identical to the production notebooks)
# ====================================================================
df_cc = pd.read_csv(f"{ROOT}/Data/CC/CC_Hz_data.txt", sep=r"\s+")
z_cc = df_cc["z_cc"].to_numpy(float)
H_obs = df_cc["H_cc"].to_numpy(float)
H_err = df_cc["H_cc_err"].to_numpy(float)

df_sne = pd.read_csv(f"{ROOT}/Data/SNe/Pantheon+SH0ES.dat", sep=r"\s+",
                     comment="#", engine="python")
raw_cov = np.loadtxt(f"{ROOT}/Data/SNe/Pantheon+SH0ES_STAT+SYS.cov", skiprows=1)
Cov_sne_full = raw_cov.reshape((1701, 1701))
mask = df_sne["zCMB"].values > 0.01
df_cut = df_sne[mask].reset_index(drop=True)
Cov_sne = Cov_sne_full[mask, :][:, mask]
z_sne = df_cut["zCMB"].values
mB_obs = df_cut["m_b_corr"].values
Cinv_sne = np.linalg.inv(Cov_sne)
L_sne = np.linalg.cholesky(Cov_sne)
one = np.ones_like(z_sne)
Cinv1 = Cinv_sne @ one
oneCinv1 = one @ Cinv1

df_bao = pd.read_csv(f"{ROOT}/Data/BAO/desi_gaussian_bao_ALL_GCcomb_mean.txt",
                     sep=r"\s+", comment="#", names=["z", "value", "quantity"])
z_bao = df_bao["z"].values
d_obs = df_bao["value"].values
q_bao = df_bao["quantity"].str.lower().values
cov_bao = np.loadtxt(f"{ROOT}/Data/BAO/desi_gaussian_bao_ALL_GCcomb_cov.txt")
Cinv_bao = np.linalg.inv(cov_bao)
L_bao = np.linalg.cholesky(cov_bao)
mask_dv = np.array(["dv" in q for q in q_bao])
mask_dm = np.array(["dm" in q for q in q_bao])
mask_dh = np.array(["dh" in q for q in q_bao])

df_frb = pd.read_csv(f"{ROOT}/Data/FRB/frb_catalog.txt", sep=r"\t", engine="python")
df_frb = df_frb.sort_values(by="Redshift", ascending=True).reset_index(drop=True)
z_frb = df_frb["Redshift"].to_numpy(float)
DM_ext_real = (df_frb["DM_obs"].to_numpy(float)
               - (df_frb["DM_MW_ISM_(NE2001)"].to_numpy(float) + 50.0))
sigma_cos_frb = 173.8 * z_frb**0.4

N_TOT = len(z_cc) + len(z_sne) + len(z_bao) + len(z_frb)

# ====================================================================
# Constants
# ====================================================================
c_kms = 299792.458
c_si = 2.99792458e8
G_si = 6.67430e-11
mp_si = 1.67262192369e-27
km_to_pc = 1.0 / 3.0857e13
m2_to_cm2 = 1e4
Mpc_to_cm = 3.0857e24
f_d, chi_e, Neff = 0.83, 0.875, 3.046
DM_PREF = (3.0 * c_si * f_d * chi_e / (8.0 * np.pi * G_si * mp_si)
           * (km_to_pc / (m2_to_cm2 * Mpc_to_cm)))

# redshift grids: GR models use [0, 2.5], f(R) models [0, 4] (as in production)
Z_GR = np.linspace(0.0, 2.5, 300)
Z_FR = np.linspace(0.0, 4.0, 300)

# ====================================================================
# f(R) H(z) interpolators (production grids)
# ====================================================================
_HS = None
_ST = None


def load_fr_grids():
    global _HS, _ST
    H0p = np.linspace(50.0, 90.0, 30)
    Omp = np.linspace(0.1, 0.6, 30)
    a_hs = np.load(f"{ROOT}/Codes/MCMC-Models/MG-based/Hu-Sawicki/CC/H_array_HS_n1_1k_grid.npy")
    _HS = RegularGridInterpolator((H0p, Omp, np.linspace(0.0, 300.0, 1000), Z_FR),
                                  a_hs, method="linear", bounds_error=False,
                                  fill_value=np.nan)
    a_st = np.load(f"{ROOT}/Codes/MCMC-Models/MG-based/Starobinsky/CC/H_array_ST_n1_30_grid.npy")
    _ST = RegularGridInterpolator((H0p, Omp, np.linspace(1.54, 10.0, 30), Z_FR),
                                  a_st, method="linear", bounds_error=False,
                                  fill_value=np.nan)


def _interp(itp, z, H0, Om0, p):
    pts = np.column_stack([np.full(len(z), H0), np.full(len(z), Om0),
                           np.full(len(z), p), z])
    return itp(pts)


# ====================================================================
# Models: H on the model's own z grid and at arbitrary z
# ====================================================================
def E_w(z, Om0, w0, wa):
    return np.sqrt(Om0 * (1 + z)**3 + (1 - Om0) * (1 + z)**(3 * (1 + w0 + wa))
                   * np.exp(-3 * wa * z / (1 + z)))


MODELS = {
    # name: (parameter names (excluding MB), bounds)
    "LCDM": (["H0", "Om0", "Ob", "sh", "emu"],
             [(50, 90), (0.1, 0.6), (0.01, 0.1), (0.2, 2.0), (20, 200)]),
    "wCDM": (["H0", "Om0", "Ob", "w", "sh", "emu"],
             [(50, 90), (0.1, 0.6), (0.01, 0.1), (-3, 1), (0.2, 2.0), (20, 200)]),
    "CPL": (["H0", "Om0", "Ob", "w0", "wa", "sh", "emu"],
            [(50, 90), (0.1, 0.6), (0.01, 0.1), (-3, 1), (-3, 2), (0.2, 2.0), (20, 200)]),
    "HS": (["H0", "Om0", "Ob", "mu", "sh", "emu"],
           [(50, 90), (0.1, 0.6), (0.01, 0.1), (0, 300), (0.2, 2.0), (20, 200)]),
    "ST": (["H0", "Om0", "Ob", "lbd", "sh", "emu"],
           [(50, 90), (0.1, 0.6), (0.01, 0.1), (1.54, 10), (0.2, 2.0), (20, 200)]),
}


def split(model, th):
    """Return (H0, Om0, Ob, extra-tuple, sigma_host, exp_mu)."""
    if model == "LCDM":
        H0, Om0, Ob, sh, emu = th
        ex = ()
    elif model == "CPL":
        H0, Om0, Ob, w0, wa, sh, emu = th
        ex = (w0, wa)
    else:
        H0, Om0, Ob, x, sh, emu = th
        ex = (x,)
    return H0, Om0, Ob, ex, sh, emu


def H_of(model, z, H0, Om0, ex):
    if model == "LCDM":
        return H0 * E_w(z, Om0, -1.0, 0.0)
    if model == "wCDM":
        return H0 * E_w(z, Om0, ex[0], 0.0)
    if model == "CPL":
        return H0 * E_w(z, Om0, ex[0], ex[1])
    if model == "HS":
        return _interp(_HS, z, H0, Om0, ex[0])
    if model == "ST":
        return _interp(_ST, z, H0, Om0, ex[0])


def predictions(model, H0, Om0, Ob, ex):
    """Theory vectors: H(z_cc), mu(z_sne) [without M_B], d_BAO, <DM_cos>(z_frb)."""
    zg = Z_FR if model in ("HS", "ST") else Z_GR
    Hg = H_of(model, zg, H0, Om0, ex)
    if not np.all(np.isfinite(Hg)) or np.any(Hg <= 0):
        return None
    # CC
    Hcc = H_of(model, z_cc, H0, Om0, ex)
    # SNe: dL via cumulative trapezoid on the model grid
    I = cumulative_trapezoid(1.0 / Hg, zg, initial=0.0)
    dL = c_kms * (1 + z_sne) * np.interp(z_sne, zg, I)
    mu = 5.0 * np.log10(dL) + 25.0
    # BAO (comoving distance with explicit trapezoid, as in production)
    integ = 1.0 / Hg
    Ib = np.empty_like(zg)
    Ib[0] = 0.0
    Ib[1:] = np.cumsum(0.5 * (integ[1:] + integ[:-1]) * np.diff(zg))
    dC = c_kms * np.interp(z_bao, zg, Ib)
    Hb = H_of(model, z_bao, H0, Om0, ex)
    h2 = (H0 / 100.0)**2
    rd = (147.05 * ((Ob * h2 / 0.02236)**-0.13) * ((Om0 * h2 / 0.1432)**-0.23)
          * ((Neff / 3.04)**-0.1))
    dm = dC / rd
    dh = (c_kms / Hb) / rd
    dv = (z_bao * dm**2 * dh)**(1 / 3)
    d = np.empty_like(z_bao)
    d[mask_dv] = dv[mask_dv]
    d[mask_dm] = dm[mask_dm]
    d[mask_dh] = dh[mask_dh]
    # FRB mean cosmic DM
    Id = cumulative_trapezoid((1 + zg) / Hg, zg, initial=0.0)
    DMc = DM_PREF * Ob * H0**2 * np.interp(z_frb, zg, Id)
    out = (Hcc, mu, d, DMc)
    if not all(np.all(np.isfinite(v)) for v in out):
        return None
    return out


def P_DM_ext(DM_ext, DM_cos_mean, sigma_host, exp_mu, npts=300):
    """Production FRB marginal likelihood (truncated Gaussian x log-normal)."""
    mu_host = np.log(exp_mu)
    upper = max(DM_ext.max(), 5.0 * exp_mu)
    g = np.linspace(1e-2, upper, npts)
    DMc = DM_ext[:, None] - g[None, :] / (1 + z_frb[:, None])
    sc = sigma_cos_frb[:, None]
    m = DM_cos_mean[:, None]
    y = (DMc - m) / sc
    phi = np.exp(-0.5 * y**2) / np.sqrt(2 * np.pi)
    norm = 0.5 * (1 + erf(m / (sc * np.sqrt(2))))
    Pc = phi / sc / norm
    Pc[DMc < 0] = 0.0
    Ph = np.exp(-0.5 * ((np.log(g) - mu_host) / sigma_host)**2) / (np.sqrt(2 * np.pi) * sigma_host * g)
    return np.maximum(trapezoid(Pc * Ph[None, :], x=g, axis=1), 1e-100)


def lnlike(model, th, data, return_MB=False):
    H0, Om0, Ob, ex, sh, emu = split(model, th)
    pr = predictions(model, H0, Om0, Ob, ex)
    if pr is None:
        return -np.inf
    Hcc, mu, d, DMc = pr
    Hd, mBd, dd, DMd = data
    lnL = -0.5 * np.sum(((Hd - Hcc) / H_err)**2)
    # SNe with M_B profiled analytically (exact maximum over M_B)
    r = mBd - mu
    MB = (Cinv1 @ r) / oneCinv1
    MB = min(max(MB, -20.0), -18.0)
    rr = r - MB
    lnL += -0.5 * (rr @ Cinv_sne @ rr)
    rb = dd - d
    lnL += -0.5 * (rb @ Cinv_bao @ rb)
    lnL += np.sum(np.log(P_DM_ext(DMd, DMc, sh, emu)))
    if return_MB:
        return lnL, MB
    return lnL


# ====================================================================
# Maximum-likelihood fits (bounded, multi-start)
# ====================================================================
def fit(model, data, starts):
    names, bounds = MODELS[model]
    lo = np.array([b[0] for b in bounds], float)
    hi = np.array([b[1] for b in bounds], float)

    def to_u(x):
        return (np.asarray(x) - lo) / (hi - lo)

    def to_x(u):
        return lo + np.clip(u, 0, 1) * (hi - lo)

    def nll(u):
        v = lnlike(model, to_x(u), data)
        return 1e30 if not np.isfinite(v) else -v

    best = (np.inf, None)
    ub = [(0.0, 1.0)] * len(lo)
    for s in starts:
        u = np.clip(to_u(s), 0, 1)
        r = minimize(nll, u, method="Powell", bounds=ub,
                     options={"xtol": 1e-7, "ftol": 1e-10, "maxfev": 20000})
        r = minimize(nll, r.x, method="Nelder-Mead", bounds=ub,
                     options={"xatol": 1e-8, "fatol": 1e-9, "maxfev": 20000,
                              "adaptive": True})
        r = minimize(nll, r.x, method="Powell", bounds=ub,
                     options={"xtol": 1e-8, "ftol": 1e-12, "maxfev": 20000})
        if r.fun < best[0]:
            best = (r.fun, to_x(r.x))
    return -best[0], best[1]


def starts_for(model, base):
    """Starting points: the LCDM best fit `base` plus a grid in the extra parameter(s)."""
    H0, Om0, Ob, sh, emu = base
    if model == "LCDM":
        return [base, [68.0, 0.30, 0.050, 0.85, 100.0]]
    if model == "wCDM":
        return [[H0, Om0, Ob, w, sh, emu] for w in (-1.0, -0.88, -1.12)]
    if model == "CPL":
        return [[H0, Om0, Ob, w0, wa, sh, emu]
                for (w0, wa) in ((-1.0, 0.0), (-0.84, -0.3), (-0.9, 0.3), (-1.1, -0.5))]
    if model == "HS":
        return [[H0, Om0, Ob, m, sh, emu] for m in (0.0, 30.0, 75.0, 140.0)]
    if model == "ST":
        return [[H0, Om0, Ob, l, sh, emu] for l in (1.6, 2.1, 3.5, 9.9)]


ALT = ["wCDM", "CPL", "HS", "ST"]


def fit_all(data, base_start=None):
    res = {}
    lnL0, th0 = fit("LCDM", data, starts_for("LCDM", base_start)
                    if base_start is not None else
                    [[68.0, 0.30, 0.050, 0.85, 100.0], [70.5, 0.305, 0.050, 0.83, 104.0]])
    res["LCDM"] = (lnL0, th0)
    for m in ALT:
        lnL, th = fit(m, data, starts_for(m, th0))
        res[m] = (lnL, th)
    return res


# ====================================================================
# Mock generation under the LCDM null
# ====================================================================
def make_mock(theta_null, MB_null, rng):
    H0, Om0, Ob, sh, emu = theta_null
    Hcc, mu, d, DMc = predictions("LCDM", H0, Om0, Ob, ())
    Hm = Hcc + H_err * rng.standard_normal(len(z_cc))
    mBm = mu + MB_null + L_sne @ rng.standard_normal(len(z_sne))
    dm = d + L_bao @ rng.standard_normal(len(z_bao))
    a = (0.0 - DMc) / sigma_cos_frb
    dcos = stats.truncnorm.rvs(a, np.inf, loc=DMc, scale=sigma_cos_frb, random_state=rng)
    dhost = np.exp(np.log(emu) + sh * rng.standard_normal(len(z_frb)))
    DMm = dcos + dhost / (1 + z_frb)
    return (Hm, mBm, dm, DMm)


# ====================================================================
# Drivers
# ====================================================================
REAL = (H_obs, mB_obs, d_obs, DM_ext_real)


def run_observed():
    load_fr_grids()
    t = time.time()
    res = fit_all(REAL)
    lnL0, th0 = res["LCDM"]
    _, MB0 = lnlike("LCDM", th0, REAL, return_MB=True)
    print(f"[observed] done in {time.time()-t:.1f} s")
    out = {"theta_LCDM": th0, "MB_LCDM": MB0}
    for m, (lnL, th) in res.items():
        print(f"  {m:5s} lnLmax = {lnL:.4f}   theta = {np.round(th, 4)}")
        out[f"lnL_{m}"] = lnL
        out[f"theta_{m}"] = th
        if m != "LCDM":
            print(f"        Lambda_LR = {2*(lnL-lnL0):.4f}")
    np.savez(os.path.join(OUT, "observed_fits.npz"), **out)


def _one_mock(args):
    i, seed, th0, MB0 = args
    rng = np.random.default_rng([seed, i])
    data = make_mock(th0, MB0, rng)
    res = fit_all(data, base_start=th0)
    return i, np.array([res[m][0] for m in ["LCDM"] + ALT]), \
        np.concatenate([res[m][1] for m in ["LCDM"] + ALT])


def run_mocks(n, nproc, seed, start=0):
    from multiprocessing import get_context
    obs = np.load(os.path.join(OUT, "observed_fits.npz"))
    th0, MB0 = obs["theta_LCDM"], float(obs["MB_LCDM"])
    load_fr_grids()  # loaded before fork -> shared copy-on-write
    tag = f"seed{seed}_{start}_{start+n}"
    fn = os.path.join(OUT, f"mocks_{tag}.npz")
    tasks = [(i, seed, th0, MB0) for i in range(start, start + n)]
    lnL = np.full((n, 1 + len(ALT)), np.nan)
    TH = [None] * n
    t = time.time()
    done = 0
    with get_context("fork").Pool(nproc) as pool:
        for i, l, th in pool.imap_unordered(_one_mock, tasks, chunksize=1):
            lnL[i - start] = l
            TH[i - start] = th
            done += 1
            if done % max(1, n // 50) == 0 or done == n:
                el = time.time() - t
                print(f"  {done}/{n} mocks  {el/60:.1f} min  ETA {el/done*(n-done)/60:.1f} min",
                      flush=True)
                np.savez(fn, lnL=lnL, models=np.array(["LCDM"] + ALT), seed=seed,
                         start=start, theta_null=th0, MB_null=MB0)
    np.savez(fn, lnL=lnL, models=np.array(["LCDM"] + ALT), seed=seed, start=start,
             theta_null=th0, MB_null=MB0,
             theta=np.array([x for x in TH if x is not None]))
    print(f"[mocks] saved {fn}")


def run_summary():
    import glob
    obs = np.load(os.path.join(OUT, "observed_fits.npz"))
    files = sorted(glob.glob(os.path.join(OUT, "mocks_*.npz")))
    L = np.concatenate([np.load(f)["lnL"] for f in files])
    L = L[np.all(np.isfinite(L), axis=1)]
    n = len(L)
    print(f"{n} valid mocks from {len(files)} file(s)")
    dk = {"wCDM": 1, "CPL": 2, "HS": 1, "ST": 1}
    lines = []
    for j, m in enumerate(ALT, start=1):
        lam = 2 * (L[:, j] - L[:, 0])
        lam_obs = 2 * (float(obs[f"lnL_{m}"]) - float(obs["lnL_LCDM"]))
        k = int(np.sum(lam >= lam_obs))
        p_emp = (k + 1) / (n + 1)
        p_wilks = stats.chi2.sf(lam_obs, dk[m])
        p_chern = 0.5 * stats.chi2.sf(lam_obs, 1) if dk[m] == 1 else np.nan
        lp = np.clip(lam, 0, None)
        frac0 = np.mean(lam < 1e-3)
        # effective scaled-chi2 fit to the positive part: lam | lam>0 ~ a * chi2_nu
        pos = lam[lam > 1e-3]
        nu_eff, _, a_eff = stats.chi2.fit(pos, floc=0)
        p_fit = (1 - frac0) * stats.chi2.sf(lam_obs / a_eff, nu_eff)
        q = np.percentile(lam, [50, 90, 95, 99])
        s = (f"{m:5s} Lobs={lam_obs:7.3f}  N>=obs={k:4d}/{n}  p_emp={p_emp:.2e}  "
             f"p_Wilks(chi2_{dk[m]})={p_wilks:.2e}  p_Chernoff={p_chern:.2e}  "
             f"P(L=0)={frac0:.3f}  fit: nu={nu_eff:.2f} a={a_eff:.2f} p_fit={p_fit:.2e}  "
             f"mean={lam.mean():.3f}  q50/90/95/99={np.round(q, 2)}  min={lam.min():.3f}")
        print(s)
        lines.append(s)
    with open(os.path.join(OUT, "lrt_calibration_summary.txt"), "w") as f:
        f.write(f"{n} mocks\n" + "\n".join(lines) + "\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["observed", "mocks", "summary", "timing"])
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--nproc", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--start", type=int, default=0)
    a = ap.parse_args()
    if a.mode == "observed":
        run_observed()
    elif a.mode == "mocks":
        run_mocks(a.n, a.nproc, a.seed, a.start)
    elif a.mode == "summary":
        run_summary()
    elif a.mode == "timing":
        load_fr_grids()
        for m, th in [("LCDM", [70.5, 0.305, 0.05, 0.83, 104]),
                      ("HS", [70.3, 0.282, 0.0514, 77.5, 0.857, 96.7])]:
            t = time.time()
            for _ in range(200):
                v = lnlike(m, th, REAL)
            print(m, v, f"{(time.time()-t)/200*1e3:.2f} ms/eval")
        # check HS(mu=0) and ST(lbd=10) vs LCDM
        th = [70.46, 0.305, 0.0501, 0.833, 104.2]
        print("LCDM", lnlike("LCDM", th, REAL))
        print("HS mu=0", lnlike("HS", th[:3] + [0.0] + th[3:], REAL))
        print("ST lbd=10", lnlike("ST", th[:3] + [10.0] + th[3:], REAL))
