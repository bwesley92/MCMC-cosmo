# Data

Every file in this folder is third-party observational data, redistributed
here **only** for reproducibility, in the exact form it is read by the
notebooks in `Codes/`. All credit belongs to the original authors and
collaborations below — please cite them (not this repository) when using
this data.

## Cosmic Chronometers — `CC/CC_Hz_data.txt`

35 `H(z)` measurements (`z`, `H(z)`, `sigma_H`), each compiled from the
individual differential-age study that reported it:

- Simon, Verde & Jimenez (2005), *Constraints on the redshift dependence of
  the dark energy potential*, [arXiv:astro-ph/0412269](https://arxiv.org/abs/astro-ph/0412269)
- Stern, Jimenez, Verde, Kamionkowski & Stanford (2010), *Cosmic
  chronometers: constraining the equation of state of dark energy. I: H(z)
  measurements*, JCAP 02, 008,
  [DOI:10.1088/1475-7516/2010/02/008](https://doi.org/10.1088/1475-7516/2010/02/008)
- Moresco et al. (2012), *Improved constraints on the expansion rate of the
  Universe up to z~1.1 from the spectroscopic evolution of cosmic
  chronometers*, JCAP 08, 006,
  [DOI:10.1088/1475-7516/2012/08/006](https://doi.org/10.1088/1475-7516/2012/08/006)
- Zhang, Zhang, Yuan, Liu, Zhang & Sun (2014), *Four new observational H(z)
  data from luminous red galaxies in the Sloan Digital Sky Survey data
  release seven*, RAA 14, 1221,
  [DOI:10.1088/1674-4527/14/10/002](https://doi.org/10.1088/1674-4527/14/10/002)
- Moresco (2015), *Raising the bar: new constraints on the Hubble
  parameter with cosmic chronometers at z~2*,
  [arXiv:1503.01116](https://arxiv.org/abs/1503.01116)
- Moresco et al. (2016), *A 6% measurement of the Hubble parameter at
  z~0.45: direct evidence of the epoch of cosmic re-acceleration*,
  [arXiv:1601.01701](https://arxiv.org/abs/1601.01701)
- Ratsimbazafy et al. (2017), *Age-dating luminous red galaxies observed
  with the Southern African Large Telescope*, MNRAS 467, 3239,
  [DOI:10.1093/mnras/stx301](https://doi.org/10.1093/mnras/stx301)
- Borghi, Moresco & Cimatti (2022), *Toward a Better Understanding of
  Cosmic Chronometers: A New Measurement of H(z) at z~0.7*, ApJL 928,
  [DOI:10.3847/2041-8213/ac3fb2](https://doi.org/10.3847/2041-8213/ac3fb2)
- Jimenez, Moresco, Verde & Wandelt (2023), *Cosmic chronometers with
  photometry: a new path to H(z)*, JCAP 11, 047,
  [DOI:10.1088/1475-7516/2023/11/047](https://doi.org/10.1088/1475-7516/2023/11/047)
- Jiao, Borghi, Moresco & Zhang (2023), *New Observational H(z) Data from
  Full-spectrum Fitting of Cosmic Chronometers in the LEGA-C Survey*,
  [arXiv:2205.05701](https://arxiv.org/abs/2205.05701)
- Tomasetti et al. (2023), *A new measurement of the expansion history of
  the Universe at z=1.26 with cosmic chronometers in VANDELS*,
  [arXiv:2305.16387](https://arxiv.org/abs/2305.16387)

## Type Ia Supernovae — `SNe/Pantheon+SH0ES.dat`, `Pantheon+SH0ES_STAT+SYS.cov`

The Pantheon+ compilation (1701 light curves) with the SH0ES Cepheid
distance-anchor calibration:

- Scolnic et al. (2022), *The Pantheon+ Analysis: The Full Data Set and
  Light-curve Release*, ApJ 938, 113,
  [DOI:10.3847/1538-4357/ac8b7a](https://doi.org/10.3847/1538-4357/ac8b7a)
- Riess et al. (2022), *A Comprehensive Measurement of the Local Value of
  the Hubble Constant... from HST and the SH0ES Team*, ApJL 934, L7,
  [DOI:10.3847/2041-8213/ac5c5b](https://doi.org/10.3847/2041-8213/ac5c5b)

Official public data release:
[github.com/PantheonPlusSH0ES/DataRelease](https://github.com/PantheonPlusSH0ES/DataRelease)

## Baryon Acoustic Oscillations — `BAO/desi_gaussian_bao_ALL_GCcomb_mean.txt`, `..._cov.txt`

13 Gaussianized BAO distance measurements from DESI Data Release 2:

- Abdul Karim et al. (2025) [DESI Collaboration], *DESI DR2 results. II.
  Measurements of baryon acoustic oscillations and cosmological
  constraints*, Phys. Rev. D 112, 8,
  [DOI:10.1103/tr6y-kpc6](https://doi.org/10.1103/tr6y-kpc6)

These exact files are the DESI DR2 all-tracer Gaussian BAO data vector and
covariance, mirrored (unmodified) from
[CobayaSampler/bao_data](https://github.com/CobayaSampler/bao_data/tree/master/desi_bao_dr2)
(`desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_{mean,cov}.txt`), a git
repository of likelihood inputs maintained for use with the Cobaya
sampling framework. The official DESI collaboration data portal is
[data.desi.lbl.gov](https://data.desi.lbl.gov).

## Fast Radio Bursts — `FRB/frb_catalog.txt`

104 localized FRBs with spectroscopically confirmed host-galaxy redshifts.
Each row carries its own `Reference` column crediting the specific
localization study; the sample as a whole draws on several facilities and
surveys, principally:

- Amiri et al. (2025) [CHIME/FRB Collaboration], *A Catalog of Local
  Universe Fast Radio Bursts from CHIME/FRB and the KKO*, ApJS 280,
  [DOI:10.3847/1538-4365/addbda](https://doi.org/10.3847/1538-4365/addbda)
- Shannon et al. (2025), *The commensal real-time ASKAP fast transient
  incoherent-sum survey (CRAFT)*, PASA 42,
  [DOI:10.1017/pasa.2025.8](https://doi.org/10.1017/pasa.2025.8)
- Sharma et al. (2024), Law et al. (2024) & Connor et al. (2025)
  [DSA-110], [arXiv:2409.16964](https://arxiv.org/abs/2409.16964),
  [arXiv:2307.03344](https://arxiv.org/abs/2307.03344),
  [arXiv:2409.16952](https://arxiv.org/abs/2409.16952)

The `DM_MW_ISM_(NE2001)` column is computed from the Galactic
electron-density model of Cordes & Lazio (2002),
[arXiv:astro-ph/0207156](https://arxiv.org/abs/astro-ph/0207156) (NE2001).
