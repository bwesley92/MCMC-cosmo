# MCMC R2-AB Cosmological Model

A comprehensive Python implementation for constraining cosmological parameters using Markov Chain Monte Carlo (MCMC) methods with the R2-AB modified gravity model.

## Overview

This project implements a Bayesian parameter estimation framework to constrain cosmological parameters using multiple observational datasets:
- **SNe Ia**: Type Ia Supernovae data (Pantheon+SH0ES)
- **CC**: Cosmic Chronometer data
- **RSD**: Redshift Space Distortions data

The code solves modified Friedmann equations for the R2-AB gravity model and performs statistical analysis to determine the best-fit parameters and their uncertainties.

## Features

- **Modified Gravity Model**: R2-AB gravity implementation with hyperbolic tangent function
- **Multi-Dataset Analysis**: Simultaneous fitting of SNe, CC, and RSD observational data
- **MCMC Sampling**: Efficient parameter space exploration using `emcee` ensemble sampler
- **Parallel Processing**: Utilizes multiprocessing for faster computation
- **Visualization**: Corner plots for posterior distributions

## Requirements

### Python Dependencies

```bash
pip install numpy pandas scipy matplotlib emcee corner
```

### Required Libraries
- `numpy` - Numerical computations
- `pandas` - Data handling
- `scipy` - Integration and optimization
- `matplotlib` - Plotting and visualization
- `emcee` - MCMC ensemble sampler
- `corner` - Corner plot generation

## Data Requirements

The code expects the following data files:

1. **Pantheon+SH0ES SNe Data**:
   - Path: `/home/usuario/Downloads/DataRelease-main/Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES.dat`
   - Covariance matrix: `/home/usuario/Downloads/DataRelease-main/Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES_STAT+SYS.cov`

2. **Cosmic Chronometer Data**: Embedded in the code
3. **RSD Data**: Embedded in the code

**Note**: Update the file paths in lines 90-98 to match your data locations.

## Model Parameters

The code fits for 5 cosmological parameters:

| Parameter | Description | Prior Range |
|-----------|-------------|-------------|
| `H0` | Hubble constant (km/s/Mpc) | 54 - 76 |
| `Ω_m0` | Matter density parameter | 0.1 - 0.5 |
| `σ8(0)` | Amplitude of matter fluctuations | 0.7 - 0.9 |
| `Mb` | Absolute magnitude of SNe | -20.2 - -19.0 |
| `b` | AB model parameter | 1.6 - 12.0 |

## Physical Models

### 1. Hubble Parameter H(t)
Solves the modified Friedmann equation for the R2-AB model with initial conditions set at redshift z ≈ 4 (t = 0.2).

### 2. Growth Function fs8(z)
Computes the growth rate of cosmic structure formation including modifications from the R2-AB model.

### 3. Luminosity Distance
Calculates the luminosity distance to supernovae using numerical integration of the inverse Hubble parameter.

## Usage

### Basic Execution

```python
python "MCMC_AB_model (cópia).py"
```

### MCMC Configuration

Key parameters in the code:
```python
nwalkers = 40        # Number of MCMC walkers
niter = 10000        # Number of production iterations
burnin = 1500        # Burn-in iterations to discard
initial = [70, 0.3, 0.8, -19.253, 5]  # Initial parameter guess
```

### Workflow

1. **Data Loading**: Reads observational datasets
2. **Model Definition**: Sets up R2-AB gravity equations
3. **Likelihood Computation**: Calculates chi-squared for each dataset
4. **MCMC Sampling**: 
   - Burn-in phase (1500 iterations)
   - Production phase (10000 iterations)
5. **Results Output**:
   - `Results_SNe+CC+RSD_AB.txt` - MCMC chains
   - `log_pro_SNe+CC+RSD_AB.txt` - Log probability values
   - `MCMC_AB_fs8+CC_SNIa.pdf` - Corner plot visualization

## Output Files

After running the code, you'll get:

- **Results_SNe+CC+RSD_AB.txt**: Complete MCMC chain samples (all walkers, all iterations)
- **log_pro_SNe+CC+RSD_AB.txt**: Log probability values for each sample
- **MCMC_AB_fs8+CC_SNIa.pdf**: Corner plot showing parameter distributions and correlations

## Key Functions

### Data Analysis Functions

- `chiCC()` - Chi-squared for Cosmic Chronometer data
- `chiRSD()` - Chi-squared for RSD data  
- `chi2_SN()` - Chi-squared for SNe data
- `chi_tot()` - Combined chi-squared statistic

### Model Functions

- `Hubble()` - ODE system for H(t) evolution
- `solHCC()` - Solves for H(z) at specific redshifts
- `contrast_AB()` - ODE system for density contrast
- `fs8_AB()` - Computes growth rate fs8
- `mag_MG()` - Calculates SNe apparent magnitude

### MCMC Functions

- `lnlike()` - Log-likelihood function
- `lnprior()` - Log-prior probability
- `lnprob()` - Log-posterior probability

## Performance Optimization

The code uses several optimizations:
- **Parallel Processing**: Utilizes `multiprocessing.Pool` for parallel walker evolution
- **Conditional Evaluations**: Checks for numerical stability (e.g., `abs(alpha) < 15`)
- **Efficient Integration**: Uses `LSODA` solver with relative tolerance of 10^-6

## Numerical Considerations

- **Integration Method**: Uses `scipy.integrate.solve_ivp` with LSODA method for stiff ODEs
- **Numerical Stability**: Implements conditional checks to avoid overflow in hyperbolic functions
- **Initial Conditions**: Set at early times (z ≈ 4) assuming GR conditions

## Customization

### Modifying Priors

Edit the `lnprior()` function (lines 400-404):
```python
if 54 < H0 < 76 and 0.1 < O_m0 < 0.5 and ...:
```

### Changing Dataset Weights

Modify the `chi_tot()` function (lines 380-384) to add weights:
```python
return w_cc*chicc + w_rsd*chirsd + w_sn*chisn
```

### Adjusting MCMC Settings

Modify variables around line 415:
- Increase `nwalkers` for better exploration
- Increase `niter` for better convergence
- Adjust `burnin` based on chain diagnostics

## Troubleshooting

### Common Issues

1. **File Not Found**: Update data file paths to your system
2. **Slow Performance**: Reduce `niter` for testing or use more CPU cores
3. **Poor Convergence**: Increase `burnin` or adjust initial parameters
4. **Memory Issues**: Reduce `nwalkers` or `niter`

### Convergence Diagnostics

Check the output to ensure:
- Acceptance fraction is between 0.2-0.5
- Chains have reached steady state
- Corner plots show smooth posterior distributions

## References

This implementation is based on modified gravity cosmology research. Key concepts:

- **R2-AB Model**: A specific form of f(R) gravity with hyperbolic tangent transition
- **Growth Rate**: fs8 quantifies structure formation
- **Cosmic Chronometers**: Determines H(z) from differential ages of galaxies
- **Pantheon+**: Latest Type Ia Supernovae compilation

## Author Notes

- Original Jupyter notebook: `MCMC_Fernanda_Completo.ipynb`
- Google Colab origin: `https://colab.research.google.com/drive/1eAf9QibntI5McERMHSL3-NTmJFQ56VHX`

## License

Please cite appropriate papers if using this code for research purposes.

## Additional Notes

- The code assumes flat universe geometry (k=0.125 regularization parameter)
- Speed of light: c = 1 (natural units used internally)
- Radiation density: Ωr0 = 0 (negligible in late-time universe approximation)
- Delta parameter: 10^-7 (controls transition scale)

---

**Last Updated**: January 2026
