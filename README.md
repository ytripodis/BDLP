# Bayesian Dynamic Latent Process Modeling in R

## Overview
This repository contains R scripts implementing Bayesian Dynamic Latent Process modeling, including data simulation, Kalman filtering, Kalman smoothing, and Gibbs sampling, tailored specifically for modeling longitudinal multivariate data with latent dynamic processes.

## Repository Contents
- `DataGenerate.R`: Simulate data for Bayesian dynamic latent process models.
- `KalmanFilter.R`: Kalman filter algorithm for latent state estimation.
- `KalmanSmoother.R`: Kalman smoothing for enhanced state estimation.
- `GibbsSampler.R`: Gibbs sampler for Bayesian inference.
- `GibbsSampleRun.R`: Example execution of the Gibbs sampler.
- `make_positive_definite.R`: Utility for numerical stability.
- `safe_inverse.R`: Robust matrix inversion tool.

## Installation
```R
install.packages(c("MASS", "MCMCpack", "Matrix"))
# BDLP
Bayesian Dynamic Latent Processes
# Generate simulated data
Rscript DataGenerate.R

# Run Gibbs sampling
Rscript GibbsSampleRun.R
Expected Output
After execution, posterior samples and diagnostic summaries are saved, including parameter estimates, credible intervals, and convergence diagnostics.

References
Shumway & Stoffer (2017): Time Series Analysis and its Applications.

Durbin & Koopman (2012): Time Series Analysis by State Space Methods.

Taddé et al. (2020): Biometrics, "Dynamic modeling of multivariate dimensions and their temporal relationships using latent processes."

License
This repository is available under the MIT License.

Citation
Please cite this repository as:

Yorghos Tripodis, (2025). Bayesian Dynamic Latent Process Modeling. GitHub repository, https://github.com/[your-github-user]/[ytripodis].
Contact
[Yorghos Tripodis]

[Boston University School of Public Health]

Email: yorghos@bu.edu

GitHub: your-github-profile



---

This structured README and MIT license will maximize usability, accessibility, and academic impact of your repository.






