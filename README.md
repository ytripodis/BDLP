# Bayesian Dynamic Latent Process Modeling in R

## Overview
This repository contains R scripts implementing Bayesian Dynamic Latent Process modeling, including data simulation, Kalman filtering, Kalman smoothing, and Gibbs sampling, tailored specifically for modeling longitudinal multivariate data with latent dynamic processes.

## Repository Contents
- `DataGenerate.R`: Simulate data for Bayesian dynamic latent process models.
- `KalmanFilter.R`: Kalman filter algorithm for latent state estimation.
- `KalmanSmoother.R`: Kalman smoother for enhanced state estimation.
- `GibbsSampler.R`: Gibbs sampler for Bayesian inference.
- `GibbsSampleRun.R`: Example execution of the Gibbs sampler.
- `make_positive_definite.R`: Utility for numerical stability.
- `safe_inverse.R`: Robust matrix inversion tool.

## Installation

Install required R packages:

install.packages(c("MASS", "MCMCpack", "Matrix"))

## Usage

### Generate Simulated Data
 
Rscript DataGenerate.R

### Run Gibbs Sampling
 
Rscript GibbsSampleRun.R

### Expected Output
After running the Gibbs sampler, you will get:

A saved RData object: GibbsSamplerResults.RData

Posterior means and 95% credible intervals for:

Transition matrix A

Coefficients beta

Loadings P

Latent process noise Sigma_eta

Measurement noise Sigma_eps

## References
Shumway & Stoffer (2017). Time Series Analysis and its Applications.

Durbin & Koopman (2012). Time Series Analysis by State Space Methods.

Taddé et al. (2020). Biometrics, "Dynamic modeling of multivariate dimensions and their temporal relationships using latent processes."

## License
This project is licensed under the MIT License.

## Citation
If you use this code, please cite it as:

Yorghos Tripodis. (2025). Bayesian Dynamic Latent Process Modeling. GitHub repository: https://github.com/ytripodis/BDLP

## Contact
Yorghos Tripodis
Boston University School of Public Health
📧 yorghos@bu.edu
🔗 GitHub Profile: ytripodis




