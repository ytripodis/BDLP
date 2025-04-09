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

```r
install.packages(c("MASS", "MCMCpack", "Matrix"))





