# Gibbs Sample Run

# Load necessary packages
library(MASS)
library(MCMCpack)
library(Matrix)  # For nearPD function

source("C:/Users/yorghos/Dropbox/Grants/SSM_grant/Bordeaux/R codes/GibbsSampler_v2/DataGenerate.R")

# Define directory where all functions are stored
path <- file.path("C:/Users/yorghos/Dropbox/Grants/SSM_grant/Bordeaux/R codes/GibbsSampler_v2")
# Set working directory
setwd("C:\\Users\\yorghos\\Dropbox\\Grants\\SSM_grant\\Bordeaux\\Data")

# Load the generated data from CSV files
Y_data <- read.csv("Y_data.csv", header = FALSE)
X_data <- read.csv("X_data.csv", header = FALSE)
Lambda_data <- read.csv("Lambda_data.csv", header = FALSE)

# Reshape data into arrays
N <- 100  # Number of individuals
J <- 10   # Number of time points
D <- 2    # Number of latent processes
K <- 4    # Number of observed markers

delta <- 1  # Time interval
iterations <- 1000  # Number of MCMC iterations
burn_in <- 500  # Burn-in period

# Convert to arrays
Y <- array(as.numeric(as.matrix(Y_data)), dim = c(N, J, K))
X <- array(as.numeric(as.matrix(X_data)), dim = c(N, J, 2))
Lambda <- array(as.numeric(as.matrix(Lambda_data)), dim = c(N, J, D))

# Define C_A and C_P (Structural Constraints)
C_A <- matrix(c(0, 1, 
                0, 0), D, D, byrow = TRUE)  

C_P <- matrix(c(1, 0,
                1, 0,
                0, 1,
                0, 1), K, D, byrow = TRUE)  # Estimate P[1,1], P[2,2], and P[3,2] only

# Initialize Parameters for the Gibbs Sampler

# More stable initialization of A
A <- matrix(c(0.0, 0.1,
              0.0, 0.0), D, D, byrow = TRUE) # Reduced values for stability


# More stable initialization of P
P <- matrix(c(0.6, 0.0,
              0.64,0.0,
              0.0, 0.8,
              0.0, 0.64), K, D, byrow = TRUE)

# Initialize noise covariance matrices with smaller values
Sigma_eta <- diag(1.0, D)  # Reduced variance for latent process noise
Sigma_eps <- diag(1.0, K)  # Reduced variance for measurement noise

# Initialize Beta
beta <- matrix(0, 2, 1)

# Priors

A_prior_mean <- rep(0, length(which(as.vector(C_A)==1)))  # Mean vector of priors for the elements of A
A_prior_cov <- diag(5, length(which(as.vector(C_A)==1)))  # Prior covariance matrix, smaller than before for better control

beta_prior_mean <- rep(0, 2)
beta_prior_cov <- diag(2) * 0.5  # Increased regularization

P_prior_mean <- matrix(0, K, D)
P_prior_cov <- diag(sum(C_P)) * 0.5 

Sigma_eta_prior_scale <- 2 # LKJ prior parameter (choose ν > 0; ν=1 is uniform)

# Load the safe_inverse and make_positive_definite functions
source("C:/Users/yorghos/Dropbox/Grants/SSM_grant/Bordeaux/R codes/GibbsSampler_v2/safe_inverse.R")  
source("C:/Users/yorghos/Dropbox/Grants/SSM_grant/Bordeaux/R codes/GibbsSampler_v2/make_positive_definite.R")  
source("C:/Users/yorghos/Dropbox/Grants/SSM_grant/Bordeaux/R codes/GibbsSampler_v2/KalmanFilter.R")
source("C:/Users/yorghos/Dropbox/Grants/SSM_grant/Bordeaux/R codes/GibbsSampler_v2/KalmanSmoother.R")
source("C:/Users/yorghos/Dropbox/Grants/SSM_grant/Bordeaux/R codes/GibbsSampler_v2/GibbsSampler.R")

# Run the Gibbs Sampler
results <- tryCatch({
  start_time <- Sys.time()

  # Run the Gibbs sampler
  results <- gibbs_sampler(Y = Y, X = X, A = A, P = P, C_A = C_A, C_P = C_P,
                           beta = beta, Sigma_eta = Sigma_eta, Sigma_eps = Sigma_eps,
                           A_prior_mean=A_prior_mean, A_prior_cov=A_prior_cov, beta_prior_mean=beta_prior_mean, 
                           P_prior_mean=P_prior_mean, P_prior_cov = P_prior_cov, Sigma_eta_prior_scale = Sigma_eta_prior_scale,
                           iterations = iterations, burn_in = burn_in, delta = delta)
  
  
  if (is.null(results)) stop("Gibbs sampler returned NULL unexpectedly.")
  if (!is.list(results)) stop("Gibbs sampler returned an object that is not a list.")
  
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  results$elapsed_time <- elapsed_time
  
  cat("Gibbs sampling completed successfully in", round(elapsed_time, 2), "seconds.\n")
  results
  
}, error = function(e) {
  cat("An error occurred during Gibbs sampling:\n", e$message, "\n")
  traceback()
  NULL
})

# Save results if sampling was successful
if (!is.null(results)) {
  save(results, file = "GibbsSamplerResults.RData")
  cat("Gibbs sampling completed successfully.\n")
} else {
  cat("Gibbs sampling failed. Check the error message above.\n")
}

# Calculate posterior means for key parameters:
A_mean <- apply(results$A_samples, c(1,2), mean)
beta_mean <- colMeans(results$beta_samples)
P_mean <- apply(results$P_samples, c(1,2), mean)
Sigma_eta_mean <- apply(results$Sigma_eta_samples, c(1,2), mean)
Sigma_eps_mean <- apply(results$Sigma_eps_samples, c(1,2), mean)

#Compute 95% credible intervals to quantify uncertainty:
A_CI <- apply(results$A_samples, c(1,2), quantile, probs = c(0.025, 0.975))
beta_CI <- apply(results$beta_samples, 2, quantile, probs = c(0.025, 0.975))

