# Load necessary packages
library(MASS)
# Set working directory
setwd("C:/Users/yorghos/Dropbox/Grants/SSM_grant/Bordeaux/Data")
# Set seed for reproducibility
set.seed(123)

# Parameters
N <- 100       # Number of individuals
J <- 10        # Number of time points
D <- 2         # Dimension of the latent process
K <- 4         # Dimension of the observed data

delta <- 1  # time between observations

# Define variance-covariance matrices
Sigma_eta <- matrix(c(1.00, 0.05,
                      0.05,1.00),D,D, byrow = TRUE)  # Variance-covariance of the Gaussian process
Sigma_eps <- diag(c(0.1, 0.2,0.3,0.4))  # Variance-covariance of the measurement error

# Simulated parameters
A <- matrix(c(0.0, 0.1,
              0.0, 0.0), D, D, byrow = TRUE)
beta <- matrix(c(0, -0.3), 2, 1)
P <- matrix(c(0.60, 0.00,
              0.64, 0.00,
              0.00, 0.80, 
              0.00, 0.36), K, D, byrow = TRUE)

# Generate covariates X
#X <- array(rnorm(N * J * 2), dim = c(N, J, 2))
X <- array(0, dim = c(N, J, 2))
X[, , 1] <- 1
X[, , 2] <- rnorm(N * J)

# Storage for latent processes and observed data
Lambda <- array(0, dim = c(N, J, D))
Y <- array(0, dim = c(N, J, K))

# Generate latent processes and observed data
for (i in 1:N) {
  # Initialize the latent process at time t=0
  Lambda[i, 1, ] <- rnorm(D)
  for (j in 2:J) {
    # Covariate effect
    covariate_effect <- as.numeric(t(X[i, j, ]) %*% beta)
    
    # Latent process update using the Markov process equation
    eta <- mvrnorm(1, mu = rep(0, D), Sigma = delta * Sigma_eta)
    Lambda[i, j, ] <- covariate_effect + (A + diag(D)) %*% Lambda[i, j - 1, ] + eta
    
    # Measurement model
    epsilon <- mvrnorm(1, mu = rep(0, K), Sigma = Sigma_eps)
    Y[i, j, ] <- P %*% Lambda[i, j, ] + epsilon
  }
}

# Save generated data as CSV files for compatibility
write.table(as.data.frame(as.matrix(Y)), "Y_data.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(as.data.frame(as.matrix(X)), "X_data.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(as.data.frame(as.matrix(Lambda)), "Lambda_data.csv", sep = ",", row.names = FALSE, col.names = FALSE)

cat("Data generation completed successfully.\n")
