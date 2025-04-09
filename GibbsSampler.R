library(MASS)
library(MCMCpack)
library(Matrix)

gibbs_sampler <- function(Y, X, A, P, C_A, C_P, beta, Sigma_eta, Sigma_eps, iterations, burn_in, delta,
                          A_prior_mean  = as.vector(diag(D)), A_prior_cov = diag(D^2) * 10,
                          beta_prior_mean = rep(0, 2), beta_prior_cov = diag(2) * 10,
                          P_prior_mean = matrix(0, K, D), 
                          P_prior_cov = as.matrix(bdiag(replicate(K, diag(D) * 10, simplify = FALSE))),
                          Sigma_eta_prior_scale = 2,
                          Sigma_eps_prior_df = K + 2, Sigma_eps_prior_scale = diag(K)) {
  N <- dim(Y)[1]
  J <- dim(Y)[2]
  D <- dim(A)[1]
  K <- dim(Y)[3]
  
  ridge <- 1e-1
  
  # Storage for samples
  num_samples <- iterations - burn_in
  if (num_samples <= 0) {
    stop("Number of samples to store is zero or negative. Check 'iterations' and 'burn_in'.")
  }
  
  A_samples <- array(NA, dim = c(D, D, num_samples))
  beta_samples <- matrix(NA, nrow = num_samples, ncol = length(beta))
  P_samples <- array(NA, dim = c(K, D, num_samples))
  Sigma_eta_samples <- array(NA, dim = c(D, D, num_samples))
  Sigma_eps_samples <- array(NA, dim = c(K, K, num_samples))
  

  for (iter in 1:iterations) {
    # cat("Iteration ", iter, "\n")
    
    # Unified Step: Joint Sampling of A and beta using transformed Z from Y via P
    tryCatch({
      # Step 1: Construct Z_ij and Delta Z_ij
      Z_array <- array(0, dim = c(N, J, D))
      PtP_inv_Pt <- safe_inverse(t(P) %*% P) %*% t(P)
      
      for (i in 1:N) {
        for (j in 1:J) {
          Z_array[i, j, ] <- as.numeric(PtP_inv_Pt %*% Y[i, j, ])
        }
      }
      
      # Flattened matrices for regression
      Z_mat <- matrix(0, N * (J - 1), D)
      dZ_mat <- matrix(0, N * (J - 1), D)
      X_flat <- matrix(0, N * (J - 1), length(beta))
      
      row_idx <- 1
      for (i in 1:N) {
        for (j in 2:J) {
          Z_mat[row_idx, ] <- Z_array[i, j - 1, ]
          dZ_mat[row_idx, ] <- Z_array[i, j, ] - Z_array[i, j - 1, ]
          X_flat[row_idx, ] <- X[i, j, ]
          row_idx <- row_idx + 1
        }
      }
      
      # Step 2: Compute Omega = δ * Sigma_eta + 2 * Sigma_eps + A Σ_eps A^T
      Omega <- delta * Sigma_eta + 2 * PtP_inv_Pt%*%Sigma_eps%*%t(PtP_inv_Pt) + A %*% PtP_inv_Pt%*%Sigma_eps%*%t(PtP_inv_Pt) %*% t(A)
      Omega <- make_positive_definite(Omega)
      
      # Step 3: Construct regression design matrix W
      C_vec <- as.vector(C_A)
      C_id <- which(C_vec == 1)
      K_A <- length(C_id)
      
      W_mat <- matrix(0, N * (J - 1) * D, K_A + length(beta))
      dZ_vec <- as.vector(t(dZ_mat))
      
      row_ptr <- 1
      for (i in 1:nrow(Z_mat)) {
        Z_i <- Z_mat[i, ]
        W_A <- kronecker(t(Z_i), diag(D))[, C_id, drop = FALSE]
        W_beta <- X_flat[i, ]
        W_row <- cbind(W_A, matrix(rep(W_beta, each = D), nrow = D, byrow = TRUE))
        W_mat[row_ptr:(row_ptr + D - 1), ] <- W_row
        row_ptr <- row_ptr + D
      }
      
      # Step 4: Posterior update for theta = [vec(A), beta]
      Omega_inv <- safe_inverse(Omega)
      Sigma_inv_big <- kronecker(diag(N * (J - 1)), Omega_inv)

      Wt_Oinv_W <- t(W_mat) %*% Sigma_inv_big %*% W_mat
      Wt_Oinv_dZ <- t(W_mat) %*% Sigma_inv_big %*% dZ_vec
      
      
      # Priors
      mu_theta <- c(A_prior_mean, beta_prior_mean)
      V_theta_inv <- diag(c(1 / diag(A_prior_cov), 1 / diag(beta_prior_cov)))
      
      V_star <- safe_inverse(Wt_Oinv_W + V_theta_inv)
      theta_hat <- V_star %*% (Wt_Oinv_dZ + V_theta_inv %*% mu_theta)
      
      # Sample theta
      theta_sample <- mvrnorm(1, theta_hat, make_positive_definite(V_star))
      
      # Step 5: Update A
      A_vec <- rep(0, D * D)
      A_vec[C_id] <- theta_sample[1:K_A]
      A <- matrix(A_vec, D, D)
      
      # Step 6: Update beta
      beta <- theta_sample[(K_A + 1):length(theta_sample)]
      
    }, error = function(e) {
      cat("Error in Joint Sampling of A and beta from transformed Z:", e$message, "\n")
    })
    
       # Step 2: Calculation of Kalman Filter and Smoother Draw Lambda (Backward Sampling)
    #cat("Step 2","\n")
    tryCatch({
      
      if (!all(dim(P) == c(K, D))) stop("P matrix has incorrect dimensions")
      if (!all(dim(A) == c(D, D))) stop("A matrix has incorrect dimensions")
      
      # Run Kalman Filter with updated A and beta
      kf_result <- kalman_filter(Y, X, A, P, beta, Sigma_eta, Sigma_eps)
      
      # Run Kalman Smoother with filtered results
      ks_result <- kalman_smoother(kf_result$Lambda_filtered, kf_result$V_filtered, A, Sigma_eta)

      Lambda_smoothed <- ks_result$Lambda_smoothed
      V_smoothed <- ks_result$V_smoothed
 
      if (is.null(Lambda_smoothed) || is.null(V_smoothed)) {
        stop("Kalman Smoother failed to produce valid results.")
      }
    }, error = function(e) {
      cat("Error in Step 2:", e$message, "\n")
    })
    
    
    # Step 2: Draw Lambda (Backward Sampling)
    #cat("Step 2","\n")
    tryCatch({
      Lambda <- array(rnorm(N * J * D), dim = c(N, J, D))
      
      #cat("Initial Lambda (before update):\n")
      #print(Lambda)
      
      for (i in 1:N) {
        Lambda[i, J, ] <- mvrnorm(1, Lambda_smoothed[i, J, ], make_positive_definite(V_smoothed[i, J, , ],ridge = ridge))
        for (j in (J-1):1) {
          G <- V_smoothed[i, j, , ] %*% t(diag(D) + A) %*% 
            safe_inverse((diag(D) + A) %*% V_smoothed[i, j, , ] %*% t(diag(D) + A) + Sigma_eta)
          Lambda_mean_post <- Lambda_smoothed[i, j, ] + G %*% (Lambda[i, j + 1, ] - (diag(D) + A) %*% Lambda_smoothed[i, j, ])
          V_post <- V_smoothed[i, j, , ] - G %*% ((diag(D) + A) %*% V_smoothed[i, j, , ] %*% t(diag(D) + A) + Sigma_eta - V_smoothed[i, j + 1, , ]) %*% t(G)
          Lambda[i, j, ] <- mvrnorm(1, Lambda_mean_post, make_positive_definite(V_post))
        }
      }
      
      # Regularize Lambda
      Lambda[is.na(Lambda) | is.infinite(Lambda)] <- 0
      Lambda[abs(Lambda) > 1000] <- sign(Lambda[abs(Lambda) > 1000]) * 1000
      Lambda[abs(Lambda) < 0.0001] <- sign(Lambda[abs(Lambda) < 0.0001]) * 0.0001
      
    }, error = function(e) {
      cat("Error in Step 2:", e$message, "\n")
    })
    
    if (is.null(Lambda) || is.null(V_smoothed)) {
      stop("Kalman Smoother failed to produce valid results.")
    }
     
     # Step 3: Robust sampling of P

    tryCatch({
      for (k in 1:K) {
        # Identify relevant indices to estimate for factor-loading row k
        relevant_indices <- which(C_P[k, ] != 0)
        
        # Extract relevant submatrices based on identified indices
        Lambda_relevant <- Lambda[, , relevant_indices, drop = FALSE]
        P_prior_cov_sub <- P_prior_cov[relevant_indices, relevant_indices, drop = FALSE]
        
        # Reshape Lambda matrix for regression
        Lambda_mat <- matrix(Lambda_relevant, ncol = length(relevant_indices))
        Y_k_vec <- as.vector(Y[, , k])
        
        # Compute posterior covariance
        Sigma_eps_k_inv <- 1 / Sigma_eps[k, k]
        P_cov_k_inv <- solve(t(Lambda_mat) %*% Lambda_mat * Sigma_eps_k_inv + safe_inverse(P_prior_cov_sub))
        
        # Compute posterior mean
        P_mean_k <- P_cov_k_inv %*% (t(Lambda_mat) %*% Y_k_vec * Sigma_eps_k_inv + solve(P_prior_cov_sub) %*% P_prior_mean[k, relevant_indices])
        
        # Sample from posterior
        P_sample_k <- mvrnorm(1, mu = P_mean_k, Sigma = P_cov_k_inv)
        
        # Update row k of P
        P[k, relevant_indices] <- P_sample_k
      }
      
      # Enforce identifiability constraints
      # Normalize columns of P to have unit Euclidean norm
      for (d in 1:D) {
        P[, d] <- P[, d] / sqrt(sum(P[, d]^2))
      }
      
      # Impose positivity constraints
      P[P < 0 & C_P == 1] <- abs(P[P < 0 & C_P == 1])
      }, error = function(e) {
        cat("dimensions of P_prior_cov_sub ",dim(P_prior_cov_sub),"\n")
        cat("Lambda_mat matrix=: ",Lambda_mat,"\n")
        cat("Error in Iteration ", iter, "\n")
      cat("Error in Step 3:", e$message, "\n")
    })
    


    # Step 4: Estimation of Measurement Error Covariance Matrix Sigma_eps
    
    # Required inputs:
    # Y: observed data array (N x J x K)
    # Lambda: latent process array (N x J x D)
    # P: current estimate of factor-loading matrix (K x D)
    # Sigma_eps_prior_df: prior degrees of freedom for inverse Gamma distribution
    # Sigma_eps_prior_scale: prior scale parameter matrix for inverse Gamma distribution (K x K)
    tryCatch({ 
    Sigma_eps_new <- matrix(0, K, K)
    
    for (k in 1:K) {
      residuals_k <- numeric(N * J)
      
      # Calculate residuals for observed marker k
      for (i in 1:N) {
        for (j in 1:J) {
          predicted_value <- sum(P[k, ] * Lambda[i, j, ])
          residuals_k[(i - 1) * J + j] <- Y[i, j, k] - predicted_value
        }
      }
      
      # Posterior parameters for inverse gamma distribution
      posterior_shape <- (N * J + Sigma_eps_prior_df) / 2
      posterior_scale <- (sum(residuals_k^2) + Sigma_eps_prior_scale[k, k]) / 2
      
      # Sample from posterior inverse gamma distribution
      Sigma_eps_new[k, k] <- 1 / rgamma(1, shape = posterior_shape, rate = posterior_scale)
      
      # Regularize Sigma_eps_new
      if (is.na(Sigma_eps_new[k, k])|is.infinite(Sigma_eps_new[k, k])) {
        Sigma_eps_new[k, k] <- 0.0001
      }
      if (Sigma_eps_new[k, k]>1000) {
        Sigma_eps_new[k, k] <- 1000
      }
      if (Sigma_eps_new[k, k]<0.0001) {
        Sigma_eps_new[k, k] <- 0.0001
      }
    }
    Sigma_eps <- Sigma_eps_new
    }, error = function(e) {
      cat("Error in Step 4:", e$message, "\n")
    })
  
    # Step 5: Estimation of Correlation Matrix Sigma_eta
    # Purpose: Estimate the correlation matrix Sigma_eta of the latent
    #          process increments within a Bayesian Gibbs sampling framework.
    # Constraints: Sigma_eta must have ones on diagonal, symmetric, and positive definite.
    # Methodology: Posterior sampling using LKJ prior (Lewandowski et al., 2009)
    #              and Fisher Z-transform for stable numerical sampling.
    # Inputs:
    #   - Lambda: Current latent state estimates from Kalman smoother (N x J x D).
    #   - A: Current temporal influence (transition) matrix (D x D).
    #   - nu: LKJ prior hyperparameter controlling correlation strength (nu=1 uniform).
    # Outputs:
    #   - Sigma_eta: Updated correlation matrix estimate (D x D).
    # Dependencies: Requires function make_positive_definite() to ensure numerical stability.
    #--
    tryCatch({
      # Compute residuals for Sigma_eta
      eta_resid <- do.call(rbind, lapply(1:N, function(i) {
        diff_Lambda <- diff(Lambda[i,,])  # Difference between consecutive time points
        Lambda_A_prod <- Lambda[i,-J,] %*% t(A)
        diff_Lambda - Lambda_A_prod
      }))
      
      # Check for numerical issues
      if (anyNA(eta_resid) || any(is.nan(eta_resid)) || any(is.infinite(eta_resid))) {
        stop("eta_resid contains NA or NaN values. Check your Lambda calculations.")
      }
      
      # Compute sample correlation matrix
      Sigma_eta_cor <- cor(eta_resid)
      
      # Regularize initial correlation matrix if necessary
      Sigma_eta_cor <- make_positive_definite(Sigma_eta_cor, ridge = 1e-3)
      
      # Ensure symmetry and diagonals exactly equal to 1
      Sigma_eta_cor <- (Sigma_eta_cor + t(Sigma_eta_cor)) / 2
      diag(Sigma_eta_cor) <- 1
      
      # Extract lower-triangular elements for Fisher Z-transform
      Z_eta <- atanh(Sigma_eta_cor[lower.tri(Sigma_eta_cor)])
      
      # LKJ prior parameter (choose ν > 0; ν=1 is uniform)
      nu <- Sigma_eta_prior_scale
      
      # Posterior variance of Fisher Z-transform (depends on ν and data size)
      posterior_sd <- sqrt(1 / (N*(J - 1) - 3 + 2*(nu - 1)))
      
      # Sample posterior Z values with LKJ prior approximation
      Z_eta_sampled <- rnorm(length(Z_eta), mean = Z_eta, sd = posterior_sd)
      
      # Transform back to correlation scale
      Sigma_eta_new_cor <- diag(D)
      Sigma_eta_new_cor[lower.tri(Sigma_eta_new_cor)] <- tanh(Z_eta_sampled)
      Sigma_eta_new_cor[upper.tri(Sigma_eta_new_cor)] <- t(Sigma_eta_new_cor)[upper.tri(Sigma_eta_new_cor)]
      
      # Ensure positive definiteness and symmetry again
      Sigma_eta <- make_positive_definite(Sigma_eta_new_cor, ridge = 1e-3)
      
      # Explicitly set diagonals to exactly 1
      diag(Sigma_eta) <- 1
      
    }, error = function(e) {
      cat("Error in Sigma_eta sampling:", e$message, "\n")
    })
    
    
    # Save samples after burn-in
    if (iter > burn_in) {
      idx <- iter - burn_in
      A_samples[, , idx] <- A
      beta_samples[idx, ] <- beta
      P_samples[, , idx] <- P
      Sigma_eta_samples[, , idx] <- Sigma_eta
      Sigma_eps_samples[, , idx] <- Sigma_eps
    }
  } 
  # cat("Ready to return results. Samples stored successfully.\n")
  if (iter == iterations) {
    return(list(
      A_samples = A_samples,
      beta_samples = beta_samples,
      P_samples = P_samples,
      Sigma_eta_samples = Sigma_eta_samples,
      Sigma_eps_samples = Sigma_eps_samples
    ))
  }
}

