# Kalman Smoother
kalman_smoother <- function(Lambda_filtered, V_filtered, A, Sigma_eta) {

  N <- dim(Lambda_filtered)[1]
  J <- dim(Lambda_filtered)[2]
  D <- dim(A)[1]
  
  ridge <- 1e-1
  max_ridge <- 10
  eigen_threshold <- 1e-6  # Minimum eigenvalue threshold for numerical stability
  
  # Storage for smoothed means and covariances
  Lambda_smoothed <- Lambda_filtered
  V_smoothed <- V_filtered
  
  for (i in 1:N) {
    for (j in (J - 1):1) {  # Backward recursion
      
      # Smoothed estimates
      V_pred <- make_positive_definite((diag(D) + A) %*% V_filtered[i, j, , ] %*% t(diag(D) + A) + Sigma_eta, ridge = ridge)

      G <- V_filtered[i, j, , ] %*% t(diag(D) + A) %*% safe_inverse(V_pred)
      
        # Update smoothed means and covariances
      predicted_value <- (diag(D) + A) %*% Lambda_filtered[i, j, ]
      smoothing_residual <- matrix(as.numeric(Lambda_smoothed[i, j + 1, ] - predicted_value), ncol = 1)
 
      smoothing_update <-G %*%smoothing_residual

      # Update Lambda_smoothed calculation
      Lambda_smoothed[i, j, ] <- Lambda_filtered[i, j, ] + smoothing_update
      
      V_smoothed[i, j, , ] <- V_filtered[i, j, , ] + G %*% (V_smoothed[i, j + 1, , ] - V_pred) %*% t(G)
      
      # Ensure symmetry in V_smoothed
      V_smoothed[i, j, , ] <- make_positive_definite(V_smoothed[i, j, , ], ridge = ridge)
    }
  }
  
  if (any(is.nan(Lambda_smoothed)) || any(is.infinite(Lambda_smoothed))) {
    warning("Kalman Smoother: Lambda_smoothed contains NA or infinite values. Regularizing.")
    Lambda_filtered[is.na(Lambda_smoothed) | is.infinite(Lambda_smoothed)] <- 0
  }
  if (any(is.nan(V_smoothed)) || any(is.infinite(V_smoothed))) {
    warning("Kalman Smoother: V_smoothed contains NA or infinite values. Regularizing.")
    V_filtered[is.na(V_smoothed) | is.infinite(V_smoothed)] <- 0
  }
  if (is.null(Lambda_smoothed)) stop("Lambda_smoothed is NULL")
  if (is.null(V_smoothed)) stop("V_smoothed is NULL")
  
  return(list(Lambda_smoothed = Lambda_smoothed, V_smoothed = V_smoothed))
}
