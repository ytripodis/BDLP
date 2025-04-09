kalman_filter <- function(Y, X, A, P, beta, Sigma_eta, Sigma_eps) {
  # cat("Kalman Filter: Function is initiated succesfully.\n")
  N <- dim(Y)[1]
  J <- dim(Y)[2]
  D <- dim(A)[1]
  
  ridge <- 1e-1  # Small regularization constant
  
  # Storage for filtered means and covariances
  Lambda_filtered <- array(0, dim = c(N, J, D))
  V_filtered <- array(0, dim = c(N, J, D, D))
  
  # Initial state
  Lambda_filtered[, 1, ] <- 0  # Assuming zero initialization
  V_filtered[, 1, , ] <- diag(D)  # Identity matrix as initial covariance

  for (i in 1:N) {

    for (j in 2:J) {
       # Prediction Step
      Lambda_pred <- (diag(D) + A) %*% matrix(Lambda_filtered[i, j - 1, ], ncol = 1)
      covariate_effect_value <- as.numeric(t(X[i, j, ]) %*% beta)
      covariate_effect_matrix <- matrix(rep(covariate_effect_value, D), nrow = D, ncol = 1)
      
      Lambda_pred <- Lambda_pred + covariate_effect_matrix
      
      # Check dimension compatibility
      if (ncol(P) != nrow(Lambda_pred)) {
        cat("ERROR: Mismatch in dimensions: P is ", dim(P), " and Lambda_pred is ", dim(Lambda_pred), "\n")
        stop("Mismatch in dimensions.")
      }
      
      # Observation Prediction
      Y_pred <- P %*% Lambda_pred
      
      # Predictive Covariance
      V_pred <- (diag(D) + A) %*% V_filtered[i, j - 1, , ] %*% t(diag(D) + A) + Sigma_eta
      V_pred <- make_positive_definite(V_pred, ridge = ridge)
      
      if (is.vector(Lambda_pred)) {
        Lambda_pred <- matrix(Lambda_pred, ncol = 1)
      }

      if (ncol(P) != nrow(Lambda_pred)) {
        cat("ERROR: Mismatch in dimensions: P is ", dim(P), " and Lambda_pred is ", dim(Lambda_pred), "\n")
        stop("Mismatch in dimensions.")
      }
      

      Y_pred <- P %*% Lambda_pred
      S <- make_positive_definite(P %*% V_pred %*% t(P) + Sigma_eps, ridge = ridge)
      
      K <- V_pred %*% t(P) %*% safe_inverse(S,ridge = ridge)
      V_filtered[i, j, , ] <- V_pred-K %*% P%*% V_pred
      Lambda_filtered[i, j, ] <- Lambda_pred + K %*% (Y[i, j, ] - Y_pred)

    }
  }
  if (any(is.nan(Lambda_filtered)) || any(is.infinite(Lambda_filtered))) {
    warning("Kalman Filter: Lambda_filtered contains NA or infinite values. Regularizing.")
    Lambda_filtered[is.na(Lambda_filtered) | is.infinite(Lambda_filtered)] <- 0
  }
  if (any(is.nan(V_filtered)) || any(is.infinite(V_filtered))) {
    warning("Kalman Filter: V_filtered contains NA or infinite values. Regularizing.")
    V_filtered[is.na(V_filtered) | is.infinite(V_filtered)] <- 0
  }
  if (is.null(Lambda_filtered)) {
    stop("Error: Lambda_filtered is not properly calculated. Check matrix dimensions and operations.\n")
  } 
  if (is.null(V_filtered)) {
    stop("Error: V_filtered is not properly calculated. Check matrix dimensions and operations.\n")
  }

  return(list(Lambda_filtered = Lambda_filtered, V_filtered = V_filtered))

  }
