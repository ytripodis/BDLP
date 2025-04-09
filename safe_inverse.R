# a function that safely inverts a matrix
safe_inverse <- function(mat, ridge = 1e-1, max_attempts = 5, condition_threshold = 1e8) {
  attempt <- 0

  repeat {
    attempt <- attempt + 1
    if (attempt > max_attempts) {
      break
    }

    if (is.null(dim(mat))) {
    }

    # Check for numerical issues
    if (any(is.nan(mat)) || any(is.infinite(mat))) {
      mat[is.na(mat) | is.infinite(mat)] <- 0
    }

    # Compute condition number and monitor stability
    condition_number <- kappa(mat)

    if (condition_number < condition_threshold) {
      # Try direct inversion with regularization
      inv_mat <- try(solve(mat + diag(ridge, nrow(mat))), silent = TRUE)

      if (inherits(inv_mat, "try-error") || any(is.nan(inv_mat)) || any(is.infinite(inv_mat))) {
       } else {
        return(inv_mat)
      }
    } else {
 
      # If matrix is symmetric, use eigenvalue decomposition
      if (isSymmetric(mat)) {
        eigen_decomp <- eigen(mat + diag(ridge, nrow(mat)), symmetric = TRUE)
        inv_mat <- eigen_decomp$vectors %*% diag(1 / (eigen_decomp$values + ridge)) %*% t(eigen_decomp$vectors)
        return(inv_mat)
      } else {
        # General matrix inversion using QR decomposition
        qr_decomp <- qr(mat + diag(ridge, nrow(mat)))
        inv_mat <- try(solve(qr_decomp), silent = TRUE)

        if (!inherits(inv_mat, "try-error")) {
          return(inv_mat)
        } else {
         }
      }
    }

    # If inversion failed, apply stronger regularization
    ridge <- ridge * 1.5  # Instead of multiplying by 10, use a gradual increase
   }

  # If all attempts fail, use pseudoinverse
  return(MASS::ginv(mat))
}
