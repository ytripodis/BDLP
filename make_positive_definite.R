# a function that enforces positive definiteness

make_positive_definite <- function(mat, ridge = 1e-1, eigen_threshold = 1e-5) {
  
  if (any(is.nan(mat)) || any(is.infinite(mat))) {
   mat[is.na(mat) | is.infinite(mat)] <- 0
  }
  
  mat <- (mat + t(mat)) / 2  # Ensure symmetry
  
  eig <- eigen(mat, symmetric = TRUE)
  
  # Adjust eigenvalues
  eig$values[eig$values < eigen_threshold] <- eigen_threshold
  
  # Reconstruct matrix
  mat_pd <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  mat_pd <- mat_pd + diag(ridge, nrow(mat_pd))  # Regularization
  mat_pd <- (mat_pd + t(mat_pd)) / 2  # Ensure symmetry again
  
  return(mat_pd)
}

