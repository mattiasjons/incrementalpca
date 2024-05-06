#' Incremental Mean and Variance
#'
#' @param X Matrix like
#' @param last_mean A number.
#' @param last_variance A number.
#' @param last_sample_count A number.
#' @param sample_weight vector like.
#' @returns A list of of updated mean, updated variance and updated sample count
incremental_mean_and_var <- function(X, last_mean, last_variance, last_sample_count, sample_weight=NULL) {
  last_sum <- last_mean * last_sample_count
  X_nan_mask <- is.na(X)

  sum_func <- if (any(X_nan_mask)) sum else sum

  if (!is.null(sample_weight)) {
    new_sum <- matrix(sample_weight, nrow = 1) %*% ifelse(X_nan_mask, 0, X)
    new_sample_count <- sum(sample_weight) - colSums(ifelse(X_nan_mask, 1, 0) * sample_weight)
  } else {
    new_sum <- apply(X, 2, sum_func, na.rm = TRUE)
    n_samples <- nrow(X)
    new_sample_count <- n_samples - colSums(X_nan_mask)
  }

  updated_sample_count <- last_sample_count + new_sample_count
  updated_mean <- (last_sum + new_sum) / updated_sample_count

  if (is.null(last_variance)) {
    updated_variance <- NULL
  } else {
    T <- new_sum / new_sample_count
    temp <- sweep(X, 2, T, "-")
    if (!is.null(sample_weight)) {
      correction <- matrix(sample_weight, nrow = 1) %*% ifelse(X_nan_mask, 0, temp)
      temp <- temp^2
      new_unnormalized_variance <- matrix(sample_weight, nrow = 1) %*% ifelse(X_nan_mask, 0, temp)
    } else {
      correction <- apply(temp, 2, sum_func, na.rm = TRUE)
      temp <- temp^2
      new_unnormalized_variance <- apply(temp, 2, sum_func, na.rm = TRUE)
    }

    new_unnormalized_variance <- new_unnormalized_variance - (correction^2 / new_sample_count)
    last_unnormalized_variance <- last_variance * last_sample_count
    last_over_new_count <- last_sample_count / new_sample_count
    updated_unnormalized_variance <- last_unnormalized_variance + new_unnormalized_variance +
      (last_over_new_count / updated_sample_count) * ((last_sum / last_over_new_count - new_sum)^2)

    zeros <- last_sample_count == 0
    updated_unnormalized_variance[zeros] <- new_unnormalized_variance[zeros]
    updated_variance <- updated_unnormalized_variance / updated_sample_count
  }

  return(list(updated_mean, updated_variance, updated_sample_count))
}



#' Adjust signs of singular vectors from an SVD decomposition
#'
#' This function adjusts the signs of the singular vectors \( U \) and \( V \) matrices from an
#' SVD decomposition to ensure consistent sign orientation. This adjustment is based on maximizing
#' the absolute values of the components in either the rows or columns of the \( U \) or \( V \)
#' matrices, respectively.
#'
#' @param u Numeric matrix representing the \( U \) matrix from an SVD decomposition.
#' @param v Numeric matrix representing the \( V \) matrix from an SVD decomposition.
#' @param u_based_decision Logical scalar indicating if the decision to adjust signs should be
#'        based on \( U \) matrix. If `TRUE`, signs are adjusted based on the columns of \( U \);
#'        if `FALSE`, the adjustment is based on the rows of \( V \).
#'
#' @return A list containing the adjusted \( U \) and \( V \) matrices.
#'
#' @details The function operates by first determining whether the adjustment should be based on
#' the \( U \) or \( V \) matrix. If `u_based_decision` is `TRUE`, the function will transpose \( U \)
#' and calculate the column with the maximum absolute value for each row in the transposed \( U \).
#' Then, it generates sign adjustments based on the signs of these maximal values. Both \( U \) and
#' \( V \) matrices are then adjusted by these signs. If `u_based_decision` is `FALSE`, a similar
#' process is applied but based on rows of \( V \) instead.
#'
#' @examples
#' # SVD decomposition of a matrix
#' A <- matrix(c(1, 2, 3, 2, -1, 4, 3, 6, -4), ncol = 3)
#' svd_A <- svd(A)
#' # Adjusting the signs of U and V
#' adjusted_svd <- svd_flip(svd_A$u, svd_A$v)
#' adjusted_u <- adjusted_svd[[1]]
#' adjusted_v <- adjusted_svd[[2]]
#'
#' @seealso \code{\link[base]{svd}}
#'
svd_flip <- function(u, v, u_based_decision=TRUE) {
  if (u_based_decision) {
    # Transpose u for easier manipulation of columns
    u_t <- t(u)

    # Find the indices of the maximum absolute values in each column of u
    max_abs_u_cols <- apply(u_t, 1, function(col) which.max(abs(col)))

    # Generate signs based on the values at these indices
    signs <- sign(sapply(max_abs_u_cols, function(idx, col) col[idx], col = u_t))

    # Adjust u and v by these signs
    u_adjusted <- u * signs
    v_adjusted <- t(t(v) * signs) # Need to transpose v twice for correct dimensionality
  } else {
    # Find the indices of the maximum absolute values in each row of v
    max_abs_v_rows <- apply(v, 1, function(row) which.max(abs(row)))

    # Generate signs based on the values at these indices
    signs <- sign(sapply(max_abs_v_rows, function(idx, row) row[idx], row = v))

    # Adjust u and v by these signs
    u_adjusted <- u * signs
    v_adjusted <- t(t(v) * signs) # Apply signs to each column of v
  }

  return(list(u_adjusted, v_adjusted))
}


