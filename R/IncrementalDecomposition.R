library(R6)
#' @title Incremental Decomposition Class
#' @description R6 class for performing and storing the results of an incremental
#'   principal component analysis (PCA) or other similar decompositions.
#'   Heavily influenced by (and sometimes and sometimes a port of)
#'   \href{https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.IncrementalPCA.html}{Scikit Learn's Incremental PCA implementation}
#'   @examples
#'   d = 20
#'   rho <- runif(d)
#'   sigma <- d:1
#'   df <- mvrnorm(100, rep(0, d), t(rho%*%diag(sigma))%*%t(rho))
#'
#'   inc.pca <- IncrementalDecomposition$new(df[1:20,], 10)
#'   image(inc.pca$get_covariance(), ylim=c(1, 0))
#'   image(inc.pca$get_precision(), ylim=c(1, 0))
#'
#'   inc.pca$partial_fit(df[21:50,], lambda = 0.98)
#'   image(inc.pca$get_precision(), ylim=c(1, 0))
#'

#' @export
IncrementalDecomposition <- R6Class("IncrementalDecomposition",
                         public = list(
                           #' @field n_components Integer.
                           #'   The number of principal components to retain in the model.
                           n_components = NULL,

                           #' @field components_ Matrix.
                           #'   The matrix of principal components. Each column represents a principal component.
                           components_ = NULL,

                           #' @field noise_variance_ Numeric.
                           #'   An estimate of the noise variance in the data.
                           noise_variance_ = NULL,

                           #' @field explained_variance_ Numeric vector.
                           #'   The amount of variance explained by each of the principal components.
                           explained_variance_ = NULL,

                           #' @field n_samples_seen_ Integer.
                           #'   The cumulative number of samples processed by the model.
                           n_samples_seen_ = NULL,

                           #' @field singular_values_ Numeric vector.
                           #'   Singular values associated with each principal component.
                           singular_values_ = NULL,

                           #' @field mean_ Numeric vector.
                           #'   The mean of each feature in the dataset, computed incrementally.
                           mean_ = NULL,

                           #' @field var_ Numeric vector.
                           #'   The variance of each feature in the dataset, computed incrementally.
                           var_ = NULL,

                           #' @field explained_variance_ratio_ Numeric vector.
                           #'   The proportion of the total variance explained by each principal component.
                           explained_variance_ratio_ = NULL,

                           #' @description initialize Initialize the Incremental Decomposition object.
                           #' @param X A matrix or data frame of data from which to compute the decomposition.
                           #' @param n_components The number of components to compute.
                           #' @details This method initializes the decomposition object, sets the number of
                           #'   components, computes initial statistics from the data (mean and variance),
                           #'   and performs an initial update of the singular value decomposition (SVD) based on `X`.
                           #'   The fields such as mean, variance, number of samples seen are initialized.
                           #'   The method sets up the object state necessary for further incremental updates.
                           initialize = function(X, n_components) {
                             self$n_components = n_components

                             n_samples <- nrow(X)
                             n_features <- ncol(X)

                             self$n_samples_seen_ <- n_samples
                             self$mean_ <- colMeans(X)
                             self$var_ <- apply(X, 2, var)

                             X <- sweep(X, 2, self$mean_, "-")
                             private$update_svd(X, self$mean_, self$var_, n_features, n_samples)
                           },

                           #' @description Partially fit the model to new data
                           #' @param X A matrix or data frame of new data samples to integrate into the model.
                           #' @param lambda A scaling factor applied to the singular values during the update.
                           #' @details This method updates the decomposition model incrementally with new data. It
                           #'   recalculates the mean and variance incrementally and adjusts the data matrix before
                           #'   performing an SVD update. The method integrates the new data `X` into the existing model by:
                           #'   - Updating the cumulative statistics (mean and variance).
                           #'   - Adjusting the new data by removing the batch mean.
                           #'   - Applying a correction to adjust for the cumulative mean difference.
                           #'   - Scaling the existing components by `lambda` times the singular values.
                           #'   - Calling the internal function `update_svd` to perform the singular value decomposition on the adjusted data.
                           #' @return None; updates the model's state in place.
                           partial_fit = function(X, lambda) {

                             n_samples <- nrow(X)
                             n_features <- ncol(X)

                             updated_stats <- incremental_mean_and_var(X, self$mean_, self$var_, rep(self$n_samples_seen_, ncol(X)))
                             col_mean = updated_stats[[1]]
                             col_var = updated_stats[[2]]
                             n_total_samples <- updated_stats[[3]][1]

                             col_batch_mean <- colMeans(X)
                             X <- sweep(X, 2, col_batch_mean, "-")
                             mean_correction <- sqrt((self$n_samples_seen_ / n_total_samples) * n_samples) * (self$mean_ - col_batch_mean)
                             X <- rbind(lambda * self$singular_values_ * self$components_, X, mean_correction)
                             private$update_svd(X, col_mean, col_var, n_features, n_total_samples)
                           },

                           #' @description Calculate the precision matrix of the model
                           #' @details This method computes the precision matrix using the components and explained variance
                           #'   from the model. Adjustments are made for noise variance, and the method involves:
                           #'   - Subtracting the noise variance from the explained variance and zeroing negative results.
                           #'   - Calculating an initial precision matrix and then adjusting its diagonal based on the
                           #'     inverse of the adjusted explained variance.
                           #'   - Refining the precision matrix by solving and scaling operations, and adjusting the diagonal
                           #'     again to account for noise variance.
                           #' @return A matrix representing the precision of the model.
                           get_precision = function(){
                             components_ <- self$components_
                             exp_var <- self$explained_variance_

                             exp_var_diff <- exp_var - self$noise_variance_
                             exp_var_diff[exp_var <= self$noise_variance_] <- 0.0

                             precision <- (components_ %*% t(components_)) / self$noise_variance_

                             diag(precision) <- diag(precision) + 1.0 / exp_var_diff

                             precision <- t(components_) %*% solve(precision) %*% components_
                             precision <- precision / -(self$noise_variance_^2)

                             # Add to diagonal again
                             diag(precision) <- diag(precision) + 1.0 / self$noise_variance_
                             return(precision)
                           },

                           #' @description Calculate the covariance matrix of the model
                           #' @details This method computes the covariance matrix using the model's principal components
                           #'   and explained variance. It adjusts for noise variance by:
                           #'   - Subtracting the noise variance from the explained variance and setting negative results to zero.
                           #'   - Multiplying the principal components by the adjusted variance to compute the covariance.
                           #'   - Adding the noise variance to the diagonal of the covariance matrix to account for data noise.
                           #' @return A matrix representing the covariance of the model.
                           get_covariance = function() {
                             components_ <- self$components_
                             exp_var <- self$explained_variance_

                             exp_var_diff <- exp_var - self$noise_variance_
                             # Adjust exp_var_diff: set to 0 where exp_var <= inc.pca$noise_variance_
                             exp_var_diff[exp_var <= self$noise_variance_] <- 0

                             # Compute covariance
                             cov <- t(components_) %*% diag(exp_var_diff) %*% components_

                             # Add noise variance to the diagonal
                             diag(cov) <- diag(cov) + self$noise_variance_

                             return(cov)
                           }),
                         private=list(
                           update_svd = function(X, col_mean, col_var, n_features, n_total_samples) {
                             svd_result <- svd(X, nu=self$n_components, nv=self$n_components)

                             U <- svd_result$u
                             S <- svd_result$d
                             Vt <- t(svd_result$v)

                             flipped = svd_flip(U, Vt, u_based_decision=F)
                             U <- flipped[[1]]
                             Vt <- flipped[[2]]

                             explained_variance <- (S^2) / (n_total_samples - 1)
                             explained_variance_ratio <- (S^2) / sum(col_var * n_total_samples)

                             self$n_samples_seen_ <- n_total_samples
                             self$components_ <- Vt[1:self$n_components,]
                             self$singular_values_ <- S[1:self$n_components]
                             self$mean_ <- col_mean
                             self$var_ <- col_var
                             self$explained_variance_ <- explained_variance[1:self$n_components]
                             self$explained_variance_ratio_ <- explained_variance_ratio[1:self$n_components]
                             if (self$n_components %in% c(n_total_samples, n_features)) {
                               self$noise_variance_ <- 0
                             } else {
                               self$noise_variance_ <- mean(explained_variance[(self$n_components + 1):length(explained_variance)])
                             }

                             invisible(self)
                           }
                         ))
