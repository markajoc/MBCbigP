#' @title Calculate conditional mean and covariance for multivariate Gaussian
#'   distributions.
#'
#' @description Calculate conditional mean and covariance for multivariate
#'   Gaussian distributions. Consider partitioning into A and B. This function
#'   returns parameters for the distribution of A given B, which is itself
#'   another Gaussian distribution.
#'
#' @param mean The mean vector for A.
#' @param cov The covariance between A and B.
#' @param precision_cond The precision matrix for B.
#' @param x_cond The value of B on which to condition (only needed for the
#'   mean).
#' @param mean_cond The mean vector for B.
#' @param sigma_cond The covariance matrix for B, optional alternative to
#'   \code{precision_cond}.
#'
#' @return The relevant parameter for the distribution of A given B.

mean_conditional <-
function (mean, cov, precision_cond, x_cond, mean_cond, sigma_cond)
{
  precision_cond <- if (missing(precision_cond))
    solve(sigma_cond)
  else precision_cond
  if (!identical(dim(cov)[2L], ncol(precision_cond)))
    cov <- t(cov)
  t(mean + cov %*% precision_cond %*% t(sweep(x_cond, 2, mean_cond)))
}

#' @rdname mean_conditional

sigma_conditional <-
function (sigma, cov, precision_cond, sigma_cond)
{
  precision_cond <- if (missing(precision_cond))
    solve(sigma_cond)
  else precision_cond
  if (!identical(dim(cov)[2L], ncol(precision_cond)))
    cov <- t(cov)
  sigma - cov %*% precision_cond %*% t(cov)
}
