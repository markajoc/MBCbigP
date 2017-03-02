#' @title Calculate conditional mean and covariance for multivariate Gaussian
#'   distributions.
#'
#' @description Calculate conditional mean and covariance for multivariate
#'   Gaussian distributions. Consider partitioning into A and B. This function
#'   returns parameters for the distribution of B given A, which is itself
#'   another Gaussian distribution.
#'
#' @param mu_B The mean vector for B.
#' @param mu_A The mean vector for A.
#' @param x_A The points in A at which to evaluate the conditional mean of B.
#' @param sigma_AA The covariance matrix for A.
#' @param sigma_AB The cross-covariance matrix for A and B.
#'
#' @return The relevant mean parameter for the distribution of B given A.

mean_conditional <-
function (mu_B, mu_A, x_A, sigma_AA, sigma_AB)
{
  ## mu_BgivenA = mu_B + sigma_BA %*% solve(sigma_AA) %*% (x_A - mu_A)

  t(mu_B + t(sigma_AB) %*% solve(sigma_AA) %*% t(sweep(x_A, 2, mu_A)))
}

#' @rdname mean_conditional

sigma_conditional <-
function (sigma_BB, sigma_AB, sigma_AA)
{
  sigma_BB - t(sigma_AB) %*% solve(sigma_AA) %*% sigma_AB
}
