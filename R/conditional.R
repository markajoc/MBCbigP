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
#' @param x_A,x_B A dataframe or matrix
#' @param sigma_AA,sigma_BB The covariance matrix for batches A and B.
#' @param sigma_AB The cross-covariance matrix for A and B.
#' @param log Logical switch for log-likelihood.
#'
#' @return The relevant mean parameter for the distribution of B given A.

mean_conditional <-
function (mu_B, mu_A, x_A, sigma_AA, sigma_AB)
{
  t(mu_B + t(sigma_AB) %*% solve(sigma_AA) %*% t(sweep(x_A, 2, mu_A)))
}

#' @rdname mean_conditional

sigma_conditional <-
function (sigma_BB, sigma_AB, sigma_AA)
{
  sigma_BB - t(sigma_AB) %*% solve(sigma_AA) %*% sigma_AB
}

#' @rdname mean_conditional

dmvnorm_conditional <-
function(x_A, x_B, mu_A, mu_B, sigma_AA, sigma_AB, sigma_BB, log = TRUE)
{
  x <- cbind(x_A, x_B)
  mu <- c(mu_A, mu_B)
  sigma <- rbind(cbind(sigma_AA, sigma_AB), cbind(t(sigma_AB), sigma_BB))
  joint <- mvtnorm::dmvnorm(x = x, mean = mu, sigma = sigma, log = TRUE)
  marginal <- mvtnorm::dmvnorm(x = x_A, mean = mu_A, sigma = sigma_AA, log =
    TRUE)
  if (log)
    joint - marginal
  else exp(joint - marginal)
}
