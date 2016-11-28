#' @title Maximisation step for Expectation-Maximisation algorithm
#'
#' @description Maximisation step for Expectation-Maximisation algorithm
#'
#' @param x A matrix of data.
#' @param z A matrix of cluster memberships/probabilities.
#' @param groups Optional, number of groups in the mixture, inferred from
#'   \code{z} if not supplied.
#' @param p Optional, dimension of \code{x}, inferred from \code{x} if not
#'   supplied
#'
#' @return A list of parameter estimates for the mixing proportions, mean
#'   vectors, and covariance matrices.


mstep <-
function (x, z, groups = NULL, p = NULL)
{
  groups <- if (is.null(groups))
    ncol(z)
  else groups
  p <- if (is.null(p))
    ncol(x)
  else p
  mixing <- colMeans(z)
  mu <- matrix(nrow = p, ncol = groups)
  sigma <- array(dim = c(p, p, groups))
  for (k in 1:groups){
    mu[, k] <- apply(x, 2, weighted.mean, w = z[, k])
    sigma[, , k] <- cov.wt(x, wt = z[, k], method = "ML")$cov
  }
  list(mixing = mixing, mean = mu, sigma = sigma, groups = groups)
}
