#' @title Fit a mixture of multivariate Gaussians
#'
#' @description Fit a mixture of multivariate Gaussians
#'
#' @param x Data frame or a matrix
#' @param groups The number of groups/mixture components to fit.
#' @param maxiter The maximum number of iterations for the E-M algorithm.
#' @param likelihood Logical indicating whether the log-likelihood should be
#'   calculated at each step and returned (defaults to \code{TRUE}).
#'
#' @return A list containing the estimated parameters for the mixture
#'   distribution.
#'
#' @examples
#' library(mclust)
#' data(banknote)
#' mbc(x = banknote[, -1], groups = 2)

mbc <-
function (x, groups, maxiter = 500, likelihood = TRUE)
{
  x <- data.matrix(x)
  N <- nrow(x)
  p <- ncol(x)
  K <- groups

  ## Initialise z matrix.

  z <- initialise.memberships(x, K)

  ## Initialise empty means and covariance matrices.

  mu <- matrix(nrow = p, ncol = K)
  sigma <- array(dim = c(p, p, K))

  ## Initialise NULL log-likelihoods

  loglikprevious <- NULL
  loglik <- NULL

  for (times in 1:maxiter){

    ## Calculate maximum likelihood estimates for the mixing proportions, means
    ## and covariance matrices

    alpha <- colMeans(z)

    for (k in 1:K){
      mu[, k] <- apply(x, 2, weighted.mean, w = z[, k])
      sigma[, , k] <- cov.wt(x, wt = z[, k], method = "ML")$cov
    }

    ## Calculate log-likelihood.

    if (!is.null(loglik))
      loglikprevious <- loglik
    if (likelihood){
      tmp <- matrix(as.double(NA), N, K)
      for (k in 1:K){
        tmp[, k] <- alpha[k] * mvtnorm::dmvnorm(x, mu[, k], sigma[, , k])
      }
      loglik <- sum(log(rowSums(tmp)))
    }

    ## If we have a previous log-likelihood, check that we are not decreasing
    ## the log-likelihood. If we are not increasing it by much, break the loop.

    if (!is.null(loglikprevious)){
      stopifnot(loglik >= loglikprevious)
      if (!(loglik > loglikprevious + (abs(loglikprevious) * 1e-7)))
        break
    }

    ## Calculate expected z values

    for (k in 1:K){
      z[, k] <- alpha[k] * mvtnorm::dmvnorm(x = x, mean = mu[, k], sigma =
        sigma[, , k])
    }
    zrm <- rowSums(z)
    z <- apply(z, 2, `/`, zrm)
  }
  structure(list(mixing = alpha, mean = mu, sigma = sigma, groupprob = z, loglik
    = loglik), class = "mbc")
}
