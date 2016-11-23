#' @title Fit a mixture of multivariate Gaussians
#'
#' @description Fit a mixture of multivariate Gaussians
#'
#' @param x Data frame or a matrix
#' @param groups The number of groups/mixture components to fit.
#' @param maxiter The maximum number of iterations for the E-M algorithm.
#'
#' @return A list containing the estimated parameters for the mixture
#'   distribution.
#'
#' @examples
#' library(mclust)
#' data(banknote)
#' mbc(x = banknote[, -1], groups = 2)

mbc <-
function (x, groups, maxiter = 500)
{
  x <- data.matrix(x)
  N <- nrow(x)
  p <- ncol(x)
  K <- groups

  ## Initialise mixing proportions

  alpha <- alpha1 <- rep(1 / K, K)

  ## Initialise means and covariance matrices (need to be positive-definite).
  ## Variables tagged with 1 are for the parallel calculations using `mclust`
  ## for comparison.

  mu <- mu1 <- t(stats::kmeans(x, K)$centers)

  sigma <- array(dim = c(p, p, K))
  tmp <- mclust::mclustVariance("VVV", d = p, G = K)$cholsigma
  for (k in 1:K){
    sigma[, , k] <- 10 * tmp[, , k] %*% t(tmp[, , k])
  }

  ## Initialise z matrix

  z <- matrix(ncol = K, nrow = N)

  ## Initialise NULL log-likelihoods

  loglikprevious <- NULL
  loglik <- NULL

  for (times in 1:maxiter){

    ## Calculate expected z values

    for (k in 1:K){
      z[, k] <- alpha[k] * mvtnorm::dmvnorm(x = x, mean = mu[, k], sigma =
        sigma[, , k])
    }
    zrm <- rowSums(z)
    z <- apply(z, 2, `/`, zrm)

    ## Using `mclust`. Needs the upper-triangular Cholesky decomposition of the
    ## covariance matrices.

    #tmp <- apply(sigma1, 3L, chol)
    #cholsigma <- array(0, dim = c(p, p, K))
    #for (k in 1:K){
    #  cholsigma[, , k] <- tmp[, k]
    #}
    #z1 <- mclust::estep("VVV", data = x, parameters = list(pro = alpha1, mean =
    #  mu1, variance = list(cholsigma = cholsigma)))$z


    ## Calculate maximum likelihood estimates for the mixing proportions, means
    ## and covariance matrices

    alpha <- colMeans(z)

    for (k in 1:K){
      mu[, k] <- apply(x, 2, weighted.mean, w = z[, k])
      sigma[, , k] <- cov.wt(x, wt = z[, k], method = "ML")$cov
    }

    ## Using `mclust`.

    #params <- mclust::mstep("VVV", x, z1)

    #alpha1 <- params$parameters$pro
    #mu1 <- params$parameters$mean
    #sigma1 <- params$parameters$variance$sigma

    if (!is.null(loglik))
      loglikprevious <- loglik

    ## Calculate log-likelihood.

    #tmp <- matrix(as.double(NA), N, K)
    #for (k in 1:K){
    #  tmp[, k] <- alpha[k] * mvtnorm::dmvnorm(x, mu[, k], sigma[, , k])
    #}
    #loglik <- sum(log(rowSums(tmp)))

    ## If we have a previous log-likelihood, check that we are not decreasing
    ## the log-likelihood. If we are not increasing it by much, break the loop.

    if (!is.null(loglikprevious)){
      stopifnot(loglik >= loglikprevious)
      if (!(loglik > loglikprevious + (abs(loglikprevious) * 1e-7)))
        break
    }
  }
  structure(list(mixing = alpha, mean = mu, sigma = sigma, groupprob = z),
    class = "mbc")
}
