#' @title Fit a mixture of multivariate Gaussians using sequential conditioning
#'
#' @description Fit a mixture of multivariate Gaussians
#'
#' @param x Data frame or a matrix
#' @param groups The number of groups/mixture components to fit.
#' @param batches The number of batches to split the columns of \code{x} for
#'   processing.
#' @param maxiter The maximum number of iterations for the E-M algorithm.
#'
#' @return A list containing the estimated parameters for the mixture
#'   distribution.
#'
#' @examples
#' library(mclust)
#' data(banknote)
#' mbc(x = banknote[, -1], groups = 2)

mbcbigp <-
function (x, groups, batches = 3, maxiter = 200)
{
  x <- data.matrix(x)
  N <- nrow(x)
  p <- ncol(x)
  K <- groups
  Q <- batches

  batches <- sort(rep(1:Q, length.out = ncol(x)))

  ## Initialise mixing proportions

  alpha <- rep(1 / K, K)

  ## Initialise membership probabilities

  z <- matrix(nrow = N, ncol = K)

  ## Initialise mean vectors

  mu <- matrix(nrow = K, ncol = p)

  ## Initialise covariance matrices

  sigma <- array(0, dim = c(p, p, K))

  ## Do the first batch using marginal density.

  batchindex <- which(batches == 1L)

  mu[, batchindex] <- (stats::kmeans(x[, batchindex, drop = FALSE], K)$centers)

  tmp <- mclust::mclustVariance("VVV", d = p, G = K)$cholsigma
  for (k in 1:K){
    sigma[batchindex, batchindex, k] <- 10 * tmp[batchindex, batchindex, k] %*%
      t(tmp[batchindex, batchindex, k])
  }
  rm(tmp)

  for (times in 1:maxiter){

    ## Expectation

    for (k in 1:K){
      z[, k] <- alpha[k] * mvtnorm::dmvnorm(x = x[, batchindex, drop = FALSE],
        mean = mu[k, batchindex, drop = FALSE], sigma = sigma[batchindex,
        batchindex, k])
    }
    z <- apply(z, 2, `/`, rowSums(z))

    ## Maximise

    alpha <- colMeans(z)

    for (k in 1:K){
      mu[k, batchindex] <- apply(x[, batchindex, drop = FALSE], 2,
        weighted.mean, w = z[, k])
      sigma[batchindex, batchindex, k] <- cov.wt(x[, batchindex, drop = FALSE],
        wt = z[, k], method = "ML")$cov
    }
  }

  ## Other batches

  for (q in 2:Q){

    ## Expectation from previous batch

    for (k in 1:K){
      z[, k] <- alpha[k] * mvtnorm::dmvnorm(x = x[, batchindex, drop = FALSE],
        mean = mu[k, batchindex, drop = FALSE], sigma = sigma[batchindex,
        batchindex, k])
    }
    z <- apply(z, 2, `/`, rowSums(z))

    alpha <- colMeans(z)

    ## Set up batch indices

    prevbatchindex <- batchindex
    batchindex <- which(batches == q)
    batchsize <- length(batchindex)

    ## Set up arrays for the off-diagonal covariance matrix, and the conditional
    ## parameter estimates

    cov_q_qminus1 <- array(dim = c(length(prevbatchindex), batchsize,
      K))
    mu_conditional <- matrix(ncol = batchsize, nrow = K)
    sigma_conditional <- array(dim = c(batchsize, batchsize, K))

    ## Calculate the inverse of the covariance of the previous batch.

    precision <- list()
    for (k in 1:K){
      precision[[k]] <- solve(sigma[prevbatchindex, prevbatchindex, k])
    }

    ## Start iterations.

    for (times in 1:maxiter){

      ## Calculate the conditional means and covariances.

      for (k in 1:K){
        mu_conditional[k, ] <- apply(x[, batchindex, drop = FALSE], 2,
          weighted.mean, w = z[, k])
        sigma_conditional[, , k] <- cov.wt(x[, batchindex, drop = FALSE],
          wt = z[, k], method = "ML")$cov
      }

      ## Expectation.

      for (k in 1:K){
        z[, k] <- alpha[k] * mvtnorm::dmvnorm(x = x[, batchindex, drop = FALSE],
          mean = mu_conditional[k, ], sigma = sigma_conditional[, , k])
      }
      z <- apply(z, 2, `/`, rowSums(z))

      ## Calculate parameter estimates

      for (k in 1:K){

        ## Covariance between current batch and previous batch.

        sigma[prevbatchindex, batchindex, k] <-
        sigma[batchindex, prevbatchindex, k] <-
        cov_q_qminus1[, , k] <- var.wt(x[, prevbatchindex], x[, batchindex], w =
          z[, k], method = "ML")

        ## Mean vectors

        extrabit <- - cov_q_qminus1[, , k] %*% precision[[k]] %*% colSums(z[, k]
          * (x[, prevbatchindex, drop = FALSE] - mu[rep(k, N), prevbatchindex, drop =
          FALSE]))
        print(extrabit)

        mu[k, batchindex] <- mu_conditional[k, ]

        ## Covariance matrices

        sigma[batchindex, batchindex, k] <- sigma_conditional[, , k] +
          t(cov_q_qminus1[, , k]) %*% precision[[k]] %*% cov_q_qminus1[, , k]
      }
    }
  }
  structure(list(mean = mu, sigma = sigma, groupprob = z), class = "mbc")
}
