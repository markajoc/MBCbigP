#' @title Fit a mixture of multivariate Gaussians using sequential conditioning
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

  print(mu)

  tmp <- mclust::mclustVariance("VVV", d = p, G = K)$cholsigma
  for (k in 1:K){
    sigma[batchindex, batchindex, k] <- 10 * tmp[batchindex, batchindex, k] %*%
      t(tmp[batchindex, batchindex, k])
  }
  rm(tmp)

  for (times in 1:maxiter){
    for (k in 1:K){
      z[, k] <- alpha[k] * mvtnorm::dmvnorm(x = x[, batchindex, drop = FALSE],
        mean = mu[k, batchindex, drop = FALSE], sigma = sigma[batchindex,
        batchindex, k])
    }
    zrm <- rowSums(z)
    z <- apply(z, 2, `/`, zrm)

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

    for (k in 1:K){
      z[, k] <- alpha[k] * mvtnorm::dmvnorm(x = x[, batchindex, drop = FALSE],
        mean = mu[k, batchindex, drop = FALSE], sigma = sigma[batchindex,
        batchindex, k])
    }
    zrm <- rowSums(z)
    z <- apply(z, 2, `/`, zrm)

    alpha <- colMeans(z)

    prevbatchindex <- batchindex
    batchindex <- which(batches == q)
    batchsize <- length(batchindex)
    cov_q_qminus1 <- array(dim = c(length(prevbatchindex), batchsize,
      K))
    mu_conditional <- matrix(ncol = batchsize, nrow = K)
    sigma_conditional <- array(dim = c(batchsize, batchsize, K))
    precision <- list()
    for (k in 1:K){
      precision[[k]] <- solve(sigma[prevbatchindex, prevbatchindex, k])
    }
    for (times in 1:maxiter){
      for (k in 1:K){
        sigma[prevbatchindex, batchindex, k] <-
        sigma[batchindex, prevbatchindex, k] <- cov_q_qminus1[, , k] <- cov.wt(x,
          wt = z[, k], method = "ML")$cov[prevbatchindex, batchindex]

        mu_conditional[k, ] <- apply(x[, batchindex, drop = FALSE], 2,
          weighted.mean, w = z[, k])

        #extrabit <- - cov_q_qminus1[, , k] %*% precision[[k]] %*%
        #  colMeans(x[, prevbatchindex] - mu[k, prevbatchindex])
        #print(extrabit)
        mu[k, batchindex] <- mu_conditional[k, ]

        sigma_conditional[, , k] <- cov.wt(x[, batchindex, drop = FALSE],
          wt = z[, k], method = "ML")$cov
        sigma[batchindex, batchindex, k] <- sigma_conditional[, , k] +
          t(cov_q_qminus1[, , k]) %*% precision[[k]] %*% cov_q_qminus1[, , k]
      }
    }
  }
  structure(list(mean = mu, sigma = sigma, groupprob = z), class = "mbc")
}
