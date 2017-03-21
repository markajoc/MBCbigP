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
  x <- data.matrix(x)
  groups <- if (is.null(groups))
    ncol(z)
  else groups
  p <- if (is.null(p))
    ncol(x)
  else p
  pro <- colMeans(z)
  mu <- matrix(nrow = p, ncol = groups)
  sigma <- array(dim = c(p, p, groups))
  for (k in 1:groups){
    mu[, k] <- apply(x, 2, weighted.mean, w = z[, k])
    sigma[, , k] <- cov.wt(x, wt = z[, k], method = "ML")$cov
  }
  rownames(mu) <- colnames(x)
  colnames(mu) <- 1:groups
  dimnames(sigma) <- list(colnames(x), colnames(x), 1:groups)
  list(pro = pro, mean = mu, sigma = sigma, groups = groups)
}

mstep_cond <-
function (x_A, x_B, mean_A, sigma_AA, groups, z)
{
  x_A <- data.matrix(x_A)
  x_B <- data.matrix(x_B)

  ## Prepare arrays.

  mu_B <- matrix(nrow = ncol(x_B), ncol = groups)
  sigma_BB <- array(0, dim = c(ncol(x_B), ncol(x_B), groups))
  sigma_AB <- array(dim = c(ncol(x_A), ncol(x_B), groups))

  ## For each group, calculate the covariance between batches, the mean adjusted
  ## for the group probabilities (z), and the covariance matrix for the current
  ## batch (should also be adjusted for z).

  for (k in 1:groups){
    mu_B[, k] <- estimate_mu_B(x_B = x_B, z = z[, k], sigma_AB = sigma_AB[, , k]
      , sigma_AA[, , k], x_A = x_A, mu_A = mean_A[, k])
    sigma_AB[, , k] <- estimate_sigma_AB(x_A = x_A, x_B = x_B, mu_A = mean_A[, k],
      mu_B = mu_B[, k], z = z[, k])
    sigma_BB[, , k] <- estimate_sigma_BB_cathal2(x_B = x_B, mu_B = mu_B[, k],
      x_A = x_A, mu_A = mean_A[, k], sigma_AA = sigma_AA[, , k], sigma_AB =
      sigma_AB[, , k], z = z[, k])
  }
  structure(list(pro = colMeans(z), mean = mu_B, sigma = sigma_BB, cov =
    sigma_AB, groups = groups), class = "mbcparameters")
}
