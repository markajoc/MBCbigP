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
    sigma[, , k] <- var.wt(x, w = z[, k], method = "ML")
  }
  rownames(mu) <- colnames(x)
  colnames(mu) <- 1:groups
  dimnames(sigma) <- list(colnames(x), colnames(x), 1:groups)
  list(pro = pro, mean = mu, sigma = sigma, groups = groups)
}

#' @title Conditional maximisation step for Expectation-Maximisation algorithm
#'
#' @description Conditional maximisation step for Expectation-Maximisation
#'   algorithm. Maximising paramters for current batch of data conditional on
#'   the parameters from the previous batch.
#'
#' @param x_A A matrix of data for previous batch.
#' @param x_B A matrix of data for the current batch.
#' @param mean_A The mean vectors for the previous batch (considered fixed).
#' @param sigma_AA The covariance matrices for the previous batch (considered
#'   fixed).
#' @param sigma_AB Starting value for the covariance between batches A and B. If
#'   \code{NULL}, initialised at all zeros.
#' @param mean_B Starting value for the mean vectors for the current batch.
#' @param z A matrix of cluster memberships/probabilities.
#' @param groups Optional, number of groups in the mixture, inferred from
#'   \code{z} if not supplied.
#' @param p Optional, dimension of \code{x}, inferred from \code{x} if not
#'   supplied
#'
#' @return A list of parameter estimates for the mixing proportions, mean
#'   vectors, and covariance matrices.

mstep_cond <-
function (x_A, x_B, mean_A, sigma_AA, sigma_AB = NULL, mean_B = NULL,
  sigma_BB = NULL, pro = NULL, groups, z, likelihood = FALSE, method_sigma_AB =
  c("numeric", "analytic"))
{
  x_A <- data.matrix(x_A)
  x_B <- data.matrix(x_B)

  ## Prepare arrays.

  mu_B <- if (is.null(mean_B))
    matrix(nrow = ncol(x_B), ncol = groups)
  else mean_B
  sigma_BB <- if (is.null(sigma_BB))
    array(dim = c(ncol(x_B), ncol(x_B), groups))
  else sigma_BB
  sigma_AB <- if (is.null(sigma_AB))
    array(0, dim = c(ncol(x_A), ncol(x_B), groups))
  else sigma_AB

  ## For each group, calculate the covariance between batches, the mean adjusted
  ## for the group probabilities (z), and the covariance matrix for the current
  ## batch (should also be adjusted for z).

  method_sigma_AB <- match.arg(method_sigma_AB)
  analytic <- if (identical(method_sigma_AB, "analytic"))
    TRUE
  else FALSE

  for (k in 1:groups){
    mu_B[, k] <- estimate_mu_B(
      x_B = x_B,
      z = z[, k],
      sigma_AB = sigma_AB[, , k],
      sigma_AA[, , k],
      x_A = x_A,
      mu_A = mean_A[, k])
    muBgivenA <- mean_conditional(
      mu_B = mu_B[, k],
      mu_A = mean_A[, k],
      x_A = x_A,
      sigma_AA = sigma_AA[, , k],
      sigma_AB = sigma_AB[, , k])
    sigma_BB[, , k] <- estimate_sigma_BB_cathal(
      x_B = x_B,
      mu_BgivenA = muBgivenA,
      z = z[, k],
      sigma_AB = sigma_AB[, , k],
      sigma_AA = sigma_AA[, , k])
    if (analytic){
      sigma_AB[, , k] <- estimate_sigma_AB_michael(
        x_A = x_A,
        x_B = x_B,
        mu_A = mean_A[, k],
        mu_B = mu_B[, k],
        sigma_AA = sigma_AA[, , k],
        sigma_BB = sigma_BB[, , k],
        z = z[, k])
    }
  }

  if (!analytic){
    sigma_AB <- estimate_sigma_AB_corr(
      x_A = x_A,
      x_B = x_B,
      mu_A = mean_A,
      mu_B = mu_B,
      sigma_AA = sigma_AA,
      sigma_BB = sigma_BB,
      sigma_AB = sigma_AB,
      pro = pro,
      groups = groups)
  }

  dimnames(mu_B) <- list(colnames(x_B), 1:groups)
  dimnames(sigma_BB) <- list(colnames(x_B), colnames(x_B), 1:groups)
  dimnames(sigma_AB) <- list(colnames(x_A), colnames(x_B), 1:groups)

  mu_A_new <- NA * mean_A
  sigma_AA_new <- NA * sigma_AA
  for (k in 1:groups){
    mu_A_new[, k] <- estimate_mu_B(
      x_B = x_A,
      z = z[, k],
      sigma_AB = t(sigma_AB[, , k]),
      sigma_AA = sigma_BB[, , k],
      x_A = x_B,
      mu_A = mu_B[, k])
    muBgivenA <- mean_conditional(
      mu_B = mu_A_new[, k],
      mu_A = mu_B[, k],
      x_A = x_B,
      sigma_AA = sigma_BB[, , k],
      sigma_AB = t(sigma_AB[, , k]))
    sigma_AA_new[, , k] <- estimate_sigma_BB_cathal(
      x_B = x_A,
      mu_BgivenA = muBgivenA,
      z = z[, k],
      sigma_AB = t(sigma_AB[, , k]),
      sigma_AA = sigma_BB[, , k])

  }

  loglik <- if (likelihood){
    calcloglik_split(
      x_A = x_A,
      x_B = x_B,
      pro = pro,
      mean_A = mean_A,
      mean_B = mu_B,
      sigma_AA = sigma_AA,
      sigma_AB = sigma_AB,
      sigma_BB = sigma_BB,
      groups = groups)
  } else NULL
  structure(list(pro = colMeans(z), mean = mu_B, sigma = sigma_BB, cov =
    sigma_AB, meanprev = mu_A_new, sigmaprev = sigma_AA_new, groups = groups),
    class = "mbcparameters")
}
