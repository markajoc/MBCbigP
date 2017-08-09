## Formulas for parameter estimates.

estimate_sigma_BB_cathal <-
function (x_B, mu_BgivenA, z, sigma_AB, sigma_AA)
{
  w <- z / sum(z)
  crossprod(sqrt(w) * (x_B - mu_BgivenA)) + t(sigma_AB) %*% chol2inv(chol(
    sigma_AA)) %*% sigma_AB
}

estimate_sigma_BB_michael <-
function (x_B, mu_B, x_A, mu_A, sigma_AA, sigma_AB, z)
{
  w <- z / sum(z)
  wsq <- sqrt(w)
  sigma_AA_inverse <- chol2inv(chol(sigma_AA))
  x_A_cen <- wsq * sweep(x_A, 2, mu_A)
  x_B_cen <- wsq * sweep(x_B, 2, mu_B)
  W_AA <- crossprod(x_A_cen)
  W_BB <- crossprod(x_B_cen)
  W_AB <- crossprod(x_A_cen, x_B_cen)
  prec_BB <- (t(sigma_AB) %*% sigma_AA_inverse %*% W_AA %*% sigma_AA_inverse %*%
    sigma_AB - 2 * t(W_AB) %*% sigma_AA_inverse %*% sigma_AB + W_BB)
  out <- prec_BB + t(sigma_AB) %*% sigma_AA_inverse %*% sigma_AB
  rownames(out) <- colnames(out) <- colnames(x_B)
  out
}

## Maximimum likelihood estimate of part of mean vector, given fixed values of
## the rest of the mean vector, and covariance matrix.

estimate_mu_B <-
function (x_B, z, sigma_AB, sigma_AA, x_A, mu_A)
{
  x_B <- data.matrix(x_B)
  x_A <- data.matrix(x_A)
  out <- colMeans.weighted(x_B, w = z) - t(sigma_AB) %*% chol2inv(chol(
    sigma_AA)) %*% colMeans.weighted(sweep(x_A, 2, mu_A), w = z)
  if (any(is.na(out))){
    cat("\n")
    stop("conditional mean estimate contains NA/NaN")
  }
  names(out) <- colnames(x_B)
  out
}

## Numerically optimise the sigma_AB term by finding the correlations that
## maximise the likelihood, given fixed covariance matrices sigma_AA and
## sigma_BB.

estimate_sigma_AB_numerical <-
function (x_A, x_B, mu_A, mu_B, sigma_AA, sigma_BB, sigma_AB, pro, groups)
{
  ## Convert the off-diagonal block of covariance values to correlation values.

  rho_AB <- sigma_AB * NA
  for (k in 1:groups){
    rho_AB[, , k] <- diag(1 / sqrt(diag(sigma_AA[, , k]))) %*%
      sigma_AB[, , k] %*% diag(1 / sqrt(diag(sigma_BB[, , k])))
  }

  ## Objective function to optimise. Joint Gaussian likelihood of both batches
  ## as function of off-diagonal block correlations.

  f <- function (x, dims, ...){
    rho_AB <- array(x, dim = dims)
    sigma_AB <- rho_AB * NA
    for (k in 1:groups){
      sigma_AB[, , k] <- diag(sqrt(diag(sigma_AA[, , k]))) %*% rho_AB[, , k] %*%
        diag(sqrt(diag(sigma_BB[, , k])))
    }
    -calcloglik_sigma_AB(x = sigma_AB, dims = dims, x_A = x_A, x_B = x_B,
      mean_A = mu_A, mean_B = mu_B, sigma_AA = sigma_AA, sigma_BB = sigma_BB,
      pro = pro, groups = groups)
  }

  ## Use `optim` to numerically solve for optimal correlations.

  est <- optim(par = rho_AB, fn = f, dims = dim(sigma_AB), x_A = x_A, x_B = x_B,
    mu_A = mu_A, mu_B = mu_B, sigma_AA = sigma_AA, sigma_BB = sigma_BB, pro =
    pro, groups = groups)

  ## Convert correlations back to covariance values.

  rho_AB_new <- array(est$par, dim = dim(sigma_AB))
  for (k in 1:groups){
    sigma_AB[, , k] <- diag(sqrt(diag(sigma_AA[, , k]))) %*%
      rho_AB_new[, , k] %*% diag(sqrt(diag(sigma_BB[, , k])))
  }
  dimnames(sigma_AB) <- list(colnames(x_A), colnames(x_B), 1:groups)
  sigma_AB
}

## Analytic maximum likelihood for off-diagonal covariance matrix block.

estimate_sigma_AB_analytic <-
function(x_A, x_B, mu_A, mu_B, sigma_AA, sigma_BB, z)
{
  ## Calculate the required weighted cross products.

  w <- z / sum(z)
  wsq <- sqrt(w)
  W_AA <- crossprod(wsq * sweep(x_A, 2, mu_A))
  W_AB <- crossprod(wsq * sweep(x_A, 2, mu_A), wsq * sweep(x_B, 2, mu_B))

  ## Calculate the estimate.

  W_AA_inverse <- chol2inv(chol(W_AA))
  out <- sigma_AA %*% W_AA_inverse %*% W_AB

  if (any(is.na(out))){
    cat("\n")
    stop("between-batch covariance estimate contains NA/NaN")
  }

  ## Check that this estimate implies sensible correlation values in [-1, 1].

  rho_AB <- out * NA
  rho_AB <- diag(1 / sqrt(diag(sigma_AA))) %*% out %*% diag(1 / sqrt(diag(
    sigma_BB)))

  if (any(!inrange(rho_AB, c(-1, 1)))){
    cat("\n")
    stop("between-batch covariance implies invalid correlation outside [-1, 1]")
  }

  dimnames(out) <- list(colnames(x_A), colnames(x_B))
  structure(out, cor = rho_AB)
}
