## Formulas for parameter estimates.

estimate_sigma_BB_cathal <-
function (x_B, mu_BgivenA, z, sigma_AB, sigma_AA)
{
  w <- z / sum(z)
  crossprod(sqrt(w) * (x_B - mu_BgivenA)) + t(sigma_AB) %*% solve(sigma_AA) %*%
    sigma_AB
}

estimate_sigma_BB_cathal2 <-
function (x_B, mu_B, x_A, mu_A, sigma_AA, sigma_AB, z)
{
  w <- z / sum(z)
  wsq <- sqrt(w)
  sigma_AA_inverse <- solve(sigma_AA)
  sigma_BA_by_sigma_AA_inverse <- t(sigma_AB) %*% sigma_AA_inverse
  out <- crossprod(wsq * sweep(x_B, 2, mu_B)) + sigma_BA_by_sigma_AA_inverse %*%
    sigma_AB + crossprod(wsq * (sweep(x_A, 2, mu_A)) %*% t(sigma_BA_by_sigma_AA_inverse))
  out
}

estimate_sigma_BB_mark <-
function (x_B, mu_B, z)
{
  w <- z / sum(z)
  out <- crossprod(sqrt(w) * sweep(x_B, 2, mu_B))
  rownames(out) <- colnames(out) <- colnames(x_B)
  out
}

estimate_mu_B <-
function (x_B, z, sigma_AB, sigma_AA, x_A, mu_A)
{
  out <- colMeans.weighted(x_B, w = z) + t(sigma_AB) %*% solve(sigma_AA) %*%
    colMeans.weighted(sweep(x_A, 2, mu_A), w = z)
  names(out) <- colnames(x_B)
  out
}

estimate_sigma_AB <-
function (x_A, x_B, mu_A, mu_B, z)
{
  w <- z / sum(z)
  out <- crossprod(sqrt(w) * sweep(x_A, 2, mu_A), sqrt(w) * sweep(x_B, 2, mu_B))
  out <- 0 * out
  rownames(out) <- colnames(x_A)
  colnames(out) <- colnames(x_B)
  out
}

estimate_sigma_AB_simple <-
function (x_A, x_B, z)
{
  var.wt(x = x_A, y = x_B, w = z)
}

estimate_sigma_AB_corr <-
function (x_A, x_B, mu_A, mu_B, sigma_AA, sigma_BB, sigma_AB, pro, groups)
{
  rho_AB <- sigma_AB * NA
  for (k in 1:groups){
    rho_AB[, , k] <- diag(1 / sqrt(diag(sigma_AA[, , k]))) %*% sigma_AB[, , k] %*%
      diag(1 / sqrt(diag(sigma_BB[, , k])))
  }
  f <- function (x, dims, ...){
    rho_AB <- array(x, dim = dims)
    sigma_AB <- rho_AB * NA
    for (k in 1:groups){
      sigma_AB[, , k] <- diag(sqrt(diag(sigma_AA[, , k]))) %*% rho_AB[, , k] %*%
        diag(sqrt(diag(sigma_BB[, , k])))
    }
    -calcloglik_sigma_AB(x = sigma_AB, dims = dims, x_A = x_A, x_B = x_B, mean_A =
      mu_A, mean_B = mu_B, sigma_AA = sigma_AA, sigma_BB = sigma_BB, pro = pro,
      groups = groups)
  }
  #print(f(rho_AB, dim(sigma_AB)))
  est <- optim(par = rho_AB, fn = f, dims = dim(sigma_AB), x_A = x_A, x_B = x_B,
    mu_A = mu_A, mu_B = mu_B, sigma_AA = sigma_AA, sigma_BB = sigma_BB, pro =
    pro, groups = groups)#, lower = -1, upper = 1)
  rho_AB_new <- array(est$par, dim = dim(sigma_AB))
  #print(est$value)
  rho_AB_new
  for (k in 1:groups){
    sigma_AB[, , k] <- diag(sqrt(diag(sigma_AA[, , k]))) %*% rho_AB_new[, , k] %*%
      diag(sqrt(diag(sigma_BB[, , k])))
  }
  dimnames(sigma_AB) <- list(colnames(x_A), colnames(x_B), 1:groups)
  sigma_AB
}

estimate_sigma_AB_michael <-
function(x_A, x_B, mu_A, mu_B, sigma_AA, sigma_BB, z)
{
  w <- z / sum(z)
  W_AA <- crossprod(sqrt(w) * sweep(x_A, 2, mu_A))
  W_AB <- crossprod(sqrt(w) * sweep(x_A, 2, mu_A), sqrt(w) * sweep(x_B, 2,
    mu_B))
  #W_AA_inverse <- solve(W_AA)
  #out <- sigma_AA %*% W_AA_inverse %*% sigma_AA %*% sigma_AA %*% W_AB
  sigma_AA_inverse <- solve(sigma_AA)
  out <- solve(sigma_AA_inverse %*% W_AA %*% sigma_AA_inverse) %*% (sigma_AA %*%
    W_AB)
  dimnames(out) <- list(colnames(x_A), colnames(x_B))
  out
}
