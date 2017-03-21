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
  out <- crossprod(wsq * sweep(x_B, 2, mu_B)) - sigma_BA_by_sigma_AA_inverse %*%
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
  out <- colMeans.weighted(x_B, w = z) - t(sigma_AB) %*% solve(sigma_AA) %*%
  colMeans.weighted(sweep(x_A, 2, mu_A), w = z)
  names(out) <- colnames(x_B)
  out
}

estimate_sigma_AB <-
function (x_A, x_B, mu_A, mu_B, z)
{
  w <- z / sum(z)
  out <- crossprod(sqrt(w) * sweep(x_A, 2, mu_A), sqrt(w) * sweep(x_B, 2, mu_B))
  rownames(out) <- colnames(x_A)
  colnames(out) <- colnames(x_B)
  out
}
