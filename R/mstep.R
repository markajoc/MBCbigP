mstep <-
function (x, z, groups = NULL, p = NULL)
{
  K <- if (is.null(groups))
    ncol(z)
  else groups
  p <- if (is.null(p))
    ncol(x)
  else p
  mixing <- colMeans(z)
  mu <- matrix(nrow = p, ncol = K)
  sigma <- array(dim = c(p, p, K))
  for (k in 1:K){
    mu[, k] <- apply(x, 2, weighted.mean, w = z[, k])
    sigma[, , k] <- cov.wt(x, wt = z[, k], method = "ML")$cov
  }
  list(mixing = mixing, mean = mu, sigma = sigma, groups = K)
}
