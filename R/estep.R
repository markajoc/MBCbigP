#' @title Expectation step for Expectation-Maximisation algorithm
#'
#' @description Expectation step for Expectation-Maximisation algorithm
#'
#' @param x A matrix of data.
#' @param parameters A list containing the parameters of the Gaussian mixture
#'   model.
#'
#' @return A matrix of probabilities of belonging to a cluster.

estep <-
function (x, parameters)
{
  K <- parameters$groups
  alpha <- parameters$mixing
  mu <- parameters$mean
  sigma <- parameters$sigma
  z <- matrix(nrow = nrow(x), ncol = K)
  for (k in 1:K){
    z[, k] <- alpha[k] * mvtnorm::dmvnorm(x = x, mean = mu[, k], sigma =
      sigma[, , k])
  }
  apply(z, 2, `/`, rowSums(z))
}
