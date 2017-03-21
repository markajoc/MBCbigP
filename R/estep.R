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
function (x, pro, mean, sigma, groups)
{
  x <- data.matrix(x)
  z <- matrix(nrow = nrow(x), ncol = groups)
  for (k in 1:groups){
    z[, k] <- log(pro[k]) + mvtnorm::dmvnorm(x = x, mean = mean[, k], sigma =
    as.matrix(sigma[, , k]), log = TRUE)
  }
  z <- exp(z)
  z / rowSums(z)
}

estep_cond <-
function (x_A, x_B, pro, mean_A, mean_B, sigma_AA, sigma_AB, sigma_BB, groups)
{
  x <- cbind(data.matrix(x_A), data.matrix(x_B))
  mean <- rbind(data.matrix(mean_A), data.matrix(mean_B))
  z <- matrix(nrow = nrow(x), ncol = groups)
  for (k in 1:groups){
    sigma <- rbind(cbind(sigma_AA[, , k], sigma_AB[, , k]),
      cbind(t(sigma_AB[, , k]), sigma_BB[, , k]))
    z[, k] <- log(pro[k]) + mvtnorm::dmvnorm(x = x, mean = mean[, k], sigma =
    as.matrix(sigma), log = TRUE)
  }
  z <- exp(z)
  z / rowSums(z)
}
