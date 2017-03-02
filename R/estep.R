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
  x <- data.matrix(x)
  z <- matrix(nrow = nrow(x), ncol = parameters$groups)
  for (k in 1:parameters$groups){
    z[, k] <- log(parameters$pro[k]) + mvtnorm::dmvnorm(x = x, mean =
      parameters$mean[, k], sigma = as.matrix(parameters$sigma[, , k]), log =
      TRUE)
  }
  z <- exp(z)
  z / rowSums(z)
}

estep_cond <-
function (x1, x2, parameters1, parameters2)
{
  x1 <- data.matrix(x1)
  x2 <- data.matrix(x2)
  N <- nrow(x1)
  z <- matrix(nrow = N, ncol = parameters1$groups)
  for (k in 1:parameters1$groups){
    sigma_cond <- sigma_conditional(sigma = parameters1$sigma[, , k], cov =
      parameters1$cov[, , k], precision_cond = parameters1$precision2[, , k])
    mean_cond <- mean_conditional(mean = parameters1$mean[, k], cov =
      parameters1$cov[, , k], precision_cond = parameters1$precision2[, , k],
      x_cond = x2[, , drop = FALSE], mean_cond = parameters2$mean[, k],
      sigma_cond = parameters2$sigma[, , k])
    for (i in 1:N){
      z[i, k] <- log(parameters1$pro[k]) + mvtnorm::dmvnorm(x = x1[i, ], mean =
        mean_cond[i, ], sigma = as.matrix(sigma_cond), log = TRUE) +
        mvtnorm::dmvnorm(x = x2[i, ], mean = parameters2$mean[, k], sigma =
        as.matrix(parameters2$sigma[, , k]), log = TRUE)
    }
  }
  z <- exp(z)
  z / rowSums(z)
}
