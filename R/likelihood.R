#' @title Calculate log-likelihood for a Gaussian mixture model.
#'
#' @description Calculate log-likelihood for a Gaussian mixture model.
#'
#' @param x A data frame or matrix of data.
#' @param parameters A list of parameters for the model.
#' @param object A fitted model object of class \code{"mbc"}.
#' @param ... not currently used.
#'
#' @return A single numeric value.

calcloglik <-
function (x, pro, mean, sigma, groups)
{
  tmp <- matrix(as.double(NA), nrow(x), groups)
  for (k in 1:groups){
    tmp[, k] <- pro[k] * mvtnorm::dmvnorm(x, mean[, k], as.matrix(sigma[, , k]))
  }
  sum(log(rowSums(tmp)))
}

#' @rdname calcloglik

calcloglik_split <-
function(x_A, x_B, pro, mean_A, mean_B, sigma_AA, sigma_AB, sigma_BB, groups)
{
  x <- cbind(data.matrix(x_A), data.matrix(x_B))
  mean <- rbind(data.matrix(mean_A), data.matrix(mean_B))
  tmp <- matrix(as.double(NA), nrow(x), groups)
  for (k in 1:groups){
    sigma <- rbind(cbind(sigma_AA[, , k], sigma_AB[, , k]),
      cbind(t(sigma_AB[, , k]), sigma_BB[, , k]))
    tmp[, k] <- pro[k] * mvtnorm::dmvnorm(x, mean[, k], as.matrix(sigma))
  }
  sum(log(rowSums(tmp)))
}

#' @rdname calcloglik

calcloglik_joint <- function(x1, x2, parameters1, parameters2)
{
  x <- cbind(x_A, x_B)
  mu <- c(mu_A, mu_B)
  sigma <- rbind(cbind(sigma_AA, sigma_AB), cbind(t(sigma_AB), sigma_BB))
}

#' @rdname calcloglik

logLik.mbc <-
function (object, ...)
{
  calcloglik(object)
}
