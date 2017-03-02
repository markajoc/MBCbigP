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
function (x, parameters)
{
  if (is.null(parameters$pro))
    parameters$pro <- rep(1, parameters$groups)
  tmp <- matrix(as.double(NA), nrow(x), parameters$groups)
  for (k in 1:parameters$groups){
    tmp[, k] <- parameters$pro[k] * mvtnorm::dmvnorm(x, parameters$mean[, k],
      as.matrix(parameters$sigma[, , k]))
  }
  sum(log(rowSums(tmp)))
}

#' @rdname calcloglik

calcloglik_cond <-
function(x, x_cond, parameters, parameters_cond)
{
  if (is.null(parameters_cond$pro))
    parameters_cond$pro <- rep(1, parameters$groups)
  tmp <- matrix(as.double(NA), nrow(x), parameters$groups)
  for (k in 1:parameters$groups){
    for (i in 1:nrow(x)){
      tmp[i, k] <- parameters_cond$pro[k] *
        mvtnorm::dmvnorm(
          x = x[i, , drop = FALSE],
          mean = parameters$mean[, k],
          sigma = as.matrix(parameters$sigma[, , k])) *
        mvtnorm::dmvnorm(
          x = x_cond[i, , drop = FALSE],
          mean = parameters_cond$mean[i, , k, drop = TRUE],
          sigma = as.matrix(parameters_cond$sigma[, , k]))
    }
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
