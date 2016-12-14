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
  tmp <- matrix(as.double(NA), nrow(x), parameters$groups)
  for (k in 1:parameters$groups){
    tmp[, k] <- parameters$mixing[k] * mvtnorm::dmvnorm(x, parameters$mean[, k],
      as.matrix(parameters$sigma[, , k]))
  }
  sum(log(rowSums(tmp)))
}

#' @rdname calcloglik

logLik.mbc <-
function (object, ...)
{
  calcloglik(object)
}
