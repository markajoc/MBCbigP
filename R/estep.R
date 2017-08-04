#' @title Expectation step for Expectation-Maximisation algorithm
#'
#' @description Expectation step for Expectation-Maximisation algorithm
#'
#' @param x,x_A,x_B A matrix of data.
#' @param pro A vector of mixing proportions.
#' @param mean,mean_A,mean_B A 2-D array of mean vectors.
#' @param sigma,sigma_AA,sigma_AB,sigma_BB A 3-D array of covariance matrices.
#' @param groups The integer number of groups.
#' @param oldz Previous matrix of cluster probabilities.
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

#' @rdname estep

## TODO: try conditioning z matrix on more history

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
  if (any(apply(z, 2L, function(x) all(x == 0)))){
    cat("\n")
    #browser()
    stop("assigned all observations to zero in column of cluster memberships")
  }
  z / rowSums(z)
}

#' @rdname estep

estep_cond2 <-
function (x_A, x_B, pro, mean_A, mean_B, sigma_AA, sigma_AB, sigma_BB, groups,
  oldz = NULL)
{
  oldz <- if (is.null(oldz))
    matrix(1, nrow = nrow(x_A), ncol = groups)
  else oldz

  x <- cbind(data.matrix(x_A), data.matrix(x_B))
  mean <- rbind(data.matrix(mean_A), data.matrix(mean_B))
  z <- matrix(nrow = nrow(x), ncol = groups)
  for (k in 1:groups){
    sigma <- rbind(cbind(sigma_AA[, , k], sigma_AB[, , k]),
      cbind(t(sigma_AB[, , k]), sigma_BB[, , k]))
    z[, k] <- log(pro[k]) + mvtnorm::dmvnorm(x = x, mean = mean[, k], sigma =
      as.matrix(sigma), log = TRUE)
  }
  z <- (3 * z + log(oldz)) / 4
  z <- exp(z)
  if (any(apply(z, 2L, function(x) all(x == 0)))){
    cat("\n")
    stop("assigned all observations to zero in column of cluster memberships")
  }
  structure(z / rowSums(z), unscaled = z)
}
