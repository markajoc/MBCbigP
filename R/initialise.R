#' @title Initialise parameters for Expectation-Maximisation algorithm
#'
#' @description Initialise the cluster membership matrix for the
#'   Expectation-Maximisation algorithm for fitting a mixture of Gaussians.
#'
#' @param x A data frame containing continuous data.
#' @param K The number of groups.
#'
#' @return A \code{K} x \code{nrow(x)} matrix of 0/1 entries denoting cluster
#'   membership.
#'
#' @examples
#' data(mtcars)
#' initialise.memberships(mtcars, 3)

initialise.memberships <-
function (x, K)
{
  cl <- stats::kmeans(x, K)
  z <- model.matrix(~ as.factor(cl$cluster) - 1)
  colnames(z) <- 1:K
  z[, , drop = FALSE]
}
