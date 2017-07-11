#' @title Initialise parameters for Expectation-Maximisation algorithm
#'
#' @description Initialise the cluster membership matrix for the
#'   Expectation-Maximisation algorithm for fitting a mixture of Gaussians.
#'
#' @param x A data frame containing continuous data.
#' @param groups The number of groups.
#' @param method The method to give initial clustering.
#'
#' @return A \code{groups} x \code{nrow(x)} matrix of 0/1 entries denoting cluster
#'   membership.
#'
#' @examples
#' data(mtcars)
#' initialise.memberships(mtcars, 3)

initialise.memberships <-
function (x, groups, method = c("kmeans", "hclust"))
{
  x <- x[, apply(x, 2, function(obj) !any(is.na(obj))), drop = FALSE]
  cl <- stats::kmeans(x, groups)
  z <- model.matrix(~ as.factor(cl$cluster) - 1)
  colnames(z) <- 1:groups
  z[, , drop = FALSE]
}
