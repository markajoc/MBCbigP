#' @title Initialise parameters for Expectation-Maximisation algorithm
#'
#' @description Initialise the cluster membership matrix for the
#'   Expectation-Maximisation algorithm for fitting a mixture of Gaussians.
#'
#' @param x A data frame containing continuous data.
#' @param groups The number of groups.
#' @param method The method to give initial clustering, currently available:
#'   \code{"kmeans"} and \code{"hclust"}.
#' @param ... Arguments to pass to clustering functions.
#'
#' @return A \code{groups} x \code{nrow(x)} matrix of 0/1 entries denoting
#'   cluster membership.
#'
#' @examples
#' data(mtcars)
#' initialise.memberships(mtcars, 3)

initialise.memberships <-
function (x, groups, method = c("kmeans", "hclust", "mclust"), ...)
{
  x <- x[, apply(x, 2, function(obj) !any(is.na(obj))), drop = FALSE]
  method <- match.arg(method)
  if (identical(method, "kmeans")){
    cl <- stats::kmeans(x, groups, ...)
    clustering <- cl$cluster
  } else if (identical(method, "hclust")){
    d <- dist(x)
    h <- hclust(d, ...)
    clustering <- cutree(h, k = groups)
  } else if (identical(method, "mclust")){
    stop("initialisation method 'mclust' not currently supported")
  } else stop("initialisation method '", method, "' not recognised")
  z <- model.matrix(~ as.factor(clustering) - 1)
  colnames(z) <- 1:groups
  z[, , drop = FALSE]
}
