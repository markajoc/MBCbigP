#' @title Weighted column means
#'
#' @description Calculate weighted column means of a matrix or dataframe.
#'
#' @param x A numeric matrix or dataframe.
#' @param w A vector of weights, must be same length as number of rows in
#'   \code{x}.
#' @param ... passed to \code{\link{weighted.mean}}
#'
#' @return A vector of weighted means.
#'
#' @examples
#' colMeans.weighted(mtcars, runif(nrow(mtcars)))

colMeans.weighted <-
function (x, w, ...)
{
  apply(x[, , drop = FALSE], 2, weighted.mean, w = w, ...)
}
