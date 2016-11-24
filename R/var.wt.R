#' @title Calculate a weighted covariance matrix between two sets of columns
#'
#' @description Calculate a weighted covariance matrix in the same way as
#'   \code{\link{cov.wt}}, but allow two sets of columns as in
#'   \code{\link{var}}.
#'
#' @param x Data frame or a matrix
#' @param y Data frame or a matrix with the same number of rows as \code{x}. If
#'   left \code{NULL}, the covariance matrix of \code{x} is calculated.
#' @param w A set of weights.
#' @param method The type of estimate return, \code{"unbiased"} (default) or
#'   \code{"ML"}, meaning maximum likelihood.
#'
#' @return A weighted covariance matrix.
#'
#' @examples
#' wts <- (tmp <- runif(nrow(mtcars))) / sum(tmp)
#' var.wt(mtcars[, 1:5], mtcars[, 1:2], w = wts)
#' cov.wt(mtcars[, 1:5], wt = wts)$cov

var.wt <-
function (x, y = NULL, w, method = c("unbiased", "ML")){

  ## Check inputs

  if (is.data.frame(x))
    x <- as.matrix(x)
  if (is.null(y))
    y <- x
  else if (is.data.frame(y))
    y <- as.matrix(y)
  n <- nrow(x)
  stopifnot(n == nrow(y))
  w <- if (missing(w))
    rep(1, n) / n
  else w / sum(w)

  ## Calculate the required sums for weighted covariance matrix

  xcenter <- colSums(w * x)
  sqw <- sqrt(w)
  x <- sqw * sweep(x, 2, xcenter)
  ycenter <- colSums(w * y)
  y <- sqw * sweep(y, 2, ycenter)

  ## Return either the unbiased or maximum likelihood (ML) covariance matrix

  switch(match.arg(method), unbiased = crossprod(x, y) / (1 - sum(w ^ 2)),
    ML = crossprod(x, y))
}
