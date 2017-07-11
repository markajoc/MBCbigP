#' @title Calculate a weighted covariance matrix between two sets of columns
#'
#' @description Calculate a weighted covariance matrix in the same way as
#'   \code{\link{cov.wt}}, but allow two sets of columns as in
#'   \code{\link{var}}.
#'
#' @param x Data frame or a matrix
#' @param y Data frame or a matrix with the same number of rows as \code{x}. If
#'   left \code{NULL}, the covariance matrix of \code{x} is calculated.
#' @param w A vector of weights. \code{w} will be divided by \code{sum(w)} to
#'   ensure weights sum to 1.
#' @param method The type of estimate returned, \code{"unbiased"} (default) or
#'   \code{"ML"}, meaning maximum likelihood.
#' @param xcenter,ycenter Optional center values for \code{x} and \code{y}.
#'
#' @return A weighted covariance matrix.
#'
#' @examples
#' wts <- (tmp <- runif(nrow(mtcars))) / sum(tmp)
#' var.wt(mtcars[, 1:5], mtcars[, 1:2], w = wts)
#' cov.wt(mtcars[, 1:5], wt = wts)$cov

var.wt <-
function (x, y = NULL, w, method = c("unbiased", "ML"), xcenter = NULL, ycenter
  = NULL){

  ## Check inputs

  x <- data.matrix(x)

  if (is.null(y)){
    y <- x
    ycenter <- xcenter
  } else {
    y <- data.matrix(y)
  }

  n <- nrow(x)

  stopifnot(n == nrow(y))

  w <- if (missing(w))
    rep(1, n) / n
  else w / sum(w)

  ## Calculate the required sums for weighted covariance matrix

  xcenter <- if (is.null(xcenter))
    colSums(w * x)
  else xcenter
  ycenter <- if (is.null(ycenter))
    colSums(w * y)
  else ycenter

  sqw <- sqrt(w)
  x <- sqw * sweep(x, 2, xcenter)
  y <- sqw * sweep(y, 2, ycenter)

  ## Return either the unbiased or maximum likelihood (ML) covariance matrix

  out <- switch(match.arg(method), unbiased = crossprod(x, y) / (1 - sum(w ^ 2))
    , ML = crossprod(x, y))

  rownames(out) <- colnames(x)
  colnames(out) <- colnames(y)
  out
}
