var.wt <-
function (x, y = NULL, w, method = c("unbiased", "ML")){
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
  else w
  xcenter <- colSums(w * x)
  sqw <- sqrt(w)
  x <- sqw * sweep(x, 2, xcenter)
  ycenter <- colSums(w * y)
  y <- sqw * sweep(y, 2, ycenter)
  switch(match.arg(method), unbiased = crossprod(x, y) / (1 - sum(w ^ 2)),
    ML = crossprod(x, y))
}
