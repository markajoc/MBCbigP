chunk <-
function (x, chunks, each)
{
  n <- length(x)
  chunks <- if (missing(chunks))
    ceiling(n / each)
  else chunks
  ind <- sort(rep(1:chunks, length.out = n))
  lapply(as.list(unique(ind)), function (o) x[which(o == ind)])
}

cov.test <-
function (z, ya, mua, yb, mub)
{
  ya <- data.matrix(ya)
  yb <- data.matrix(yb)
  tmp <- 0
  for (i in seq_along(z)){
    tmp <- tmp + z[i] * (ya[i, ] - mua) %*% t(yb[i, ] - mub)
  }
  tmp / sum(z)
}
