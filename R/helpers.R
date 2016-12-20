chunk <- function (x, chunks, each)
{
  n <- length(x)
  chunks <- if (missing(chunks))
    ceiling(n / each)
  else chunks
  ind <- sort(rep(1:chunks, length.out = n))
  lapply(as.list(unique(ind)), function (o) x[which(o == ind)])
}
