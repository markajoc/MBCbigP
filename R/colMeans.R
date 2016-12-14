colMeans.weighted <-
function (x, w, ...)
{
  apply(x[, , drop = FALSE], 2, weighted.mean, w = w, ...)
}
