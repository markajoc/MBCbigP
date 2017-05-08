## Criteria for measuring the 'quality' of clustering solutions.

criteria_1 <-
function(mean, sigma, groups)
{
  d <- dist(t(mean))
  s <- mean(apply(sigma, 3, diag))
  max(d / s)
}
