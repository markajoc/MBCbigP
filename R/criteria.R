#' @title Criteria for measuring the quality of clustering solution.
#'
#' @description Summarise a model-based clustering solution by a single number.
#'
#' @param mean A matrix, containing the mean vectors for the clusters.
#' @param sigma An array, containing the covariance matrices for the clusters.
#' @param groups The number of clusters.
#' @param ... Not used.
#'
#' @return A number.

bestseparated <-
function(mean, sigma, groups, ...)
{
  d <- dist(t(mean))
  s <- mean(apply(sigma, 3, diag))
  mean(head(sort(d / s, decreasing = TRUE), 5))
}
