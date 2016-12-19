#' @title Fit a mixture of multivariate Gaussians
#'
#' @description Fit a mixture of multivariate Gaussians
#'
#' @param x Data frame or a matrix
#' @param groups The number of groups/mixture components to fit.
#' @param maxiter The maximum number of iterations for the E-M algorithm.
#' @param likelihood Logical indicating whether the log-likelihood should be
#'   calculated at each step and returned (defaults to \code{TRUE}).
#'
#' @return A list containing the estimated parameters for the mixture
#'   distribution.
#'
#' @examples
#' library(mclust)
#' data(banknote)
#' mbc(x = banknote[, -1], groups = 2)

mbc <-
function (x, groups = 2, maxiter = 500, likelihood = TRUE, verbose = FALSE, plot
  = FALSE)
{
  x <- data.matrix(x)
  N <- nrow(x)
  p <- ncol(x)

  ## Initialise z matrix.

  if (verbose)
    cat("\nInitialising clusters...")
  z <- initialise.memberships(x, groups)

  ## Initialise NULL log-likelihoods

  loglikprevious <- NULL
  loglik <- NULL

  if (verbose)
    cat("\nStarting E-M iterations...")

  if (plot)
    plotobj <- plotmbc(x = x, parameters = NULL, groups = groups, p = p)

  for (times in 1:maxiter){

    ## Calculate maximum likelihood estimates for the mixing proportions, means
    ## and covariance matrices

    parameters <- mstep(x, z, groups, p)

    if (plot)
      plotobj <- update.mbcplot(plotobj, parameters)

    ## Calculate log-likelihood.

    if (!is.null(loglik))
      loglikprevious <- loglik
    if (likelihood && (times %% 5 == 0)){
      loglik <- calcloglik(x, parameters)
    }

    ## If we have a previous log-likelihood, check that we are not decreasing
    ## the log-likelihood. If we are not increasing it by much, break the loop.

    if (!is.null(loglikprevious)){
      stopifnot(loglik >= loglikprevious)
      if (!(loglik > loglikprevious + (abs(loglikprevious) * 1e-7)))
        break
    }

    ## Calculate expected z values

    z <- estep(x, parameters)

    if (verbose && (times %% 5 == 0))
      cat("\nFinished iteration ", times)
  }
  if (verbose)
    cat("\n")
  parameters$z <- z
  parameters$loglik <- loglik

  if (plot)
    plotmbc(x = x, parameters = parameters, groups = groups, p = p)

  invisible(structure(parameters, class = "mbc"))
}
