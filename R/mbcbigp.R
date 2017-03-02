#' @title Fit a mixture of multivariate Gaussians using sequential conditioning
#'
#' @description Fit a mixture of multivariate Gaussians
#'
#' @param x Data frame or a matrix
#' @param groups The number of groups/mixture components to fit.
#' @param batches The batches to split the columns of \code{x} for
#'   processing. A single number gives the desired number of batches, whereas a
#'   vector of length \code{ncol(x)} is treated as an integer index of batch
#'   membership.
#' @param maxiter The maximum number of iterations for the E-M algorithm.
#' @param plot Logical, should a visualisation be produced to show progress.
#' @param likelihood Logical, should the likelihood be calculated and used as a
#'   stopping criterion.
#' @param verbose Logical, should statements about progress be printed.
#'
#' @return A list containing the estimated parameters for the mixture
#'   distribution.
#'
#' @examples
#' \dontrun{
#' data(banknote, package = "mclust")
#' mbcbigp(x = banknote[, -1], groups = 2, batches = 2)
#' mbcbigp(x = banknote[, -1], groups = 2, batches = 3)
#' }

mbcbigp <-
function (x, groups = 2, batches = 3, batchindex = NULL, maxiter = 10, plot =
  FALSE, likelihood = FALSE, verbose = TRUE)
{
  x <- data.matrix(x)
  N <- nrow(x)
  p <- ncol(x)

  if (is.null(batchindex)){
    if ((lb <- length(batches)) < 2){
      Q <- batches
      batches <- sort(rep(1:batches, length.out = p))
    } else {
      stopifnot(identical(lb, p))
    }
    batchindex <- lapply(as.list(unique(batches)), function (o) which(o ==
      batches))
  }
  Q <- length(batchindex)
  if (plot)
    plotobject <- plotmbc(x, parameters = NULL, groups = groups)

  ## Do the first batch using marginal density.

  ## Initialise membership probabilities

  z <- initialise.memberships(x[, batchindex[[1L]], drop = FALSE], groups)

  ## Initialise NULL log-likelihoods

  loglikprevious <- NULL
  loglik <- NULL

  cat("\n")

  for (times in 1:maxiter){

    cat("Batch 1 , iteration ", times, "\n")

    ## Maximise

    parameters <- mstep(x[, batchindex[[1L]], drop = FALSE], z, groups, p =
      length(batchindex[[1L]]))

    if (plot)
      plotobject <- update(plotobject, parameters)

    ## Expectation.

    z <- estep(x[, batchindex[[1L]], drop = FALSE], parameters)

    #if (plot)
    #  plot(z[, 1L], ylim = 0:1, main = "q = 1")

    ## Calculate log-likelihood.

    if (!is.null(loglik))
      loglikprevious <- loglik
    if (likelihood && (times %% 5 == 0)){
      loglik <- calcloglik(x[, batchindex[[1L]]], parameters)
    }

    ## If we have a previous log-likelihood, check that we are not decreasing
    ## the log-likelihood. If we are not increasing it by much, break the loop.

    if (!is.null(loglikprevious)){
      stopifnot(loglik >= loglikprevious)
      if (!(loglik > loglikprevious + (abs(loglikprevious) * 1e-5)))
        break
    }
  }

  ## Other batches

  for (q in 2:Q){

    ## Rename parameters from previous batch

    parameters_old <- parameters
    cat("\n")

    ## Start iterations.

    for (times in 1:maxiter){
      cat("Batch", q, ", iteration ", times, "\n")

      ## Maximisation.

      parameters <- mstep_cond(
        x1 = x[, batchindex[[q]], drop = FALSE],
        x2 = x[, batchindex[[q - 1L]], drop = FALSE],
        z = z,
        mu2 = parameters_old$mean,
        sigma2 = parameters_old$sigma,
        groups = groups,
        p = NULL)

      if (plot)
        plotobject <- update(plotobject, parameters)

      ## Expectation.

      z <- estep_cond(
        x1 = x[, batchindex[[q]], drop = FALSE],
        x2 = x[, batchindex[[q - 1L]], drop = FALSE],
        parameters1 = parameters,
        parameters2 = parameters_old)

      #if (plot)
      #  plot(z[, 1L], ylim = 0:1, main = paste0("q = ", q))
    }
  }
  invisible(structure(list(pro = colMeans(z), z = z), class = "mbc"))
}
