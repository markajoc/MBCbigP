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
function (x, groups = 2, batches = 3, batchindex = NULL, maxiter = 50, plot =
  FALSE, likelihood = FALSE, verbose = TRUE)
{
  set.seed(13432523)
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
  #if (plot)
  #  plotobject <- plotmbc(x, parameters = NULL, groups = groups)

  ## Do the first batch using marginal density.

  ## Initialise membership probabilities

  z <- initialise.memberships(x[, batchindex[[1L]], drop = FALSE], groups)

  ## Initialise NULL log-likelihoods

  loglikprevious <- NULL
  loglik <- NULL

  if (verbose)
    cat("\n")

  for (times in 1:maxiter){

    if (verbose)
      cat("Batch 1 , iteration ", times, "\n")

    ## Maximise

    parameters <- mstep(x[, batchindex[[1L]], drop = FALSE], z, groups, p =
      length(batchindex[[1L]]))

    #if (plot)
    #  plotobject <- update(plotobject, parameters)

    ## Expectation.

    z <- estep(x = x[, batchindex[[1L]], drop = FALSE], groups = groups, mean =
      parameters$mean, sigma = parameters$sigma, pro = parameters$pro)

    if (plot)
      plot(z[, 1L], ylim = 0:1, main = "q = 1")

    ## Calculate log-likelihood.

    if (!is.null(loglik))
      loglikprevious <- loglik
    if (likelihood && (times %% 5 == 0)){
      loglik <- calcloglik(x = x[, batchindex[[1L]]], pro = parameters$pro,
        mean = parameters$mean, sigma = parameters$sigma, groups = groups)
    }

    ## If we have a previous log-likelihood, check that we are not decreasing
    ## the log-likelihood. If we are not increasing it by much, break the loop.

    if (!is.null(loglikprevious)){
      stopifnot(loglik >= loglikprevious)
      if (!(loglik > loglikprevious + (abs(loglikprevious) * 1e-5)))
        break
    }
  }

  batch <- list()
  batch[[1]] <- parameters

  loglik1 <- list()
  loglik1[[1]] <- loglik
  loglik <- loglik1

  ## Other batches

  for (q in 2:Q){

    batch[[q]] <- list()
    batch[[q]]$partrace <- list()

    loglik[[q]] <- vector()

    ## Rename parameters from previous batch

    parameters_old <- parameters

    if (verbose)
      cat("\n")

    ## Start iterations.

    for (times in 1:maxiter){

      if (verbose)
        cat("Batch", q, ", iteration ", times, "\n")

      ## Maximisation.

      parameters <- mstep_cond(
        x_B = x[, batchindex[[q]], drop = FALSE],
        x_A = x[, batchindex[[q - 1L]], drop = FALSE],
        z = z,
        mean_A = parameters_old$mean,
        sigma_AA = parameters_old$sigma,
        sigma_AB = parameters$cov,
        pro = parameters$pro,
        groups = groups)

      if (likelihood){
        loglik[[q]][times] <- calcloglik_split(
          x_A = x[, batchindex[[q - 1L]], drop = FALSE],
          x_B = x[, batchindex[[q]], drop = FALSE],
          pro = parameters$pro,
          mean_A = parameters_old$mean,
          mean_B = parameters$mean,
          sigma_AA = parameters_old$sigma,
          sigma_AB = parameters$cov,
          sigma_BB = parameters$sigma,
          groups = groups)
      }

      #batch[[q]]$partrace[[times]] <- parameters


      #if (plot)
      #  plotobject <- update(plotobject, parameters)

      ## Expectation.

      z <- estep_cond(
        x_B = x[, batchindex[[q]], drop = FALSE],
        x_A = x[, batchindex[[q - 1L]], drop = FALSE],
        pro = parameters$pro,
        mean_A = parameters_old$mean,
        mean_B = parameters$mean,
        sigma_AA = parameters_old$sigma,
        sigma_AB = parameters$cov,
        sigma_BB = parameters$sigma,
        groups = groups)

      #if (plot)
      #  plot(z[, 1L], ylim = 0:1, main = paste0("q = ", q))
      if (likelihood & (times > 1)){
        if (abs(diff(c(loglik[[q]][times], loglik[[q]][times - 1]))) < 1e-3){
          print("no change in loglik")
          break
        }
      }
    }
    #batch[[q]] <- parameters

  }
  invisible(structure(list(pro = colMeans(z), z = z, batch = batch, loglik =
    loglik), class = "mbc"))
}
