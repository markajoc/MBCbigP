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
#'
#' @return A list containing the estimated parameters for the mixture
#'   distribution.
#'
#' @examples
#' library(mclust)
#' data(banknote)
#' mbcbigp(x = banknote[, -1], groups = 2, batches = 2)
#' mbcbigp(x = banknote[, -1], groups = 2, batches = 3)

mbcbigp <-
function (x, groups = 2, batches = 3, maxiter = 200)
{
  x <- data.matrix(x)
  N <- nrow(x)
  p <- ncol(x)

  if ((lb <- length(batches)) < 2){
    Q <- batches
    batches <- sort(rep(1:batches, length.out = ncol(x)))
  } else {
    stopifnot(identical(lb, p))
    Q <- length(unique(batches))
  }

  ## Initialise membership probabilities

  z <- initialise.memberships(x, groups)

  ## Do the first batch using marginal density.

  batchindex <- which(batches == 1L)
  cat("\n")

  for (times in 1:maxiter){

    cat("Batch 1 , iteration ", times, "\n")

    ## Maximise

    parameters <- mstep(x[, batchindex, drop = FALSE], z, groups, p = length(
      batchindex))

    ## Expectation.

    z <- estep(x[, batchindex, drop = FALSE], parameters)

    plot(z[, 1L], ylim = 0:1, main = "q = 1")
    Sys.sleep(0.1)
  }

  ## Other batches

  for (q in 2:Q){

    ## Set up batch indices

    prevbatchindex <- batchindex
    batchindex <- which(batches == q)
    batchsize <- length(batchindex)

    ## Rename parameters from previous batch

    parameters_old <- parameters
    cat("\n")

    ## Start iterations.

    for (times in 1:maxiter){
      cat("Batch", q, ", iteration ", times, "\n")

      ## Maximisation.

      parameters <- mstep_cond(
        x1 = x[, batchindex, drop = FALSE],
        x2 = x[, prevbatchindex, drop = FALSE],
        z = z,
        mu2 = parameters_old$mean,
        sigma2 = parameters_old$sigma,
        groups = groups,
        p = NULL)

      ## Expectation.

      z <- estep_cond(
        x1 = x[, batchindex, drop = FALSE],
        x2 = x[, prevbatchindex, drop = FALSE],
        parameters1 = parameters,
        parameters2 = parameters_old)

      plot(z[, 1L], ylim = 0:1, main = paste0("q = ", q))
      Sys.sleep(0.1)
    }
  }
  invisible(structure(list(pro = colMeans(z), z = z), class = "mbc"))
}
