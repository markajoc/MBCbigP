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
#' @param batchindex A list of column names or integer indices for the batches
#'   of \code{x}.
#' @param batchsize The desired number of columns in each batch.
#' @param maxiter The maximum number of iterations for the E-M algorithm.
#' @param plot Logical, should a visualisation be produced to show progress.
#' @param likelihood Logical, should the likelihood be calculated and used as a
#'   stopping criterion.
#' @param verbose Logical, should statements about progress be printed.
#' @param abstol Numeric, tolerance for log-likelihood stopping rule. Defaults
#'   to 1e-6.
#' @param method_sigma_AB One of \code{"analytic"} (default) or
#'   \code{"numerical"}, to choose the method of estimating the between-batch
#'   covariance for a cluster.
#' @param z Initial cluster membership probabilities. Numeric matrix.
#' @param scorefun A function to score the clustering of the data at each batch,
#'   defaults to \code{bestseparated}.
#' @param updateA If \code{FALSE} (default), at each batch the parameters for
#'   the current batch are maximised conditionally on the parameters from the
#'   previous batch. If \code{TRUE}, this is followed by also maximising the
#'   parameters from the previous batch conditionally on the current batch.
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
function (x, groups = 2, batches = NULL, batchindex = NULL, batchsize = NULL,
  maxiter = 50, plot = FALSE, likelihood = TRUE, verbose = TRUE, abstol = 1e-6,
  method_sigma_AB = c("analytic", "numerical"), z = NULL, scorefun = bestseparated,
  updateA = FALSE)
{
  ## Setting seed for development purposes. TODO: remove later.
  set.seed(13432523)

  ## Setting up some variables.

  starttime <- Sys.time()
  method_sigma_AB <- match.arg(method_sigma_AB)
  x <- data.matrix(x)
  N <- nrow(x)
  p <- ncol(x)

  if (verbose)
    cat("Setting up batches...\n")

  batchsize <- if (is.null(batchsize))
    min(ncol(x) / 2, 20)
  else batchsize

  if (!is.null(batchsize))
    batches <- round(p / batchsize)

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

  if (verbose)
    cat("\nBatch 1 \n-------------------------------------------")
  else {
    pb <- txtProgressBar(min = 0, max = Q, style = 3)
  }

  ## 'success' stores which batches defined by 'batchindex' succeed. 'score'
  ## stores the value from the function 'scorefun()'.

  success <- score <- rep(NA, Q)

  ## 'batch' is the list for storing the results from 'mbc()' for the first
  ## batch or 'mbc_cond()' for each subsequent batch.

  batch <- list()
  batch[[1L]] <- mbc(x = x[, batchindex[[1L]]], groups = groups, maxiter =
    maxiter, likelihood = likelihood, verbose = verbose, plot = plot)
  success[1L] <- 1
  score[1L] <- scorefun(mean = batch[[1]]$mean, sigma = batch[[1]]$sigma,
    groups = groups)

  if (!verbose)
    setTxtProgressBar(pb, 1)

  ## Other batches

  qprev <- 1  ## indexes 'batchindex', becomes 'q' on successful iteration.
  q <-        ## indexes 'batchindex', advances on each iteration.
  qindex <- 2 ## indexes 'batch', advances on each successful iteration.

  while (q <= Q){

    if (verbose)
      cat("\nBatch", q,"\n-------------------------------------------\n")

    tryobj <- try({

      ## Initial value for cluster membership probabilities for new batch.

      newz <- batch[[qindex - 1L]]$z
      #newz <- weightedmean.list(lapply(batch, `[[`, "z"), score[1:(q - 1)])
      #newz <- initialise.memberships(x = x[, batchindex[[q]]], groups = groups)

      ## Conditional model-based clustering.

      batch[[qindex]] <- mbc_cond(x_A = x[, batchindex[[qprev]]],
        x_B = x[, batchindex[[q]]], mean_A = batch[[qindex - 1L]]$mean,
        sigma_AA = batch[[qindex - 1L]]$sigma, z = newz, pro =
        batch[[qindex - 1L]]$pro, groups = groups, maxiter = maxiter, likelihood
        = likelihood, verbose = verbose, plot = plot, abstol = abstol,
        method_sigma_AB = method_sigma_AB, updateA = updateA)
    }, silent = FALSE)

    ## If 'mbc_cond()' throws an error, skip this batch and continue with the
    ## next batch. If not, calculate the 'scorefun()', update indices and move
    ## on.

    if (identical(class(tryobj), "try-error")){
      warning("batch ", q, " failed", immediate. = TRUE)
    } else {
      score[qindex] <- scorefun(mean = batch[[qindex]]$mean, sigma =
        batch[[qindex]]$sigma, groups = groups)
      success[qindex] <- q
      qprev <- q
      qindex <- qindex + 1L
    }
    if (!verbose)
      setTxtProgressBar(pb, q)
    q <- q + 1L
  }

  ## Extract cluster membership probabilities.

  z <- lapply(batch, function (x) x$z)

  ## Create a simple mean of membership probabilities.

  zave <- Reduce(`+`, z) / length(z)

  ## Make a list of colnames for each batch if 'batchindex' is not character.

  batchindexnames <- if (!is.character(batchindex[[1]]))
    lapply(batchindex, function (tmp) colnames(x)[tmp])
  else batchindex  

  ## Print the proportion of batches that succeeded.

  if (verbose){
    cat(paste0("\n\nFinished. Processed ", length(batch), "/", Q, " batches (",
      round(100 * length(batch) / Q),"%).\n\n"))
  } else close(pb)

  endtime <- Sys.time()

  ## Invisibly return S3 object of class 'mbcbigp'. 'z' is the final cluster
  ## membership probabilities, and 'zave' is a simple average of all membership
  ## probabilities over all successful batches. These values should be used with
  ## care. The membership probabilities for each batch are stored as 'z' in each
  ## element of the 'batch' list.

  invisible(structure(list(z = tail(z, 1)[[1]], zave = zave, batch = batch,
    batchindex = batchindexnames, success = success[!is.na(success)], score =
    score[!is.na(score)], groups = groups, time = endtime - starttime),
    class = "mbcbigp"))
}
