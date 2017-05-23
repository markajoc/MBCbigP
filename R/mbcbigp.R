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
function (x, groups = 2, batches = 3, batchindex = NULL, batchsize = NULL,
  maxiter = 50, plot = FALSE, likelihood = FALSE, verbose = TRUE, abstol = 1e-3,
  method_sigma_AB = c("numeric", "analytic"), z = NULL, scorefun = criteria_1)
{
  method_sigma_AB <- match.arg(method_sigma_AB)
  set.seed(13432523)
  x <- data.matrix(x)
  N <- nrow(x)
  p <- ncol(x)

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

  #if (plot){
  #  plotobject <- plotmbc(x[, c(batchindex[[1L]], batchindex[[2L]])],
  #     parameters = NULL, groups = groups)
  #}

  ## Do the first batch using regular model-based clustering.

  if (verbose)
    cat("\nBatch 1 \n-------------------------------------------")
  else {
    pb <- txtProgressBar(min = 0, max = Q, style = 3)
  }

  success <- score <- rep(NA, Q)

  batch <- list()
  batch[[1L]] <- mbc(x = x[, batchindex[[1L]]], groups = groups, maxiter =
    maxiter, likelihood = likelihood, verbose = verbose, plot = plot)
  success[1L] <- 1

  score[1L] <- scorefun(mean = batch[[1]]$mean, sigma = batch[[1]]$sigma,
    groups = groups)
  setTxtProgressBar(pb, 1)

  ## Other batches

  qprev <- 1 ## indexes 'batchindex', becomes 'q' on successful iteration.
  q <- ## indexes 'batchindex', advances on each iteration.
  qindex <- 2 ## indexes 'batch', advances on each successful iteration.

  while (q <= Q){

    if (verbose)
      cat("\nBatch", q,"\n-------------------------------------------\n")

    tryobj <- try({
      newz <- batch[[qindex - 1L]]$z
      #newz <- initialise.memberships(x = x[, batchindex[[q]]], groups = groups)
      batch[[qindex]] <- mbc_cond(x_A = x[, batchindex[[qprev]]],
        x_B = x[, batchindex[[q]]], mean_A = batch[[qindex - 1L]]$mean,
        sigma_AA = batch[[qindex - 1L]]$sigma, z = newz, pro =
        batch[[qindex - 1L]]$pro, groups = groups, maxiter = maxiter, likelihood
        = likelihood, verbose = verbose, plot = plot, abstol = abstol,
        method_sigma_AB = method_sigma_AB)
    }, silent = FALSE)
    if (identical(class(tryobj), "try-error")){
      warning("batch ", q, " failed")
    } else {
      score[qindex] <- scorefun(mean = batch[[qindex]]$mean, sigma =
        batch[[qindex]]$sigma, groups = groups)
      success[qindex] <- q
      qprev <- q
      qindex <- qindex + 1L
    }
    setTxtProgressBar(pb, q)
    q <- q + 1L
  }
  z <- lapply(batch, function (x) x$z)
  zave <- Reduce(`+`, z) / length(z)
  batchindexnames <- if (!is.character(batchindex[[1]]))
    lapply(batchindex, function (tmp) colnames(x)[tmp])
  if (verbose){
    cat(paste0("\n\nFinished. Processed ", length(batch), "/", Q, " batches (",
      round(100 * length(batch) / Q),"%)\n\n"))
  } else close(pb)
  invisible(structure(list(z = zave, batch = batch, batchindex =
    batchindexnames, success = success[!is.na(success)], score =
    score[!is.na(score)], groups = groups), class = "mbcbigp"))
}
