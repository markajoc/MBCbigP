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

    parameters <- mstep(x = x, z = z, groups = groups, p = p)

    if (plot)
      plotobj <- update.mbcplot(plotobj, parameters)

    ## Calculate log-likelihood.

    if (!is.null(loglik))
      loglikprevious <- loglik
    if (likelihood && (times %% 5 == 0)){
      loglik <- calcloglik(x = x, pro = parameters$pro, mean = parameters$mean,
        sigma = parameters$sigma, groups = parameters$groups)
    }

    ## If we have a previous log-likelihood, check that we are not decreasing
    ## the log-likelihood. If we are not increasing it by much, break the loop.

    if (!is.null(loglikprevious)){
      stopifnot(loglik >= loglikprevious)
      if (!(loglik > loglikprevious + (abs(loglikprevious) * 1e-7)))
        break
    }

    ## Calculate expected z values

    z <- estep(x = x, pro = parameters$pro, mean = parameters$mean, sigma =
      parameters$sigma, groups = parameters$groups)

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

mbc_cond <- function(x_A, x_B, mean_A, sigma_AA, z, pro, groups, maxiter = 500,
  likelihood = TRUE, verbose = FALSE, plot = FALSE, abstol = 1e-3)
{
  loglik <- vector()

  for (times in 1:maxiter){

    if (verbose)
      cat("Iteration ", times, "\n")

    parameters <- mstep_cond(x_B = x_B, x_A = x_A, z = z, mean_A = mean_A,
      sigma_AA = sigma_AA, sigma_AB = if (times == 1) NULL else parameters$cov,
      pro = if (times == 1) pro else parameters$pro, groups = groups)

    if (likelihood){
      loglik[times] <- calcloglik_split(x_A = x_A, x_B = x_B, pro =
        parameters$pro, mean_A = mean_A, mean_B = parameters$mean,
        sigma_AA = sigma_AA, sigma_AB = parameters$cov, sigma_BB =
        parameters$sigma, groups = groups)
    }

    z <- estep_cond(x_B = x_B, x_A = x_A, pro = parameters$pro, mean_A = mean_A,
      mean_B = parameters$mean, sigma_AA = sigma_AA, sigma_AB = parameters$cov,
      sigma_BB = parameters$sigma, groups = groups)

    if (likelihood & (times > 1)){
      if (abs(diff(c(loglik[times], loglik[times - 1]))) < abstol){
        print("no change in loglik")
        break
      }
    }
  }
  invisible(structure(list(parameters = parameters, z = z, loglik = loglik),
    class = c("mbc_cond", "mbc")))
}
