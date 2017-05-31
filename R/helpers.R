chunk <-
function (x, chunks, each)
{
  n <- length(x)
  chunks <- if (missing(chunks))
    ceiling(n / each)
  else chunks
  ind <- sort(rep(1:chunks, length.out = n))
  lapply(as.list(unique(ind)), function (o) x[which(o == ind)])
}

cov.test <-
function (z, ya, mua, yb, mub)
{
  ya <- data.matrix(ya)
  yb <- data.matrix(yb)
  tmp <- 0
  for (i in seq_along(z)){
    tmp <- tmp + z[i] * (ya[i, ] - mua) %*% t(yb[i, ] - mub)
  }
  tmp / sum(z)
}

reform_sigma <-
function (sigma_AA, sigma_AB, sigma_BB, groups){
  sigma <- array(dim = c(rep(dim(sigma_AA)[1] + dim(sigma_BB)[1], 2), groups))
  for (k in 1:groups){
    sigma[, , k] <- rbind(cbind(sigma_AA[, , k], sigma_AB[, , k]), cbind(t(
      sigma_AB[, , k]), sigma_BB[, , k]))
  }
  dimnames(sigma) <- list(c(dimnames(sigma_AA)[[1]], dimnames(sigma_BB)[[1]]),
    c(dimnames(sigma_AA)[[1]], dimnames(sigma_BB)[[1]]), 1:groups)
  sigma
}

reform_mean <-
function (mean_A, mean_B, groups){
  mean <- array(dim = c(dim(mean_A)[1] + dim(mean_B)[1], groups))
  for (k in 1:groups){
    mean[, k] <- c(mean_A[, k], mean_B[, k])
  }
  rownames(mean) <- c(rownames(mean_A), rownames(mean_B))
  mean
}

inrange <-
function (x, range)
{
  sapply(findInterval(x, range(range), rightmost.closed = TRUE), identical, 1L)
}

weightedmean.list <-
function(x, w, ...)
{
  if (!identical(length(x), length(w)))
    warning("lengths of x and w do not match")
  w <- unlist(w)
  w <- rep(w, length.out = length(x))
  if (any(w < 0))
    stop("negative weights in w")
  w <- w / sum(w)
  out <- x[[1L]] * 0
  for (i in seq_along(x)){
    out <- out + w[i] * x[[i]]
  }
  out
}

randIndex.mbcbigp <-
function(object, true, ...)
{
  clustering <- lapply(object$batch, function(obj) mclust::map(obj$z))
  unname(unlist(lapply(clustering, function(obj) flexclust::randIndex(obj, true)
    )))
}
