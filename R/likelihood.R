calcloglik <-
function (x, parameters)
{
  K <- parameters$groups
  alpha <- parameters$mixing
  mu <- parameters$mean
  sigma <- parameters$sigma
  tmp <- matrix(as.double(NA), nrow(x), K)
  for (k in 1:K){
    tmp[, k] <- alpha[k] * mvtnorm::dmvnorm(x, mu[, k], sigma[, , k])
  }
  sum(log(rowSums(tmp)))
}

logLik.mbc <-
function (object, ...)
{
  
}
