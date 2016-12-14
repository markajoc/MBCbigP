context("conditional functions")

library(mclust)
data(banknote)

o <- Mclust(banknote[, -1], 2, "VVV")

mean <- o$param$mean[, 1]
sigma <- o$param$variance$sigma[, , 1]

batcha <- 1:4

meana <- mean[batcha]
meanb <- mean[-batcha]

sigmaa <- sigma[batcha, batcha]
sigmab <- sigma[-batcha, -batcha]

covab <- sigma[batcha, -batcha]

mean_cond <- mean_conditional(mean = meanb, cov = covab, x_cond = banknote[, -1]
  [, batcha], mean_cond = meana, precision_cond = solve(sigmaa))

sigma_cond <- sigma_conditional(sigma = sigmab, cov = covab, precision_cond =
  solve(sigmaa))

test1 <- mvtnorm::dmvnorm(x = banknote[, -1], mean = mean, sigma = sigma, log =
  TRUE)

test2 <- 0 * test1

for (i in 1:nrow(banknote)){
test2[i] <- mvtnorm::dmvnorm(x = banknote[i, -1][, -batcha], mean = mean_cond[i,
  ], sigma = sigma_cond, log = TRUE) + mvtnorm::dmvnorm(x = banknote[i, -1][,
  batcha], mean = meana, sigma = sigmaa, log = TRUE)
}

test_that("conditional functions work properly", {
  expect_equal(test1, test2)
})
