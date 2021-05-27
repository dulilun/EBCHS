test_that("EB_CS returns three vectors", {
  p = 1000
  k = 7
  # the prior distribution for lambda
  alpha = 2
  beta =  10
  # lambda
  lambda = rep(0, p)
  pi_0 = 0.9
  p_0 = floor(p*pi_0)
  p_1 = p-p_0
  lambda[(p_0+1):p] = rgamma(p_1, shape = alpha, rate=1/beta)
  # Generate a Poisson RV
  J = sapply(1:p, function(x){rpois(1, lambda[x]/2)})
  X = sapply(1:p, function(x){rchisq(1, k+2*J[x])})
  qq_set = seq(0.1, 0.9, 0.1)
  out = EB_CS(X, k, qq=qq_set, method='LS', mixture = FALSE)
  E = out$E_lambda
  V = out$V_lambda
  S = out$S_lambda
  expect_equal(length(E), p)
  expect_equal(length(V), p)
  expect_equal(length(S), p)
  expect_equal(sum(V>=0), p)
})
