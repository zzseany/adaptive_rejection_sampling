###########################################################
## STAT 243 Final Project
## Adaptive Rejection Sampling Algorithm
##
## Date: 12/12/18
## Authors: Jonathan Morrell, Ziyang Zhou, Vincent Zijlmans
###########################################################

library('ars')

test_that("check if the input density is function", {
  ## std. normal dist
  expect_equal(length(ars(dnor, 100)), 100)
  
  ## exp. dist
  expect_equal(length(ars(dexp, 100, c(0.5), c(0,Inf))), 100)
  
  ## gamma dist
  expect_equal(length(ars(dgamma, 100, x0=c(10), bounds=c(0.0, Inf), shape=3.0, rate=2.0)), 100)
  
  ## case when the input density is not a function
  expect_error(ars(1, 100))
})

test_that("check if the bound is of length 2", {
  ## length 2 case
  expect_equal(length(ars(dnorm, 100, bounds = c(-100,100))), 100)
  
  ## case when length not equal to 2 
  expect_error(ars(dnorm, 100, bounds = c(-1,0,1)))
})

test_that("check if the upper bound and the lower bound are equal", {
  expect_warning(ars(dnorm,100,bounds = c(100,100)))
})

test_that("check if x_0 is between bounds", {
  ## when x_0 is at the left side of bounds
  expect_error(ars(dnorm, 100, x0 = c(-1), bounds = c(0,1)))
  
  ## when x_0 is at the right side of bounds
  expect_error(ars(dnorm, 100, x0 = c(2), bounds = c(0,1)))
  
  ## when x_0 is in between bounds
  expect_equal(length(ars(dnorm, 100, x0 = c(0.5), bounds = c(0,1))), 100)
  
})

test_that("check if the sampling distribution is close enough to the original distribution", {
  # Note we perform Kolmogorov-Smirnov Test with the null
  # hypothesis that the samples are drawn from the same
  # continuous distribution
  
  ## std. normal 
  expect_equal(ks.test(ars(dnorm, 1000), rnorm(1000))$p.value > 0.05, T)
  
  ## exp(1)
  expect_equal(ks.test(ars(dexp, 1000, x0 = 5, bounds = c(0, Inf)), rexp(1000))$p.value > 0.05, T)
  
  ## gamma(3,2)
  expect_equal(ks.test(ars(dgamma, 1000, x0 = 5, bounds = c(0, Inf), shape = 3, scale = 2),
                       rgamma(1000, shape = 3, scale = 2))$p.value > 0.05, T)
  ## unif(0,1)
  expect_equal(ks.test(ars(dunif, 1000, x0 = 0.5, bounds = c(0,1)), runif(1000))$p.value > 0.05, T)
  
  ## logistics 
  expect_equal(ks.test(ars(dlogis, 1000, x0 = 0, bounds = c(-10,10)), rlogis(1000))$p.value > 0.05, T)
  
  ## beta(3,2)
  expect_equal(ks.test(ars(dbeta, 1000, x0 = 0.5, bounds = c(0, 1), shape1 = 3, shape2 = 2), 
                       rbeta(1000, shape1 = 3, shape2 = 2))$p.value > 0.05, T)
  
  ## laplace
  library(rmutil)
  expect_equal(ks.test(ars(dlaplace, 1000, x0 = 0, bounds = c(-5,5)),
                       rlaplace(1000))$p.value > 0.05, T)
  
  ## chi(2)
  expect_equal(ks.test(ars(dchisq, 1000, x0 = 1, bounds = c(0, Inf), df = 2),
                       rchisq(1000, df = 2))$p.value > 0.05, T)
  
  ## weibull(2,1)
  expect_equal(ks.test(ars(dweibull, 1000, shape = 2, x0 = 1, bounds = c(0, Inf)),
                       rweibull(1000, shape = 2))$p.value > 0.05, T)
})

test_that("check if non-concavity", {
  
  ## simple exponential exp(x^2)
  de = function(x){
    return(exp(x^2))
  }
  expect_error(ars(de, 1000, x0 = 0, bounds = c(-5,5)))
  
  ## student t(2)
  expect_error(ars(dt, 1000, x0 = 1, bounds = c(-5,5), df = 2))
  
  ## cauchy
  expect_error(ars(dcauchy, 1000, x0 = 0, bounds = c(-5,5)))
  
  ## pareto(1,2)
  expect_error(ars(dpareto, 1000, x0 = 3, bounds = c(1, Inf), m = 1, s = 2))
  
  ## lognormal
  expect_error(ars(dlnorm, 1000, x0 = 1, bounds = c(0, Inf)))
  
  ## F dist (1,1)
  expect_error(ars(stats::df, 1000, x0 = 1, bounds = c(0, Inf), df1 = 1, df2 = 2))
})


