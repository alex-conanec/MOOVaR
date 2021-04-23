library(testthat)
test_that("unscale correctly", {

  n <- 300
  d = 4
  X <- as.data.frame(matrix(runif(n = n*d, min = -2, max = 2), nrow = n, ncol = d))

  X_sc = scale(X)
  X_sc_param = list(mu = attr(X_sc, "scaled:center"),
                          sd = attr(X_sc, "scaled:scale"))

  expect_true(all(scale_x(X_sc, mu = X_sc_param$mu, s = X_sc_param$sd,
                          reverse = TRUE) - X < 10^(-15)))


})

test_that("scale correctly with parameter estimated on train data", {
  n <- 1000
  d = 4
  X <- as.data.frame(matrix(runif(n = n*d, min = -2, max = 2), nrow = n, ncol = d))

  train = sample(seq_len(n), 0.8*n)
  X_train = X[train, ]
  X_test = X[-train, ]

  X_train_sc = scale(X_train)
  X_test_sc = scale(X_test)

  X_test_sc_param = list(mu = attr(X_test_sc, "scaled:center"),
                         sd = attr(X_test_sc, "scaled:scale"))

  scale_x(X_test, mu = X_test_sc_param$mu, s = X_test_sc_param$sd,
          reverse = FALSE) - X_test_sc

  X_train_sc_param = list(mu = attr(X_train_sc, "scaled:center"),
                          sd = attr(X_train_sc, "scaled:scale"))

  scale_x(X_test, mu = X_train_sc_param$mu, s = X_train_sc_param$sd,
          reverse = FALSE) - X_test_sc


})
