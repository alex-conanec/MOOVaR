library(dplyr)
test_that("give the right quantile conditionnal estimation for dim(Y)=1", {

  n <- 300
  X <- runif(n = n, min = -2, max = 2)
  Y <- X^2 + rnorm(n)

  sapply(seq_len(100), function(i){
    alpha <- runif(1, min=0, max=1)
    x <- runif(n = 1, min = -2, max = 2)

    q_hat <- conditionnal_quantile(X, Y, x, alpha)
    q_true <- qnorm(p = alpha, mean = x^2, sd = 1)

    (q_hat - q_true)^2
  }) %>% mean() %>% sqrt()-> MSE

  MSE
  plot(X, Y)
  abline(v=x, col = "red")
  abline(h=q_hat, col = "blue")
  abline(h=q_true, col = "green")

  MSE_ref <- 0.5
  expect_true(MSE < MSE_ref)

})


test_that("the quantile estimation converge to the true quantile when n tend to +inf", {

  n_range <- c(50, 100, 150, 200, 300, 450, 600, 800, 1000, 1500, 2200, 3000)

  sapply(n_range, function(n){
    X <- runif(n = n, min = -2, max = 2)
    Y <- X^2 + rnorm(n)

    parallel::mclapply(seq_len(80), function(i){
      alpha <- 0.2
      x <- runif(n = 1, min = -2, max = 2)

      q_hat <- conditionnal_quantile(X, Y, x, alpha, h_n = 0.1)
      q_true <- qnorm(p = alpha, mean = x^2, sd = 1)

      (q_hat - q_true)^2
    }, mc.cores = 8) %>% unlist() %>% mean() %>% sqrt()

  }) -> MSE

  plot(x=n_range, y=MSE, type = 'l')

  expect_true(MSE < MSE_ref)

})


test_that("the quantile estimation converge to the true quantile when n tend to +inf", {

  n_range <- c(50, 100, 150, 200, 300, 450, 600, 800, 1000, 1500, 3000)

  sapply(n_range, function(n){
    X <- data.frame(x1 = runif(n = n, min = -2, max = 2),
                    x2 = runif(n = n, min = -2, max = 2))
    Y <- X$x1 + X$x2 + rnorm(n)

    parallel::mclapply(seq_len(80), function(i){
      alpha <- 0.1
      x <- runif(n = 1, min = -2, max = 2)

      q_hat <- conditionnal_quantile(X, Y, x, alpha)
      q_true <- qnorm(p = alpha, mean = x$x1 + x$x2, sd = 1)

      (q_hat - q_true)^2
    }, mc.cores = 8) %>% unlist() %>% mean() %>% sqrt()

  }) -> MSE

  plot(x=n_range, y=MSE, type = 'l')

  expect_true(MSE < MSE_ref)

})


test_that("give the right quantile conditionnal estimation for dim(Y)>1", {

  n=10000
  X <- data.frame(x1 = runif(n = n, min = -2, max = 2),
                  x2 = runif(n = n, min = -2, max = 2))
  Y <- X$x1 + X$x2 + rnorm(n)

  sapply(seq_len(100), function(i){
    alpha <- 0.1
    x <- c(x1 = runif(n = 1, min = -2, max = 2),
           x2 = runif(n = 1, min = -2, max = 2))

    q_hat <- conditionnal_quantile(X, Y, x, alpha)
    q_true <- qnorm(p = alpha, mean = x$x1 + x$x2, sd = 1)

    (q_hat - q_true)^2
  }) %>% mean() %>% sqrt()-> MSE

})
