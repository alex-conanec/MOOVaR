test_that("the stability of hopt and if is it closed to the real best h_n to estimate", {

  library(tidyverse)

  set.seed(258164)
  n <- 300
  X <- runif(n = n, min = -2, max = 2)
  Y <- X^2 + rnorm(n)
  X <- as.data.frame(X)

  critere_hopt_true <- function(h){
    sapply(seq_len(NROW(X)), function(i){
      (pnorm(q = Y[i], mean = X[i,]^2, sd = 1) -
         F_n(X = X, Y = Y, y=Y[i], x=X[i,], h_n = h))^2
    }) %>% sum()
  }

  h_true <- optimize(f = critere_hopt_true, interval = c(0.01, 1))$minimum

  h_n <- sapply(seq_len(10), function(i){
    hopt(X=X, Y = Y, n_tilde = 100)$h
  })


  h_n <- unlist(h_n[1,])
  #stabilite
  library(testthat)
  expect_lt(max(h_n) - min(h_n), 0.05)

  #compare to true hopt
  expect_lt(max(h_n) - h_true, 0.05)
  expect_lt(h_true - min(h_n), 0.05)
})


test_that("the stability of hopt and if is it closed to the real best h_n to estimate", {

  library(tidyverse)
  library(condQuant)

  set.seed(258164)
  n=300
  X <- data.frame(X1 = runif(n = n, min = -2, max = 2),
                  X2 = runif(n = n, min = -2, max = 2))

  Y <- X[,1]^2 + X[,2]^2 + rnorm(n)

  critere_hopt_true <- function(h){
    parallel::mclapply(seq_len(NROW(X)), function(i){
      (pnorm(q = Y[i], mean = X[i,1]^2 + X[i,2]^2, sd = 1) -
         F_n(X = X, Y = Y, y=Y[i], x=X[i,], h_n = h))^2
    }, mc.cores = 8) %>% unlist() %>% sum()
  }

  system.time(
    h_true <- optim(par = c(0.2, 0.2), fn = critere_hopt_true)$par
    )

  system.time(
  h_n <- sapply(seq_len(10), function(i){
    hopt(X=X, Y = Y, n_tilde = 100)$h
  })
  )

  ##
  hopt(X=X, Y = Y, n_tilde = 30)
  ###

  #stabilite
  library(testthat)
  for (j in 1:NCOL(X)){
    expect_lt(max(h_n[j]) - min(h_n[j]), 0.05)

    #compare to true hopt
    expect_lt(max(h_n[j]) - h_true[j], 0.05)
    expect_lt(h_true[j] - min(h_n[j]), 0.05)
  }

})
