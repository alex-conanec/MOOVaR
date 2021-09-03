context("def_cstr_X_space")
library(tidyverse)
library(MOOVaR)

test_that("verifie que l'aire de l'espace contraint est coherent avec une distribution uniform et normal", {


  #generate data ---
  #parameter of the data
  n = 1000
  d = 2 #dim de X

  mini = rep(0, d)
  maxi = rep(10, d)


  #runif
  b = matrix(c(0, 0, 5, 5), byrow = T, ncol = d)
  B = matrix(c(4, 4, 10, 10), byrow = T, ncol = d)
  p = NROW(B) #nombre de "foyer"
  n_MC = 5000
  A_MC = prod(maxi - mini)
  X_MC = sapply(seq_len(d), function(i){
    runif(n_MC, mini[i], maxi[i])
  })
  A_true_uni = 0.95*4*4+5*5

  S = 12 #nb de scenario
  res = parallel::mclapply(seq_len(S), function(s){
    X_U = sapply(seq_len(d), function(i){
      sapply(seq_len(p), function(k){
        runif(n/p, b[k,i], B[k,i])
      })
    })

    res = def_cstr_X_space(X=X_U, alpha = 0.05)
    A_MC*sum(res$g(x=X_MC))/n_MC
  }, mc.cores = 3) %>% unlist()

  expect_lt(mean(res), A_true_uni + 2)
  expect_gt(mean(res), A_true_uni - 2)



  #rnorm
  mu = matrix(c(2, 2, 7, 7), byrow = T, ncol = d)
  sig = matrix(c(1, 1, 1, 1), byrow = T, ncol = d)

  A_true_norm = 2 * pi * qnorm(0.975, sd = 1)^2

  S = 12 #nb de scenario
  res = parallel::mclapply(seq_len(S), function(s){
    X_N = sapply(seq_len(d), function(i){
      sapply(seq_len(p), function(k){
        rnorm(n/p, mu[k,i], sig[k,i])
      })
    })

    res = def_cstr_X_space(X=X_N, alpha = 0.05)
    A_MC*sum(res$g(x=X_MC))/n_MC
  }, mc.cores = 3) %>% unlist()


  expect_lt(mean(res), A_true_norm + 2) #surestime
  expect_gt(mean(res), A_true_norm - 2)

})


test_that("it work well with quali", {

  n = 500
  d = 2 #dim de X

  mini = rep(0, d)
  maxi = rep(10, d)

  #MC
  n_MC = 5000
  A_MC = prod(maxi - mini)
  X_MC = sapply(seq_len(d), function(i){
    runif(n_MC, mini[i], maxi[i])
  })

  #runif
  b = matrix(c(0, 0, 5, 5), byrow = T, ncol = d)
  B = matrix(c(4, 4, 10, 10), byrow = T, ncol = d)
  p = NROW(B) #nombre de "foyer"

  X_U = sapply(seq_len(d), function(i){
    sapply(seq_len(p), function(k){
      runif(n/p, b[k,i], B[k,i])
    })
  })

  #d2 = 1
  d2 = 1
  lev = rep(4, d2)
  X2 <- sapply(seq_len(d2), function(k){
    as.factor(sample(x = LETTERS[seq_len(lev[k])], size = n, replace = TRUE))
  }) %>% as.data.frame()

  X = data.frame(X_U, X2)
  colnames(X) = paste0("X", 1:NCOL(X))

  res = def_cstr_X_space(X = X, alpha = 0.15)

  #test visuel
  moda = LETTERS[2] #changer jusqu'a 5
  X = data.frame(X_U, V1 = X2) %>% filter(V1 == moda)
  X_f = data.frame(X_MC, V1 = moda)
  plot(X_f[,1:2], col = c("red", "green")[as.numeric(res$g(X_f) + 1)] )
  points(X, col = "blue")


  #avec uniquement quali
  res = def_cstr_X_space(X = X2, alpha = 0.05)
  res$g(x = "A")
  res$g(x = "E")

  #d2 = 3
  d2 = 3
  lev = rep(2, d2)
  X2 <- sapply(seq_len(d2), function(k){
    as.factor(sample(x = LETTERS[seq_len(lev[k])], size = n, replace = TRUE))
  }) %>% as.data.frame()

  X = data.frame(X_U, X2)
  colnames(X) = paste0("X", 1:NCOL(X))
  res = def_cstr_X_space(X = X, alpha = 0.15)

  #test visuel
  moda = c("B", "B", "A")  #changer jusqu'a 5
  X = data.frame(X_U, X2) %>% filter(V1 == moda[1] & V2 == moda[2] & V3 == moda[3])
  X_f = data.frame(X_MC, V1 = moda[1], V2 = moda[2], V3 = moda[3])
  plot(X_f[,1:2], col = c("red", "green")[as.numeric(res$g(X_f) + 1)] )
  points(X, col = "blue")


  #avec uniquement quali
  res = def_cstr_X_space(X = X2, alpha = 0.05)
  res$g(x = c("A", "B", "A"))
  res$g(x = data.frame(V1="A", V2="B", V3="A"))
  res$g(x = c("A", "B", "C"))


})

test_that("test that the part of X class feasible is above alpha", {

  set.seed(123)
  n = 5000
  d = 3
  alpha = 0.05

  #runif
  b = matrix(c(0, 0, 0, 5, 5, 5), byrow = T, ncol = d)
  B = matrix(c(4, 4, 4, 10, 10, 10), byrow = T, ncol = d)
  p = NROW(B) #nombre de "foyer"

  X_U = sapply(seq_len(d), function(i){
    sapply(seq_len(p), function(k){
      runif(n/p, b[k,i], B[k,i])
    })
  })

  d2 = 3
  lev = rep(3, d2)
  X2 <- sapply(seq_len(d2), function(k){
    as.factor(sample(x = LETTERS[seq_len(lev[k])], size = n, replace = TRUE))
  }) %>% as.data.frame()

  X = cbind(X = X_U, X2)
  res = def_cstr_X_space(X, alpha = alpha)


  expect_gt(sum(res$g(X))/n, alpha)


})

