test_that("everything work", {

  library(tidyverse)
  library(MOOVaR)
  set.seed(1234)

  n <- 50
  p = 2

  # X quanti ----
  d1 = 3
  mini <- rep(x = -2, times = d1)
  maxi <- rep(x = 2, times = d1)
  X1 <- sapply(seq_len(d1), function(k){
    runif(n = n, min = mini[k], max = maxi[k])
  }) %>% as.data.frame()
  colnames(X1) = paste0("X", seq_len(d1))

  #X quali----
  d2 = 2
  lev = sample(x = 2:4, d2, replace = TRUE)
  X2 <- sapply(seq_len(d2), function(k){
    as.factor(sample(x = seq_len(lev[k]), size = n, replace = TRUE))
  }) %>% as.data.frame()
  colnames(X2) = paste0("X", (d1+1):(d1+d2))

  X = cbind(X1, X2)

  beta = matrix(runif(n = p*(d1 + sum(apply(X2, 2, function(x) nlevels(as.factor(x))))),
                      -1, 1), ncol = p)


  fn = function(X, beta){
    is_fac = sapply(X, is.factor)

    if (any(is_fac)){
      X_dummy = lapply(which(is_fac), function(j){
        model.matrix(~. + 0, data = X[,j, drop = FALSE])
      }) %>% bind_cols()

      X = bind_cols(X[!is_fac], X_dummy)
    }

    as.matrix(X) %*% beta
  }

  formals(fn)$beta = beta

  expect_success({
    res_nsga = NSGA(X = X,
                    fn = fn,
                    n_objective = p,
                    sens = rep("min", p),
                    k_tournament = 4,
                    n_echange = 3,
                    n_bloc = 2,
                    crossing_over_method = "bloc",
                    N = NROW(X),
                    g = NULL,
                    freq_m=0.1,
                    seed = NULL,
                    TT=30,
                    verbose = TRUE,
                    front_tracking = TRUE,
                    mutation_method = "mixte",
                    type_var = NULL
    )
  })

  expect_success( plot(res_nsga) )
  expect_success( plot(res_nsga, choice = "qlt") )
  expect_success( plot(res_nsga, choice = "evol") )

  # n_objective = p
  # sens = rep("min", n_objective)
  # k_tournament = 4
  # n_echange = 3
  # n_bloc = 2
  # crossing_over_method = "uniforme"
  # N = NROW(X)
  # g = NULL
  # # crossing_over_size = round(NCOL(X)/2)
  # freq_m=0.2
  # seed = NULL
  # B=10
  # verbose = TRUE
  # front_tracking = TRUE
  # mutation_method = "mixte"
  # type_var = NULL
  # b=1
})
