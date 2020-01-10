h_n_opt <- function(X, Y, N = NROW(X)){

  Q <- function(tt, h_n){

    F_n_t <- function(x, y){
      x_mat <- matrix(rep(unlist(x), NROW(X) - 1), ncol = NCOL(x), byrow = TRUE)
      d <- (x_mat - X[-tt, ])/h_n
      Nn <- sum(sapply(t(d), K) * as.numeric(Y[-tt] <= y)) #eventuellement un pb
      Dn <- sum(sapply(t(d), K))
      Nn/Dn
    }

    A <- function(y, x){
      (as.numeric(Y[tt] <= y) - F_n_t(x, y) )^2 #* omega(y)
    }

    library(dplyr)
    lapply(data.frame(X), function(Xi){
      data.frame(X=runif(n = N, min = min(Xi), max = max(Xi)))
    }) %>% bind_cols -> x

    parallel::mclapply(seq_len(NROW(x)), function(i){
      formals(A)$x <- x[i,]
      integrate(A, lower = -Inf, upper = +Inf)$value
    }, mc.cores = 8) %>% unlist() %>% mean()

  }

  J <- function(h_n){
    sapply(1:n, function(tt){
      Q(tt, h_n)
    }) %>% sum()
  }

  # min J
  X <- data.frame(X)
  h_n_max <- #max de h_n est 1/4 du plus grand interval des X
    max((apply(data.frame(X), MARGIN = 2, FUN = max) -
           apply(data.frame(X), MARGIN = 2, FUN = min))/4)

  n <- NROW(X)
  optimize(J, c(10^(-9), h_n_max), maximum = FALSE)$minimum

  # a voir si on ne peut pas optimiser sur un vecteur h_n si l'amplitude des X est differentes

}
