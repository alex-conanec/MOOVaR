library(tidyverse)
set.seed(1234)
n <- 300
q = 3
mini <- rep(x = -2, times = q)
maxi <- rep(x = 2, times = q)

X <- lapply(seq_len(q), function(k){
  res <- data.frame(X=runif(n = n, min = mini[k], max = maxi[k]))
  names(res) <- paste0("X", k)
  res
}) %>% bind_cols()

fn <- list(
  Y1 = function(X) unlist(X[1] + 2*X[2] + X[3]),
  Y2 = function(X) unlist(-X[1] - X[2] - 3*X[3])
)

set.seed(123)
sensibilite_nsga(X, fn)

set.seed(1234)
sensibilite_nsga(X, fn)


