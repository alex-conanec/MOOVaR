#test q = 1
library(tidyverse)

set.seed(258164)
n <- 300
X <- runif(n = n, min = -2, max = 2)
Y <- X^2 + rnorm(n)

#parametres
h_n <- 0.2

B = 100
x_range <- runif(n=B, min = -2, max = 2)
y_range <- runif(n=B, min = -2, max = 6)

parallel::mclapply(x_range, function(x){
  sapply(y_range, function(y){
    abs( F_n(X = X, Y = Y, x = x, y = y, h_n = h_n) - pnorm(y, mean = x^2, sd = 1) )
  }) %>% mean()
}, mc.cores = 7) %>% unlist() %>% mean() -> a

library(testthat)

expect_that(a, is_less_than(0.1))

#test q = 2
n <- 300
X <- data.frame(X1 = runif(n = n, min = -2, max = 2),
                X2 = runif(n = n, min = -2, max = 2))

Y <- X[,1]^2 + X[,2]^2 + rnorm(n)


# z=sapply(t(d), K)
# length(z)
# # surface <- res %>% select(-y_barre) %>% spread(key = x2, value = q_n) %>% .[, -1] %>% as.matrix()
# plot3D::scatter3D(x = X[,1], y = X[,2], z = z)
# # ,
# #                   surf = list(x = x_range, y = x_range, z = surface))
# plot3Drgl::plotrgl()
#




h_n <- c(0.2, 0.2)

B = 10
x_range <- runif(n=B, min = -2, max = 2)
y_range <- runif(n=B, min = -2, max = 6)

parallel::mclapply(x_range, function(x){
  sapply(y_range, function(y){
    abs( F_n(X = X, Y = Y,
             x = data.frame(X1 = x, X2 = x),
             y = y, h_n = h_n)
         - pnorm(y, mean = 2*x^2, sd = 1) )
  }) %>% mean()
}, mc.cores = 7) %>% unlist() %>% mean() -> b

expect_that(b, is_less_than(0.1))
expect_lt(object = b, expected = 0.1) #similaire
