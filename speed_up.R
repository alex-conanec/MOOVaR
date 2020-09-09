library(tidyverse)
library(condQuant)

set.seed(258164)
n <- 300
X <- runif(n = n, min = -2, max = 2)
Y <- X^2 + rnorm(n)

#parametres
alpha <- 0.1
h_n <- hopt(X, Y)$h

#estimation de q_n et calcul de q
x_range <- seq(-2, 2, length.out = n/2)

parallel::mclapply(x_range, function(x){
  q <- qnorm(p = alpha, mean = x^2, sd = 1)
  q_n <- conditionnal_quantile(X, Y, x, alpha, h_n = h_n)
  data.frame(q=q, q_n=q_n)
}, mc.cores = 8) %>% bind_rows()-> res

#plot
res <- cbind(x_range = x_range, res) %>%
  gather(key = key, value = value, -x_range)
data <- cbind(data.frame(X=X, Y=Y), res)

ggplot(data) +
  geom_point(aes(x=X, y=Y)) +
  geom_line(aes(x=x_range, y=value, col = key))


#q=2
X <- data.frame(X1 = runif(n = n, min = -2, max = 2),
                X2 = runif(n = n, min = -2, max = 2))
Y <- X[,1]^2 + X[,2]^2 + rnorm(n)

h_n = hopt(X, Y)$h

#test vitesse
library(microbenchmark)
library(profvis)

# conditionnal_quantile
microbenchmark(conditionnal_quantile(X, Y, x=c(-0.8, -0.8), alpha, h_n = h_n))
microbenchmark(qnorm(p = 0.1, mean = (-0.8)^2 + (-0.8)^2, sd = 1))

profvis::profvis(conditionnal_quantile(X, Y, x = c(-0.8, -0.8), alpha, h_n = h_n))

# hopt
# microbenchmark(hopt(X, Y))
profvis::profvis(hopt(X, Y))




