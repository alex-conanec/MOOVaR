set.seed(123)
library(MOOVaR)
library(tidyverse)

{# X
n <- 300
d = 4
X <- as.data.frame(matrix(runif(n = n*d, min = -2, max = 2), nrow = n, ncol = d))
colnames(X) = paste0("X", seq_len(d))

#calcul des Y selon une relation lineaire + du bruit
p = 3
beta = matrix(round(runif(n = p*d, min = -2, max = 2), 1), ncol = p, nrow = d)

Y = as.matrix(X) %*% beta
cor(Y)
epsilon1 = rnorm(n, sd = 0.5)
epsilon2 = scale(-0.5*epsilon1 + rnorm(n, sd = 0.25))
epsilon3 = scale(0.5*epsilon1 + rnorm(n, sd = 0.25))

epsilon = matrix(c(epsilon1, epsilon2, epsilon3), ncol = p, nrow = n, byrow = F)

Y = Y + epsilon}


#test
vizu = function(m){
  for (j in seq_len(p)){
    plot(Y[,j], Y[,j], col = "red", main = j)
    points(Y[,j], m(X)[,j])
  }
}

m = fit_model(X, Y, tau = 0.3, reg_method = "linear", method = "quantile")
vizu(m)

m = fit_model(X, Y, tau = 0.3, reg_method = "neural_network", method = "quantile",
                       penalty = rep(0, p))

vizu(m)

m = fit_model(X, Y, tau = 0.3, reg_method = "neural_network", method = "quantile",
                       penalty = NULL)
vizu(m)


m = fit_model(X, Y, reg_method = "linear", method = "expected")
vizu(m)

#error
m = fit_model(X, Y, reg_method = "neural_network", method = "expected")
vizu(m)
rm(list = ls())

# Sys.setenv(RETICULATE_PYTHON="/usr/local/share/.virtualenvs/r-reticulate/bin/")

library(reticulate)
library(keras)
inputs <- keras::layer_input(shape = d, name = 'aux_input')

outputs = inputs %>%
  keras::layer_dense(units = 4, activation = 'relu') %>%
  layer_dropout(rate = 0.5) %>%
  keras::layer_dense(units = 4, activation = "relu") %>%
  layer_dropout(rate = 0.5) %>%
  keras::layer_dense(units = 1)

model <- keras::keras_model(inputs = inputs, outputs = outputs)

model %>% keras::compile(
  optimizer = keras::optimizer_rmsprop(),
  loss = "mse",
  metrics = c("mean_absolute_error")
)

x = scale(X)
y = scale(Y[,1])
model %>% keras::fit(x = list(x), y = list(y), batch_size  = n/10,
                epochs = 500, validation_split = 0.2, verbose = 0, shuffle=T)

plot(y, predict(model, x))
mm = lm(y~x)

plot(y, predict(mm, data = x))

