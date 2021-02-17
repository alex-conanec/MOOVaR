library(tidyverse)
library(optisure)
set.seed(12)
n <- 300

# X quanti ----
# X
d1 = 3
mini <- rep(x = -2, times = d1)
maxi <- rep(x = 2, times = d1)
X1 <- sapply(seq_len(d1), function(k){
  runif(n = n, min = mini[k], max = maxi[k])
}) %>% as.data.frame()
colnames(X1) = paste0("X", seq_len(d1))
X = X1

#calcul des Y selon une relation lineaire + du bruit
p = 2
fn = lapply(seq_len(p), function(j){
  beta = runif(n = d1, min = -2, max = 2)
  function(X) {
    as.matrix(X) %*% beta
  }
})
names(fn) = paste0("Y", seq_len(p))

Y = lapply(fn, function(f){
  f(X) #+ rnorm(n)
}) %>% as.data.frame(row.names = seq_len(n))

cost_fn = list(c1 = function(X, Y, y_need = 1) Y[,1],
               c2 = function(X, Y, y_need = 1:2) 3 * Y[,2] - Y[,1])

g_cstr = list(
  function(x){
    x = as.data.frame(x)
    x[,1] + x[,2] + x[,3] < 2
  }
)
res = optisure(X, Y, fn = cost_fn, alpha = 0.05, g = g_cstr, X_space_csrt=T, B = 10)
plot(res)



#avec space restriction
mini = rep(0, d)
maxi = rep(10, d)

#runif
b = matrix(c(-2, -2, -2, 0, 0, 0), byrow = T, ncol = d)
B = matrix(c(0, 0, 0, 2, 2, 2), byrow = T, ncol = d)
nb_foyer = NROW(B) #nombre de "foyer"
X1 = sapply(seq_len(d), function(i){
  sapply(seq_len(nb_foyer), function(k){
    runif(n/nb_foyer, b[k,i], B[k,i])
  })
})
X = X1
plot(as.data.frame(X))

Y = lapply(fn, function(f){
  f(X) #+ rnorm(n)
}) %>% as.data.frame(row.names = seq_len(n))

res = optisure(X, Y, fn = cost_fn, alpha = 0.5, g = g_cstr, X_space_csrt=T,
               B = 20, parametric = TRUE)
plot(res)

#X quali----
d2 = 2
lev = sample(x = 2:4, d2, replace = TRUE)
X2 <- sapply(seq_len(d2), function(k){
  as.factor(sample(x = seq_len(lev[k]), size = n, replace = TRUE))
}) %>% as.data.frame()
colnames(X2) = paste0("X", (d1+1):(d1+d2))

X = cbind(X1, X2)


#calcul des Y selon une relation lineaire + du bruit
fn = lapply(seq_len(p), function(j){
  beta = runif(n = (d1 + sum(lev-1)), min = -2, max = 2)
  function(X) {
    is_num = sapply(X, is.numeric)
    X1 = X[, is_num]
    X2 = X[, !is_num]
    X2 = model.matrix(~., data = X2)[,-1]
    X = cbind(X1, X2)
    as.matrix(X) %*% beta
  }
})
names(fn) = paste0("Y", seq_len(p))

cost_fn = list(c1 = function(X, Y, y_need = 1) Y[,1],
               c2 = function(X, Y, y_need = 1:2) 3 * Y[,2] - Y[,1])

res = optisure(X, Y, fn = cost_fn, alpha = 0.05, g = g_cstr, X_space_csrt=T,
               B = 10, parametric=TRUE)


plot(res)

