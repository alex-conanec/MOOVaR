library(dplyr)
test_that("give the right quantile conditionnal estimation for dim(Y)=1", {

  #test q = 1
  n <- 300
  X <- runif(n = n, min = -2, max = 2)
  Y <- X^2 + rnorm(n)

  #parametres
  h_n <- hopt(X, Y)$h

  N = 100
  x_range <- runif(n = N, min = -2, max = 2)
  alpha_range = runif(n = N, min = 0.01, max = 0.99)
  q_n = conditionnal_quantile(X, Y, x = x_range,
                              alpha=alpha_range,
                              # alpha=0.05,
                              h_n = h_n,
                              iter_max = 15, tol = 0.01)$y
  q_true = qnorm(alpha_range, mean = x_range^2, sd = 1)

  expect_that(t.test(q_n, q_true, paired = T)$p.value, is_more_than(0.05))


  #test q = 2
  set.seed(12) #123 NS diff #12 p.value=0.001
  n <- 300
  d = 3
  X <- lapply(seq_len(d), function(j){
    data.frame(runif(n = n, min = -2, max = 2))
  }) %>% bind_cols()
  colnames(X) = paste0("X", 1:d)

  Y <- rowSums(X^2) + rnorm(n)

  x_range <- matrix(runif(n = d*N, min = -2, max = 2), ncol = d)
  alpha_range = runif(n = N, min = 0.01, max = 0.99)
  h_n <- hopt(X, Y)$h

  q_n = conditionnal_quantile(X, Y, x = x_range,
                              alpha = alpha_range,
                              h_n = h_n,
                              iter_max = 15, tol = 0.01)$y
  q_true = qnorm(alpha_range, mean = rowSums(x_range^2), sd = 1)


  expect_that(t.test(q_n, q_true, paired = T)$p.value, is_more_than(0.05))



  data.frame(
    q_n_Q1 = conditionnal_quantile(X, Y, x = x_range,
                                   alpha = 0.25,
                                   h_n = h_n,
                                   iter_max = 15, tol = 0.01)$y,
    q_n_Q3 = conditionnal_quantile(X, Y, x = x_range,
                                   alpha = 0.75,
                                   h_n = h_n,
                                   iter_max = 15, tol = 0.01)$y,
    y =  rowSums(x_range^2)
  ) %>% arrange(y) %>%
    mutate(id = 1:NROW(.)) %>%
    gather(key = type, value=q_value, -id, -y) %>%
    ggplot(aes(x=id, y=value)) +
    geom_point(aes(x = id, y = y)) +
    geom_line(aes(x = id, y = q_value, colour = type))


  data.frame(
    id = 1:length(q_true),
    q_n,
    q_true
  ) %>%
    gather(key = type, value=value, -id) %>%
    mutate(id = as.factor(id)) %>%
    ggplot(aes(x=id, y=value)) +
    geom_point(aes(colour = type))+
    geom_line(arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"))
  #
  #
  # data.frame(
  #   id = 1:length(q_true),
  #   q_n,
  #   q_true
  # ) %>%
  #   mutate(diff = q_true - q_n) %>%
  #   ggplot(aes(x = q_true, y = diff)) +
  #   geom_point()
  #

})
#
# set.seed(1234) #12
# d1 = 3
# mini <- rep(x = -2, times = d1)
# maxi <- rep(x = 2, times = d1)
# X <- sapply(seq_len(d1), function(k){
#   runif(n = n, min = mini[k], max = maxi[k])
# }) %>% as.data.frame()
# colnames(X) = paste0("X", seq_len(d1))
#
# #calcul des Y selon une relation lineaire + du bruit
# p = 2
# fn = lapply(seq_len(p), function(j){
#   beta = runif(n = d1, min = -2, max = 2)
#   print(beta)
#   function(X) {
#     as.matrix(X) %*% beta
#   }
# })
# names(fn) = paste0("Y", seq_len(p))
#
# Y = lapply(fn, function(f){
#   f(X) + rnorm(n)
# }) %>% as.data.frame(row.names = seq_len(n))
#
# train = sample(seq_len(n), size = 2*n/3)
#
# h_n = lapply(seq_len(p), function(j) hopt(X = X[train,], Y = Y[train, j])$h )
# qj_x = lapply(seq_len(p), function(j){
#   function(x){
#     conditionnal_quantile(X = X,
#                           Y = Y[,j],
#                           x = x,
#                           alpha,
#                           h_n = h_n[[j]])$y
#                           # h_n = h_n_bis[[j]])$y
#   }
# })
# Y_all=Y
# Y = Y_all[,2]
# h_n_bis = list(rep(0.2, 3), rep(0.2, 3))
# # x = data.frame(-1, -1, -1)
# plot(Y)
# points(qj_x[[1]](X), qj_x[[2]](X), col = 'green')
# h_n
#
#
# res = data.frame(
#   q_n_Q1 = conditionnal_quantile(X[train,], Y[train,2], x = X[-train,],
#                                  alpha = 0.5,
#                                  h_n = h_ks,
#                                  iter_max = 15, tol = 0.01)$y,
#   # q_n_Q3 = conditionnal_quantile(X[train,], Y[train,1], x = X[-train,],
#   #                                alpha = 0.75,
#   #                                h_n = h_n[[1]],
#   #                                iter_max = 15, tol = 0.01)$y,
#   y =  Y[-train,2]
# )
#
# sum((res$q_n_Q1 - res$y)^2)
#
# library(ks)
# h_ks = diag(ks::Hpi(x = X))
#
# res %>% arrange(y) %>%
#   mutate(id = 1:NROW(.)) %>%
#   gather(key = type, value=q_value, -id, -y) %>%
#   ggplot(aes(x=id, y=value)) +
#   geom_point(aes(x = id, y = y)) +
#   geom_line(aes(x = id, y = q_value, colour = type))
#
#
#
# points(Y, col = 'green')
# alpha=0.6
#
