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
  X <- data.frame(X1 = runif(n = n, min = -2, max = 2),
                  X2 = runif(n = n, min = -2, max = 2))

  Y <- rowSums(X^2) + rnorm(n)

  x_range <- matrix(runif(n = 2*N, min = -2, max = 2), ncol = 2)
  alpha_range = runif(n = N, min = 0.01, max = 0.99)
  h_n <- hopt(X, Y)$h

  q_n = conditionnal_quantile(X, Y, x = x_range,
                              alpha = alpha_range,
                              h_n = h_n,
                              iter_max = 15, tol = 0.01)$y
  q_true = qnorm(alpha_range, mean = rowSums(x_range^2), sd = 1)


  expect_that(t.test(q_n, q_true, paired = T)$p.value, is_more_than(0.05))


  # data.frame(
  #   id = 1:length(q_true),
  #   q_n,
  #   q_true
  # ) %>%
  #   gather(key = type, value=value, -id) %>%
  #   mutate(id = as.factor(id)) %>%
  #   ggplot(aes(x=id, y=value)) +
  #   geom_point(aes(colour = type))+
  #   geom_line(arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"))
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

