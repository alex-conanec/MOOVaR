library(tidyverse)
library(optisure)

test_that("F_n cumulative distribution fonction give similar results closed enough to the reality", {

  #test q = 1
  n <- 300
  X <- runif(n = n, min = -2, max = 2)
  Y <- X^2 + rnorm(n)

  #parametres
  h_n <- hopt(X, Y)$h

  N = 100
  x_range <- runif(n = N, min = -2, max = 2)
  y_range <- runif(n = N, min = -2, max = 6)

  p_n = F_n(X = X, Y = Y, x = x_range, y = y_range, h_n = h_n)
  p_true = pnorm(y_range, mean = x_range^2, sd = 1)


  expect_that(t.test(p_n, p_true, paired = T)$p.value, is_more_than(0.05))

  #test q = 2
  set.seed(12)
  n <- 300
  X <- data.frame(X1 = runif(n = n, min = -2, max = 2),
                  X2 = runif(n = n, min = -2, max = 2))

  Y <- rowSums(X^2) + rnorm(n)

  x_range <- matrix(runif(n = 2*N, min = -2, max = 2), ncol = 2)
  y_range <- runif(n = N, min = -2, max = 6)
  h_n <- hopt(X, Y)$h

  p_n = F_n(X = X, Y = Y, x = x_range, y = y_range, h_n = h_n)
  p_true = pnorm(y_range, mean = rowSums(x_range^2), sd = 1)

  expect_that(t.test(p_n, p_true, paired = T)$p.value, is_more_than(0.05))

  # pb = data.frame(
  #   id = 1:length(p_true),
  #   x_range,
  #   y_range,
  #   p_n,
  #   p_true
  # ) %>%
  #   mutate(diff_abs = abs(p_true - p_n),
  #          diff = p_true - p_n)
  #
  #
  # ggplot(pb, aes(x = p_true, y = diff)) +
  #   geom_point()
  #
  # ggplot(pb, aes(x = X1, y = X2)) +
  #   geom_point(aes(size = diff_abs, colour = diff_abs)) +
  #   geom_point(data = data.frame(X, Y), aes(x = X1, y = X2),
  #              colour = "red", size = .5)
  #
  # ggplot(pb, aes(x = X1, y = y_range)) +
  #   geom_point(aes(size = diff_abs, colour = diff_abs)) +
  #   geom_point(data = data.frame(X, Y), aes(x = X1, y = Y),
  #              colour = "red", size = .5)
  #
  # ggplot(pb, aes(x = X2, y = y_range)) +
  #   geom_point(aes(size = diff_abs, colour = diff_abs)) +
  #   geom_point(data = data.frame(X, Y), aes(x = X2, y = Y),
  #              colour = "red", size = .5)

})

