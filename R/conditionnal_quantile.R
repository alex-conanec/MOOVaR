#' conditionnal_quantile
#'
#' Non parametric estimation function of the conditional alpha-quantile
#'
#' @param X a numerical vector, matrix, or data.frame of explanatory variables  of dimension nxq
#' @param Y a numerical vector of response variable of size n
#' @param x the conditional value of the random variable X
#' @param alpha the alpha value of the quantile of regression
#' @param h_n the mobile windows size use to smooth the distribution function.
#' The parameter is tuned if NULL (default value).
#'
#' @param iter_max the maximum iteration allowed to find the root and inverse the CDF.
#' 15 by default.
#' @param tol the tolerance threshold below which the convergence is reached. 0.001 by default.
#'
#' @return the value of the conditional quantile
#'
#' @examples
#' library(tidyverse)
#'
#' set.seed(258164)
#' n <- 300
#' X <- runif(n = n, min = -2, max = 2)
#' Y <- X^2 + rnorm(n)
#'
#' #parametres
#' alpha <- 0.1
#' h_n <- hopt(X, Y)$h
#'
#' #estimation de q_n et calcul de q
#' x_range <- seq(-2, 2, length.out = n/3)
#'
#' res = conditionnal_quantile(X, Y, x=x_range, alpha, h_n = h_n,
#'                               iter_max = 15, tol = 0.01)
#'
#' data = data.frame(
#'   x = x_range,
#'   q = qnorm(p = alpha, mean = x_range^2, sd = 1),
#'   q_n = res$y,
#'   X,
#'   Y
#' ) %>% gather(key = key, value = value, -x, -X, -Y)
#'
#' #plot
#' ggplot(data) +
#'   geom_point(aes(x=X, y=Y)) +
#'   geom_line(aes(x=x, y=value, col = key))
#'
#' #ecart absolu a 0 dans la recherche de la racine
#' boxplot(abs(res$`f(y)`))
#'
#' #nombre d'iteration necessaire pour converger
#' ggplot(res) + geom_bar(aes(iter))
#'
#'
#' #q=3
#' set.seed(123)
#' q = 3
#' mini <- rep(x = -2, times = q)
#' maxi <- rep(x = 2, times = q)
#'
#' X <- lapply(seq_len(q), function(k){
#'   res <- data.frame(X=runif(n = n, min = mini[k], max = maxi[k]))
#'   names(res) <- paste0("X", k)
#'   res
#' }) %>% bind_cols() %>% as.matrix()
#'
#' beta = runif(3)
#' Y = X %*% beta + rnorm(NROW(X))
#'
#' h_n = hopt(X, Y)$h
#' res = conditionnal_quantile_2(X, Y, x=X[1:20,], alpha, h_n = h_n,
#'                               iter_max = 15, tol = 0.01)
#'
#' x_q = 1 #jouer avec le q ici pour changer de coupe dans la dimension
#' plot(X[,x_q], Y)
#'
#' points(data.frame(x = X[1:20, x_q], y = qnorm(alpha, mean = X[1:20,] %*% beta)) %>%
#'          arrange(x), type = "b", col = "red")
#' points(data.frame(x = X[1:20, x_q], y = res$y) %>%
#'          arrange(x), type = "b", col = "blue")
#'
#'
#' boxplot(res$y - qnorm(alpha, mean = X[1:20,] %*% beta))
#'
#' @export

conditionnal_quantile <- function(X, Y, x, alpha, h_n = NULL, iter_max = 15,
                                    tol = 0.001){

  ##bisection method

  # find the root of f
  f <- function(y, x, alpha){
    alpha - F_n(X = X, Y = Y, x = x, y = y, h_n = h_n)
  }

  x = as.matrix(x)
  q = NROW(x)

  if (q > 1 & length(alpha) == 1) alpha = rep(alpha, q)

  a = rep(max(Y), q)
  b = rep(min(Y), q)

  names(a) = seq_len(q)
  res = data.frame(matrix(NA, ncol = 3, nrow = q))
  colnames(res) = c("iter", "y", "f(y)")

  for (i in seq_len(iter_max)){

    c = (a+b)/2

    f_c = f(y = c, x, alpha)

    # tolerance
    done = abs(f_c) < tol
    row_done = as.numeric(names(c)[done])

    res[row_done, 1] = i
    res[row_done, 2] = c[done]
    res[row_done, 3] = f_c[done]

    if (all(done)) break

    a = a[!done]; b = b[!done]; c = c[!done]; x = x[!done,, drop = F]; f_c = f_c[!done]; alpha = alpha[!done]

    neg = sign(f_c) < 0
    pos = sign(f_c) > 0
    a[neg] = c[neg]
    b[pos] = c[pos]

  }

  not_done_yet = is.na(res$iter)
  res[not_done_yet, 1] = i
  res[not_done_yet, 2] = c
  res[not_done_yet, 3] = f_c

  res
}

