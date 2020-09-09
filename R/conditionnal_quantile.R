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
#' h_n <- 0.2
#'
#' #estimation de q_n et calcul de q
#' x_range <- seq(-2, 2, length.out = n/2)
#'
#' parallel::mclapply(x_range, function(x){
#'   q <- qnorm(p = alpha, mean = x^2, sd = 1)
#'   q_n <- conditionnal_quantile(X, Y, x, alpha, h_n = h_n)
#'   data.frame(q=q, q_n=q_n)
#' }, mc.cores = 8) %>% bind_rows()-> res
#'
#' #plot
#' res <- cbind(x_range = x_range, res) %>%
#'   gather(key = key, value = value, -x_range)
#' data <- cbind(data.frame(X=X, Y=Y), res)
#'
#' ggplot(data) +
#'   geom_point(aes(x=X, y=Y)) +
#'   geom_line(aes(x=x_range, y=value, col = key))
#'
#' @export

conditionnal_quantile <- function(X, Y, x, alpha, h_n = NULL){

  P <- function(y){
    (alpha - F_n(X = X, Y = Y, x = x, y = y, h_n = h_n))^2
  }

  optimize(P, c(min(Y), max(Y)), maximum = FALSE)$minimum
}


