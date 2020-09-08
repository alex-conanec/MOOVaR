#' Cumulative distribution function (cdf)
#'
#' Non parametric estimation of Cumulative distribution function (cdf).
#' This can be express as P(Y<=y|X=x)
#'
#' @param X a numerical vector, matrix, or data.frame of explanatory variables of dimension nxq
#' @param Y a numerical vector of response variable of size n
#' @param x the conditional value of the random variable X
#' @param y the value for which we want to evaluate the probability of being
#' lower or equal to the random variable Y knowing that X=x
#' @param h_n the mobile windows size use to smooth the distribution function.
#' @param k a positive integer lower than n to exclude observation from the training set in
#' order to perform cross-validation. NULL by default.
#'
#'
#' @return P(Y<=y|X=x)
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
#' h_n <- 0.2
#'
#' #estimation
#' F_n(X = X, Y = Y, x = x, y = y, h_n = h_n)
#' pnorm(y, mean = x^2, sd = 1) #true probability
#'
#' @export
F_n <- function(X, Y, y, x, h_n = NULL, k = NULL){
  q <- NCOL(x)
  X <- as.data.frame(X)

  if (!is.null(k)){
    X <- X[-k, , drop = FALSE]
    Y <- Y[-k]
    n <- n-1
  }

  if (is.null(h_n)){
    h_n <- hopt(X = X, Y = Y, n_tilde = 100)$h
  }

  d = sweep(
    -sweep(X, MARGIN = 2, FUN = "-", STATS = unlist(x)),
    MARGIN = 2, FUN = "/", STATS = unlist(h_n))

  w = pi^(-q/2)*exp(-rowSums(d^2)/2)
  Nn <- sum(w * as.numeric(Y <= y))
  Dn <- sum(w)
  Nn/Dn
}
