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
#' h_n <- hopt(X, Y)$h
#'
#' #estimation
#' F_n(X = X, Y = Y, x = x, y = y, h_n = h_n)
#' pnorm(y, mean = x^2, sd = 1) #true probability
#'
#' @export

F_n <- function(X, Y, y, x, h_n = NULL){

  X <- as.data.frame(X)
  q <- NCOL(X)
  n <- NROW(X)
  m = NROW(x)

  if (is.null(h_n)){
    h_n <- hopt(X = X, Y = Y)$h
  }

  #repeter chaque ligne de x n*m fois
  s=array(rep(t(as.matrix(x)), each = n), c(n, q, m))

  # calcul des distances
  apply(s, MARGIN = 3, function(u){
    (u - X)/h_n
  }) %>% unlist() %>% array(dim = c(n, q, m)) -> d

  # application du noyau gaussien
  # w = pi^(-q/2)*exp(-apply(X = d, MARGIN = c(2,3), FUN = function(x) sum(x^2)/2))
  w = pi^(-q/2)*exp(-apply(X = d, MARGIN = c(1,3), FUN = function(x) sum(x^2)/2))

  # calcule de l'indicatrice
  ind_y = matrix( apply(matrix(Y, ncol = m, nrow = n, byrow = FALSE), MARGIN = 1,
        function(x) as.numeric(x <= y)), nrow = m)

  # somme du produit scalaire entre w et l'indicatrice puis quotien
  Nn <- colSums(w * t(ind_y))
  Dn <- colSums(w)
  Nn/Dn
  # return(1)
}
