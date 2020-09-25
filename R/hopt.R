#' mobile windows h tuning
#'
#' Tuning of the mobile windows h to estimate the cumulative distribution function (cdf)
#'
#' @param X a numerical vector, matrix, or data.frame of explanatory variables of dimension nxq.
#' @param Y a numerical vector of response variable of size n.
#' @param h0 initialisation of the h value of the optimisation process. By default equals to rep(0.5, p).
#' @param lower a numerical vector low bownding hopt. By default equals to rep(0.01, q).
#' @param upper a numerical vector low bownding hopt. By default equals to rep(2, q).
#'
#' @return a list of the parameter optimal value h and the criteria minimum reach
#'
#' @examples
#' library(tidyverse)
#' library(condQuant)
#'
#' set.seed(258164)
#' n <- 300
#' X <- runif(n = n, min = -2, max = 2)
#' Y <- X^2 + rnorm(n)
#' X <- as.data.frame(X)
#'
#' #critere optimiser manuellement
#' critere_hopt = function(h, Xo, Yo){
#'
#'   q = NCOL(Xo)
#'   n = NROW(Xo)
#'
#'   #calcul de toutes les distances une à une ponderer par la largeur de la fenetre
#'   b = lapply(seq_len(q), function(j){
#'     A = matrix(Xo[,j], nrow = n, ncol = n, byrow = F)
#'     (A - t(A)) / h[j]
#'   })
#'
#'   # transformation de la liste sous forme de array pour utiliser apply apres
#'   d = array(unlist(b), dim = c(n, n, q))
#'
#'   # application du noyaux aux distances sotckees dans d
#'   w = (2*pi)^(-q/2)*exp(-apply(d^2, c(1,2), sum)/2)
#'
#'   # indicatrice pour chaque y en pairwise
#'   A = matrix(Yo, nrow = length(Yo), ncol = length(Yo), byrow = F)
#'   ind_y = apply(A - t(A) < 0, c(1,2), as.numeric)
#'
#'   #double (triple) somme vectorisee
#'   diag(w) = 0
#'   B = w %*% t(ind_y)
#'   diag(B) = 0
#'   sum( (t(ind_y) - B/colSums(w))^2 )
#'
#' }
#'
#' formals(critere_hopt)$Xo = X
#' formals(critere_hopt)$Yo = Y
#'
#' h_range <- seq(0.01, 0.5, 0.05)
#'
#' res <- sapply(h_range, critere_hopt)
#' plot(x = h_range, y = res, type = 'b')
#'
#' hopt(X = X, Y = Y)
#'
#' @export

hopt <- function(X, Y, h0 = NULL, lower = NULL, upper = NULL){

  X <- as.data.frame(X)
  q <- NCOL(X)

  width = apply(X, 2, max) - apply(X, 2, min)
  if (is.null(h0)) h0 = width/10
  if (is.null(upper)) upper = width
  if (is.null(lower)) lower = width/200

  critere_hopt = function(h, Xo, Yo){

    q = NCOL(Xo)
    n = NROW(Xo)

    #calcul de toutes les distances une à une ponderer par la largeur de la fenetre
    b = lapply(seq_len(q), function(j){
      A = matrix(Xo[,j], nrow = n, ncol = n, byrow = F)
      (A - t(A)) / h[j]
    })

    # transformation de la liste sous forme de array pour utiliser apply apres
    d = array(unlist(b), dim = c(n, n, q))

    # application du noyaux aux distances sotckees dans d
    w = (2*pi)^(-q/2)*exp(-apply(d^2, c(1,2), sum)/2)

    # indicatrice pour chaque y en pairwise
    A = matrix(Yo, nrow = length(Yo), ncol = length(Yo), byrow = F)
    ind_y = apply(A - t(A) < 0, c(1,2), as.numeric)

    #double (triple) somme vectorisee
    diag(w) = 0
    B = w %*% t(ind_y)
    diag(B) = 0
    sum( (t(ind_y) - B/colSums(w))^2 )

  }

  formals(critere_hopt)$Xo = X
  formals(critere_hopt)$Yo = Y

  if (q > 1){
    res <- optim(par = h0, fn = critere_hopt)
    # res <- optim(par = h0, fn = critere_hopt, lower = lower, upper = upper,
    #              method = "L-BFGS-B" )
    list(h = res$par, min_criteria = res$value)

  }else{
    res <- optimize(f = critere_hopt, interval = c(lower, upper), maximum = FALSE)
    list(h = res$minimum, min_criteria = res$objective)
  }
}

