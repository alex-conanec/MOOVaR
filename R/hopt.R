#' mobile windows h tuning
#'
#' Tuning of the mobile windows h to estimate the cumulative distribution function (cdf)
#'
#' @param X a numerical vector, matrix, or data.frame of explanatory variables of dimension nxq.
#' @param Y a numerical vector of response variable of size n.
#' @param n_tilde a positive integer lower than n, to perform a cross-validation
#' with a less expessive cost than with leave-one-out method.
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
#' critere_hopt <- function(h, E){
#'   X_tilde <- X[E, , drop = FALSE]
#'   Y_tilde <- Y[E]
#'   n_range <- seq_len(n_tilde)
#'   parallel::mclapply(n_range, function(k){
#'     sapply(n_range[-k], function(i){
#'       ( as.numeric(Y_tilde[k] <= Y_tilde[i]) -
#'           F_n(X = X, Y = Y,
#'               x = X_tilde[k,], y = Y_tilde[i],
#'               h_n = h, k = E[k]) )^2
#'     }) %>% mean()
#'   }, mc.cores = 8) %>% unlist() %>% sum()
#' }
#'
#' n_tilde <- 50 #play with it to find a trade-off between time computing and h stability
#' E <- sample(x = seq_len(n), n_tilde)
#' h_range <- seq(0.01, 0.8, 0.1)
#'
#' res <- sapply(h_range, critere_hopt, E = E)
#' plot(x = h_range, y = res, type = 'b')
#'
#' hopt(X = X, Y = Y, n_tilde = 100)

#' @export

hopt <- function(X, Y, n_tilde = NULL, h0 = rep(0.5, NCOL(X)),
                 lower = rep(0.01, NCOL(X)), upper = rep(2, NCOL(X))){

  X <- as.data.frame(X)
  q <- NCOL(X)
  n <- NROW(X)

  if (is.null(n_tilde)) n_tilde <- n

  critere_hopt <- function(h, E){
    X_tilde <- X[E, , drop = FALSE]
    Y_tilde <- Y[E]
    n_range <- seq_len(n_tilde)
    parallel::mclapply(n_range, function(k){
      sapply(n_range[-k], function(i){
        ( as.numeric(Y_tilde[k] <= Y_tilde[i]) -
            F_n(X = X, Y = Y,
                x = X_tilde[k,], y = Y_tilde[i],
                h_n = h, k = E[k]) )^2
      }) %>% mean()
    }, mc.cores = 8) %>% unlist() %>% sum()
  }

  E <- sample(x = seq_len(n), n_tilde)
  formals(critere_hopt)$E <- E


  if (q > 1){
    res <- optim(par = h0,
          # lower = lower, upper = upper,
          fn = critere_hopt)
    list(h = res$par, min_criteria = res$value)

  }else{
    res <- optimize(f = critere_hopt, interval = c(lower, upper), maximum = FALSE)
    list(h = res$minimum, min_criteria = res$objective)
  }
}
