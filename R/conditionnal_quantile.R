#' @export
conditionnal_quantile <- function(X, Y, x, alpha,
                                  h_n = 10*(max(X)-min(X))/length(X)){

  F_n <- function(x, y){
    x_mat <- matrix(rep(unlist(x), NROW(X)), ncol = NCOL(x), byrow = TRUE)
    d <- (x_mat - X)/h_n
    Nn <- sum(sapply(t(d), K) * as.numeric(Y <= y))
    Dn <- sum(sapply(t(d), K))
    Nn/Dn
  }

  P <- function(y){
    (alpha - F_n(x = x, y = y))^2
  }

  optimize(P, c(min(Y), max(Y)), maximum = FALSE)$minimum
}
