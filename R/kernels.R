#' @export
omega <- function(u) as.numeric(abs(u) <= 1)/2

#' @export
K <- function(x){
  k <- length(x)
  # sqrt(2*pi)^(-k) * exp(1/k * -sum(x^2) / 2)
  sqrt(2*pi)^(-k) * exp(-sum(x^2) / 2)
}
