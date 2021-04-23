#' @export
quant_reg_global_risk = function(X, Y, tau, floor = 0.1, max_iter = 20, ...){

  p = NCOL(Y)
  n = NROW(Y)
  lambda = rep(0, p)
  lambda_prec = rep(1, p)

  eta = function(t, initial_eta = 1.5, k = 0.05){
    initial_eta * exp(-k*t)
  }

  for (t in seq_len(max_iter)){
    m = fit_model(X, Y, tau + lambda, method = "quantile", ...)

    #eval
    delta = Y - m(X)
    cat(colSums(delta > 0)/n, "\n")
    L = sum(apply(delta > 0, 1, all))/n
    cat(L, "\n")

    #evolution lambda
    if (p > 2){
      S = sapply(seq_len(p), function(j){
        sum(apply(delta[, -j, drop = FALSE] > 0, 1, all))/n
      })
      S = S/sum(S)
      lambda = lambda + (L - (1 - tau)) * S * eta(t)
    }else{
      lambda = lambda + (L - (1 - tau)) * eta(t)
    }


    #avoid learn too high repartition = avoid estimation pb
    lambda[lambda + tau < floor] = floor - tau

    if (all(abs(lambda - lambda_prec) < tau * 0.001)){
      break
    }
    lambda_prec = lambda

  }


  list(model = m, globale_confiance = L, tau = tau + lambda)
}
