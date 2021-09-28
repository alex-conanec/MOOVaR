#' quant_reg_global_risk
#'
#' @param X data.frame of the dependent variables.
#' @param Y data.frame of the response variables.
#' @param tau vector of size NCOL(Y) containing the risk in ]0,1[ of each Y.
#' @param floor threshold to prevent estimate hard quantile with not enought data
#' @param max_iter integer to stop search if not converged. 20 by default.
#' @param path_tracking path where to write the step of the running function.

#' @export
quant_reg_global_risk = function(X, Y, tau, allowed_dependence,
                                 initial_eta = 2.5, k = 0.2,
                                 floor = 10/NROW(X), max_iter = 20, tune = FALSE,
                                 path_tracking = NULL){

  p = NCOL(Y)
  lambda = rep(0, p)
  lambda_prec = rep(1, p)

  if (tune){
    a = optim(par = c(2.5, 0.2), fn = tune_step)
    a$value
    initial_eta = a$par[1]
    k = a$par[2]
  }

  eta = function(t){
    initial_eta * exp(-k*t)
  }

  res = list(globale_confiance =  ifelse(1 - tau > 0.5, 0, 1))

  mask = lapply(1:p, function(j){
    which(!duplicated(X[,c(TRUE, allowed_dependence[,j]), drop = FALSE]) &
      !is.na(Y[,j,T]))
  }) %>% Reduce(f = intersect, x = .)

  n = length(mask)

  X_test = X[,-1,F]
  for (i in which(!apply(is.na(allowed_dependence), 1, any))){
    a = X[mask,-1,F][,i,T]
    X_test[is.na(X_test[,i,T]),i] = a[!is.na(a)][1]
  }

  for (t in seq_len(max_iter)){
    m = fit_model(X=X[mask,], Y=Y[mask,],
                  tau = tau + lambda,
                  method = "quantile",
                  allowed_dependence = allowed_dependence,
                  path_tracking = path_tracking)

    #eval
    delta = as.matrix(Y[mask,]) - as.matrix(m(X_test[mask,]))
    cat(colSums(delta > 0)/n, "\n")
    L = sum(apply(delta > 0, 1, all))/n

    if (abs(L - 1 + tau) < abs(res$globale_confiance - 1 + tau) ){
      cat("ecart:", abs(L - 1 + tau), "\n")
      res = list(globale_confiance = L, tau = tau + lambda)
    }

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
    lambda[lambda + tau > 1 - floor] = 1 - floor - tau

    if (all(abs(lambda - lambda_prec) < tau * 0.001) |
        abs(res$globale_confiance - 1 + tau) < 0.005){
      break
    }
    lambda_prec = lambda

  }

  res$m = fit_model(X, Y,
                    tau = res$tau,
                    method = "quantile",
                    allowed_dependence = allowed_dependence,
                    path_tracking = path_tracking)
  res
}

tune_step = function(x){
  L = quant_reg_global_risk(X, Y, tau, allowed_dependence,
                            initial_eta = x[1], k = x[2],
                            floor = 5/NROW(X), max_iter = 5, tune = FALSE,
                            path_tracking = NULL)$globale_confiance
  abs(L - 1 + tau)
}
