#' optisure
#'
#' @examples
#' # generate data ----
#' library(tidyverse)
#' library(optisure)
#' set.seed(1234)
#' n <- 300
#'
#' # X quanti ----
#' # X
#' d1 = 3
#' mini <- rep(x = -2, times = d1)
#' maxi <- rep(x = 2, times = d1)
#' X1 <- sapply(seq_len(d1), function(k){
#'   runif(n = n, min = mini[k], max = maxi[k])
#' }) %>% as.data.frame()
#' colnames(X1) = paste0("X", seq_len(d1))
#' X = X1
#'
#' #calcul des Y selon une relation lineaire + du bruit
#' p = 2
#' fn = lapply(seq_len(p), function(j){
#'   beta = runif(n = d1, min = -2, max = 2)
#'   function(X) {
#'     as.matrix(X) %*% beta
#'   }
#' })
#' names(fn) = paste0("Y", seq_len(p))
#'
#' Y = lapply(fn, function(f){
#'   f(X) + rnorm(n)
#' }) %>% as.data.frame(row.names = seq_len(n))
#'
#' # Opti
#' res = optisure(X, Y, N=20, B=10, alpha = 0.5)
#'
#' plot(ellipse_var(Y$Y1, Y$Y2, alpha = 0.05), type = 'l', xlab = "Y1", ylab = "Y2",
#'      xlim = c(-8,8), ylim = c(-10,6))
#' points(Y, col = "grey")
#' points(res$Y, type = 'b', col = "blue", lwd = 3)
#'
#' #X quali----
#' d2 = 2
#' lev = sample(x = 2:4, d2, replace = TRUE)
#' X2 <- sapply(seq_len(d2), function(k){
#'   as.factor(sample(x = seq_len(lev[k]), size = n, replace = TRUE))
#' }) %>% as.data.frame()
#' colnames(X2) = paste0("X", (d1+1):(d1+d2))
#'
#' X = cbind(X1, X2)
#' head(X)
#'
#' #calcul des Y selon une relation lineaire + du bruit
#' fn = lapply(seq_len(p), function(j){
#'   beta = runif(n = (d1 + sum(lev-1)), min = -2, max = 2)
#'   function(X) {
#'     is_num = sapply(X, is.numeric)
#'     X1 = X[, is_num]
#'     X2 = X[, !is_num]
#'     X2 = model.matrix(~., data = X2)[,-1]
#'     X = cbind(X1, X2)
#'     as.matrix(X) %*% beta
#'   }
#' })
#' names(fn) = paste0("Y", seq_len(p))
#'
#' #opti
#' res = optisure(X, Y, N=20, B=10, alpha = 0.5)
#'
#' plot(ellipse_var(Y$Y1, Y$Y2, alpha = 0.05), type = 'l', xlab = "Y1", ylab = "Y2",
#'      xlim = c(-8,8), ylim = c(-10,6))
#' points(Y, col = "grey")
#' points(res$Y, type = 'b', col = "blue", lwd = 3)
#'
#' @importFrom magrittr %>%
#' @export
optisure <- function(X, Y, sens = rep("min", NCOL(Y)),
                     utility_risk = rep("quantile", NCOL(Y)),
                     g = NULL, X_space_csrt = TRUE,
                     tau = 0.5, globale_tau = FALSE, reg_method = "linear", #rajouté
                     # parametric = FALSE, deg = 3, #a virer
                     N = NULL, B = NULL, seed = NULL,
                     ...){

  library(reticulate)
  d = NCOL(X)
  p = NCOL(Y)
  n = NROW(X)
  col_names_Y = colnames(Y)

  #sens optimisation
  Y = as.data.frame(Y)
  # Y = Y * sapply(sens, function(x) ifelse(x == "min", -1, 1))
  Y = sweep(Y, 2, sapply(sens, function(x) ifelse(x == "min", -1, 1)), "*")

  quantile_utility_idx = utility_risk == rep("quantile", NCOL(Y))

  if (any(quantile_utility_idx)){
    if (globale_tau){
      qrgr = quant_reg_global_risk(X = X, Y[,quantile_utility_idx, drop = FALSE],
                                   tau = tau, reg_method = reg_method) #pas sur de la syntax pour transmettre reg_method
      m_quantile = qrgr$m
    }else{
      m_quantile = fit_model(X, Y=Y[,quantile_utility_idx, drop = FALSE],
                             tau = tau, reg_method = reg_method,
                             method = "quantile", ...)
      # m_quantile = fit_model(X, Y=Y[,quantile_utility_idx, drop = FALSE],
      #                        tau = tau, reg_method = reg_method,
      #                        method = "quantile", penalty = rep(0,p))
    }
  }

  if (any(!quantile_utility_idx)){
    m_expect = fit_model(X, Y = Y[,!quantile_utility_idx, drop = FALSE],
                         reg_method = reg_method, method = "expected")
  }

  if (all(quantile_utility_idx)){
    m = m_quantile
  }else if (all(!quantile_utility_idx)) {
    m = m_expect
  }else{
    m = function(X){
      res = cbind(m_quantile(X), m_expect(X))[, c(which(quantile_utility_idx),
                                                  which(!quantile_utility_idx))]
      colnames(res) = colnames(Y)
      res
    }
  }


  #tunage nsga parametres
  if (is.null(N)) N = 50 ############
  if (is.null(B)) B = 20

  if (X_space_csrt){
    cat("Train belonging to density constraint \n")
    g$cstr_X_space = def_cstr_X_space(X)$g
  }


  #appelle de NSGA ou optim si NCOL(Y)==1
  res = NSGA(X = X, fn = m, n_objective = p, g=g,
             sens = rep("max", p),
             N = N, seed = seed,
             B = B, verbose = TRUE)

  #mise en forme des resultats
  if (globale_tau){
    res$tau_j = qrgr$tau
    res$globale_confiance = qrgr$globale_confiance
  }else{
    res$globale_confiance = sum(apply(Y - m(X) > 0, 1, all))/n
    res$tau_j = rep(tau, p)
  }


  # Y = Y * sapply(sens, function(x) ifelse(x == "min", -1, 1))
  # res$Y = res$Y * sapply(sens, function(x) ifelse(x == "min", -1, 1))
  Y = sweep(Y, 2, sapply(sens, function(x) ifelse(x == "min", -1, 1)), "*")
  res$Y = sweep(res$Y, 2, sapply(sens, function(x) ifelse(x == "min", -1, 1)), "*")

  res$Y0 = Y
  colnames(res$Y) = colnames(Y)
  res$tau = tau

  res
}
