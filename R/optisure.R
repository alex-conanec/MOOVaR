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
#' res = optisure(X, Y, N=20, TT=10, alpha = 0.5)
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
#' res = optisure(X, Y, N=20, TT=10, alpha = 0.5)
#'
#' plot(ellipse_var(Y$Y1, Y$Y2, alpha = 0.05), type = 'l', xlab = "Y1", ylab = "Y2",
#'      xlim = c(-8,8), ylim = c(-10,6))
#' points(Y, col = "grey")
#' points(res$Y, type = 'b', col = "blue", lwd = 3)
#'
#' @importFrom magrittr %>%
#' @export
optisure <- function(X, Y, sens = rep("max", NCOL(Y)),
                     quantile_utility_idx = rep(TRUE, NCOL(Y)),
                     g = NULL, X_space_csrt = TRUE,
                     tau = rep(0.5, NCOL(Y)),
                     globale_tau = rep(F, NCOL(Y)),
                     alpha = 0.15,
                     # reg_method = "linear", #rajouté
                     N = NULL, TT = NULL, seed = NULL,
                     updateProgress = NULL,
                     path_tracking = NULL,
                     ...){

  # library(reticulate)
  d = NCOL(X)
  p = NCOL(Y)
  n = NROW(X)
  col_names_Y = colnames(Y)
  for (i in which(!sapply(X, is.numeric))) X[,i] = as.factor(X[,i,drop=TRUE])

  #sens optimisation
  Y = as.data.frame(Y)
  Y = sweep(Y, 2, sapply(sens, function(x) ifelse(x == "min", -1, 1)), "*")

  #tracking
  if (is.function(updateProgress)) {
    updateProgress(detail = "Training model")
  }

  tracking_msg(path_tracking, msg = "Training model")

  if (any(quantile_utility_idx)){
    if (!any(quantile_utility_idx[!globale_tau])){#global risk management only
      glob_tau = unique(tau[quantile_utility_idx])
      if (length(glob_tau) > 1){
        stop("tau must be the same for objectif treating in globale risk")
      }
      qrgr = quant_reg_global_risk(X = X, Y[,quantile_utility_idx, drop = FALSE],
                                   tau = glob_tau, path_tracking = path_tracking)
      m_quantile = qrgr$m
    }else if (!any(quantile_utility_idx[globale_tau])){#individual risk management only
      m_quantile = fit_model(X, Y = Y[,quantile_utility_idx, drop = FALSE],
                             tau = tau[quantile_utility_idx],
                             method = "quantile", path_tracking = path_tracking)
    }else{#gobal AND individual risk management
      #global risk
      global_risk_idx = quantile_utility_idx & globale_tau
      glob_tau = unique(tau[global_risk_idx])
      if (length(glob_tau) > 1){
        stop("tau must be the same for objectif treating in globale risk")
      }
      qrgr = quant_reg_global_risk(X = X, Y = Y[,global_risk_idx, drop = FALSE],
                                   tau = glob_tau, path_tracking = path_tracking)
      m_quantile_glob = qrgr$m

      #individual risk
      ind_risk_idx = quantile_utility_idx & !globale_tau
      m_quantile_ind = fit_model(X, Y = Y[,ind_risk_idx, drop = FALSE],
                                 tau = tau[ind_risk_idx],
                                 method = "quantile",
                                 path_tracking = path_tracking)
    }
  }

  if (any(!quantile_utility_idx)){
    m_expect = fit_model(X, Y = Y[,!quantile_utility_idx, drop = FALSE],
                         method = "expected", path_tracking = path_tracking)
  }

  if (all(quantile_utility_idx)){
    if (any(quantile_utility_idx[globale_tau]) & any(!quantile_utility_idx[globale_tau])){
      global_quantile = quantile_utility_idx & globale_tau
      m = function(X){
        res = cbind(m_quantile_glob(X), m_quantile_ind(X))[, c(which(global_quantile),
                                                               which(!global_quantile))]
        colnames(res) = colnames(Y)
        res
      }
    }else{
      m = m_quantile
    }
  }else if (all(!quantile_utility_idx)) {
    m = m_expect
  }else{
    if (any(quantile_utility_idx[globale_tau]) & any(quantile_utility_idx[!globale_tau])){
      global_quantile = quantile_utility_idx & globale_tau
      ind_quantile = quantile_utility_idx & !globale_tau

      m = function(X){
        res = cbind(m_quantile_glob(X), m_quantile_ind(X), m_expect(X))[
          , c(which(global_quantile),
              which(ind_quantile),
              which(!quantile_utility_idx))]
        colnames(res) = colnames(Y)
        res
      }
    }else{
      m = function(X){
        res = cbind(m_quantile(X), m_expect(X))[, c(which(quantile_utility_idx),
                                                    which(!quantile_utility_idx))]
        colnames(res) = colnames(Y)
        res
      }
    }

  }


  #tunage nsga parametres
  if (is.null(TT)) TT = 20
  if (is.null(N)) N = round(NROW(X)/2)

  if (X_space_csrt){
    #tracking
    if (is.function(updateProgress)) {
      updateProgress(detail = "Decision constraint training")
    }
    tracking_msg(path_tracking, msg = "Decision constraint training")

    cat("Train belonging to density constraint \n")
    X_space_csrt = def_cstr_X_space(X, alpha = alpha,
                                    path_tracking = path_tracking)
    g$cstr_X_space = X_space_csrt$g
  }

  if (any(sapply(X, is.numeric))){
    mutation_method = "mixte"
  }else{
    mutation_method = "simple"
  }

  #appelle de NSGA
  res = NSGA(
    X = X,
    fn = m,
    n_objective = p,
    sens = rep("max", p),
    g = g,

    #k tournois
    k_tournament = 10,

    #crossing over
    crossing_over_method = "uniforme",
    n_echange = 2,
    n_bloc = 2,

    #mutation
    mutation_method = mutation_method,
    freq_m = 0.3,

    #budget
    N = N,
    TT = TT,

    #verbose
    verbose = TRUE,
    front_tracking = TRUE,
    seed = seed,
    updateProgress = updateProgress,
    path_tracking = path_tracking
  )

  if (is.function(updateProgress)) {
    updateProgress(detail = "end")
  }
  tracking_msg(path_tracking, msg = "End")

  #mise en forme des resultats ##### a modifier aussi
  # if (globale_tau){
  #   res$tau_j = qrgr$tau
  #   res$globale_confiance = qrgr$globale_confiance
  # }else{
  #   res$globale_confiance = sum(apply(Y - m(X) > 0, 1, all))/n
  #   res$tau_j = rep(tau, p)
  # }


  # Y = Y * sapply(sens, function(x) ifelse(x == "min", -1, 1))
  # res$Y = res$Y * sapply(sens, function(x) ifelse(x == "min", -1, 1))
  Y = sweep(Y, 2, sapply(sens, function(x) ifelse(x == "min", -1, 1)), "*")
  res$Y = sweep(res$Y, 2, sapply(sens, function(x) ifelse(x == "min", -1, 1)), "*")

  res$Y0 = Y
  colnames(res$Y) = colnames(Y)
  res$tau = tau
  res$X_space_csrt = X_space_csrt

  res
}
