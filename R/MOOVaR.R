#' MOOVaR
#'
#' @param X matrix or data.frame of n observation of the decision variable of size d.
#' @param Y matrix or data.frame of n observation of the p objectifs.
#' @param sens vector of size p containing either "max" (by default) or "min"
#'  to choose how to optimize each objectif.
#' @param quantile_utility_idx vector of size p containing boolean to choose the utility method.
#' TRUE and therefore quantile method is used by default.
#' @param g list of constraint given as function of X. NULL by default
#' @param X_space_csrt boolean to choose either a constraint should be added to
#'  check if each solution belongs to the decision space.
#' @param tau vector of size p containing the risk in ]0,1[ of each objectif.
#' @param globale_tau vector of size p boolean indicating either the objectif risk must ne manage with other.
#' @param alpha probability to deal with the trade-off false positive false
#'  negatif. Only if X_space_csrt=TRUE.
#' @param updateProgress function to follow the progression of the running function
#' @param path_tracking path where to write the step of the running function.
#'
#'
#' @examples
#' library(tidyverse)
#' library(MOOVaR)
#' set.seed(12)
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
#' p = 4
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
#'
#' g_cstr = list(
#'   function(x){
#'     x = as.data.frame(x)
#'     x[,1] + x[,2] + x[,3] < 2
#'   }
#' )
#'
#'
#' res = MOOVaR(X, Y,
#'              sens = rep("max", NCOL(Y)),
#'              quantile_utility_idx = c(T,T,T,F),
#'              tau = rep(0.2, p),
#'              globale_tau = c(F, T, T, F),
#'              g = g_cstr,
#'              X_space_csrt=T, TT = 10)
#'
#'
#' plot(res)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @export
MOOVaR <- function(X, Y, sens = rep("max", NCOL(Y)),
                   quantile_utility_idx = rep(TRUE, NCOL(Y)),
                   optim_method = c("real_ind", "nsga", "none"),
                   g = NULL, X_space_csrt = FALSE,
                   tau = rep(0.5, NCOL(Y)),
                   globale_tau = rep(F, NCOL(Y)),
                   alpha = 0.15,
                   updateProgress = NULL,
                   path_tracking = NULL,
                   mutation_method = "simple",
                   allowed_dependence = matrix(TRUE, nrow = NCOL(X), ncol = NCOL(Y)),
                   seed_R2 = NULL,
                   ...){

  # library(reticulate)
  d = NCOL(X)-1
  p = NCOL(Y)
  n = NROW(X)
  col_names_Y = colnames(Y)
  is_fac = !sapply(X, is.numeric)
  for (i in which(is_fac)) X[,i] = as.factor(X[,i,drop=TRUE])

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
      qrgr = quant_reg_global_risk(X = X, Y=Y[,quantile_utility_idx, drop = FALSE],
                                   allowed_dependence = allowed_dependence[,quantile_utility_idx, drop = FALSE],
                                   tau = glob_tau, path_tracking = path_tracking)
      m_quantile = qrgr$m
    }else if (!any(quantile_utility_idx[globale_tau]) |
              sum(quantile_utility_idx[globale_tau]) == 1){#individual risk management only
      m_quantile = fit_model(X, Y = Y[,quantile_utility_idx, drop = FALSE],
                             tau = tau[quantile_utility_idx],
                             allowed_dependence = allowed_dependence[,quantile_utility_idx, drop = FALSE],
                             method = "quantile", path_tracking = path_tracking)
    }else{#gobal AND individual risk management
      #global risk
      global_risk_idx = quantile_utility_idx & globale_tau
      glob_tau = unique(tau[global_risk_idx])
      if (length(glob_tau) > 1){
        stop("tau must be the same for objectif treating in globale risk")
      }
      qrgr = quant_reg_global_risk(X = X, Y = Y[,global_risk_idx, drop = FALSE],
                                   allowed_dependence = allowed_dependence[,global_risk_idx, drop = FALSE],
                                   tau = glob_tau, path_tracking = path_tracking)
      m_quantile_glob = qrgr$m

      #individual risk
      ind_risk_idx = quantile_utility_idx & !globale_tau
      m_quantile_ind = fit_model(X, Y = Y[,ind_risk_idx, drop = FALSE],
                                 tau = tau[ind_risk_idx],
                                 allowed_dependence = allowed_dependence[,ind_risk_idx, drop = FALSE],
                                 method = "quantile",
                                 path_tracking = path_tracking)
    }
  }

  if (any(!quantile_utility_idx)){
    m_expect = fit_model(X, Y = Y[,!quantile_utility_idx, drop = FALSE],
                         allowed_dependence = allowed_dependence[,!quantile_utility_idx, drop = FALSE],
                         method = "expected", path_tracking = path_tracking)
  }

  if (all(quantile_utility_idx)){
    if (any(quantile_utility_idx[globale_tau]) & any(quantile_utility_idx[!globale_tau])){
      m = function(X){
        cbind(m_quantile_glob(X), m_quantile_ind(X)) %>%
          as.data.frame() %>% select(!!colnames(Y))
      }
      beta = cbind(formals(m_quantile_glob)$beta,
                   formals(m_quantile_ind)$beta)
    }else{
      m = m_quantile
      beta = formals(m)$beta
    }
  }else if (all(!quantile_utility_idx)) {
    m = m_expect
    beta = formals(m)$beta
  }else{
    if (any(quantile_utility_idx[globale_tau]) & any(quantile_utility_idx[!globale_tau])){
      global_quantile = quantile_utility_idx & globale_tau
      ind_quantile = quantile_utility_idx & !globale_tau

      m = function(X){
        cbind(m_quantile_glob(X), m_quantile_ind(X), m_expect(X)) %>%
          as.data.frame() %>% select(!!colnames(Y))
      }

      beta = cbind(formals(m_quantile_glob)$beta,
                   formals(m_quantile_ind)$beta,
                   formals(m_expect)$beta)

    }else{
      m = function(X){
        cbind(m_quantile(X), m_expect(X)) %>%
          as.data.frame() %>% select(!!colnames(Y))
      }

      beta = cbind(formals(m_quantile)$beta,
                   formals(m_expect)$beta)
    }

  }


  mask = apply(!is.na(Y), 1, all)
  X = X[mask,]
  Y = Y[mask,]

  summary(X)
  #remove id
  X = X[,-1,F]

  if (optim_method == "real_ind"){
    time_deb = Sys.time()
    Y_obj <- m(X) %>% as.data.frame() %>% mutate(id = 1:NROW(.))
    Y_obj$rank = dominance_ranking(Y_obj[,-NCOL(Y_obj)], rep("max", p))
    Y_obj = crowding_distance(Y_obj)
    Y_obj = Y_obj %>% filter(rank == 1, crowding_distance > 0)

    res = list(X = X[Y_obj$id,], Y = Y_obj[,1:p], time = Sys.time() - time_deb,
               fn = m, g = g, X0 = X, all_front = NULL, sens = rep("max", p))
    class(res) = "nsga"

    Y = sweep(Y, 2, sapply(sens, function(x) ifelse(x == "min", -1, 1)), "*")
    res$Y = sweep(res$Y, 2, sapply(sens, function(x) ifelse(x == "min", -1, 1)), "*")
    res$Y0 = Y
    colnames(res$Y) = colnames(Y)

  }else if (optim_method == "nsga"){

    if (X_space_csrt){
      #tracking
      if (is.function(updateProgress)) {
        updateProgress(detail = "Decision constraint training")
      }
      tracking_msg(path_tracking, msg = "Decision constraint training")

      cat("Train belonging to density constraint \n")
      cstr_X_space = def_cstr_X_space(X, alpha = alpha,
                                      path_tracking = path_tracking)
      g$cstr_X_space = cstr_X_space$g
    }


    if (any(sapply(X, is.numeric))){
      mutation_method = "mixte"
    }else{
      mutation_method = "simple"
    }

    res = NSGA(
      X = X,
      fn = m,
      n_objective = p,
      sens = rep("max", p),
      g = g,
      mutation_method = mutation_method,
      updateProgress = updateProgress,
      path_tracking = path_tracking,
      ...
    )
    Y = sweep(Y, 2, sapply(sens, function(x) ifelse(x == "min", -1, 1)), "*")
    res$Y = sweep(res$Y, 2, sapply(sens, function(x) ifelse(x == "min", -1, 1)), "*")
    res$Y0 = Y
    colnames(res$Y) = colnames(Y)
  }else{
    warning("optim_method must be in c('nsga', 'real_ind'), no optimisation done")
    res=list()
  }

  if (is.function(updateProgress)) {
    updateProgress(detail = "end")
  }
  tracking_msg(path_tracking, msg = "End")

  #mise en forme des resultats
  if (any(quantile_utility_idx)){
    res$tau_j = rep("-", p)
    res$tau_j = tau[quantile_utility_idx]
    mask = apply(!is.na(Y), 1, all)
    res$globale_confiance = rep("-", p)
    if (NROW(X)>0){
      res$globale_confiance[quantile_utility_idx] =
        rep(sum(apply( Y[mask,quantile_utility_idx] -
                        m(X[mask,])[,quantile_utility_idx] > 0, 1, all))/n, sum(quantile_utility_idx))
    }
    if (sum(quantile_utility_idx[globale_tau]) >= 2){
      res$tau_j[quantile_utility_idx & globale_tau] = round(qrgr$tau, 3)
      res$globale_confiance[quantile_utility_idx & globale_tau] = round(qrgr$globale_confiance, 3)
    }
  }else{
    res$tau_j = rep("-", p)
    res$globale_confiance = rep("-", p)
  }

  res$tau = tau
  # res$X_space_csrt = cstr_X_space
  res$beta = beta
  res$m = m
  # res$R2 = R2

  res
}
