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
optisure <- function(X, Y, sens = rep("min", length(Y)),
                     g = NULL, X_space_csrt = TRUE,
                     tau = 0.5, globale_tau = FALSE, reg_method = "linear", #rajouté
                     # parametric = FALSE, deg = 3, #a virer
                     N = NULL, B = NULL, seed = NULL,
                     penalty = NULL){

  d = NCOL(X)
  p = NCOL(Y)
  n = NROW(X)
  col_names_Y = colnames(Y)
  if (is.null(penalty)){
    penalty = rep(NA, p)
  }else if (length(penalty) < p){
    stop("The penalty must be a vector of size p")
  }

  #sens quantile
  Y = Y * sapply(sens, function(x) ifelse(x == "min", -1, 1))

  fit_model = function(tau){
    if (length(tau) == 1 ) tau = rep(tau, p)
    if (reg_method == "linear"){
      X_mat = model.matrix(~., X)[,-1]
      beta_rq = sapply(seq_len(p), function(j){
        vide = capture.output(cv_rq <- hqreg::cv.hqreg(X_mat, Y[,j],
                                                       method = "quantile",
                                                       tau = tau[j],
                                                       alpha = 0,
                                                       ncores = 1, nfolds = 10, seed = 123))
        cv_rq$fit$beta[,which(cv_rq$lambda == cv_rq$lambda.1se)]
      })
      m = function(X) model.matrix(~., X) %*% beta_rq

    }else if (reg_method == "neural_network"){

      penalty_cv = function(X, y, tau){

        prop_train = 0.8
        train_idx = sample(seq_len(n), n*prop_train)
        test_X = X[-train_idx,]
        test_y = y[-train_idx]
        train_X = X[train_idx,]
        train_y = y[train_idx]

        MSE = function(penalty){
          vide = capture.output(m <- qrnn::qrnn2.fit(x = as.matrix(train_X),
                          y = as.matrix(train_y),
                          n.hidden=2, n.hidden2=2,
                          tau = tau, penalty=penalty))

          sum((qrnn::qrnn2.predict(as.matrix(test_X), m) - test_y)^2)
        }

        optimize(f = MSE, lower = 0, upper = 1)$minimum

      }

      for (j in seq_len(p)){
        if (is.na(penalty[j])){
          cat("cross-validate determining ridge penalty of the neural network",
              j, "\n")
          penalty[j] = penalty_cv(X, Y[,j], tau)
        }
      }

      cat("train neural network \n")
      m_j = lapply(seq_len(p), function(j){
        vide = capture.output(m <- qrnn::qrnn2.fit(
          x = as.matrix(X), y = as.matrix(Y[,j, drop=FALSE]),
          n.hidden=2, n.hidden2=2, tau = tau[j], penalty = penalty[j]))
        m
      })

      m = function(X){

        #be sure that NROW(X) > 1 sinon ca fou le bordel
        res = lapply(m_j, function(m){
          qrnn::qrnn2.predict(as.matrix(X), m)
        }) %>% as.data.frame()
        colnames(res) = colnames(Y)
        res
      }

    }

    m
  }

  quant_reg_global_risk = function(X, Y, tau, plafond = 0.95, max_iter = 20){

    d = NCOL(X_mat)
    lambda = rep(0, p)
    lambda_prec = rep(1, p)


    for (it in seq_len(max_iter)){
      m = fit_model(tau + lambda)

      #eval
      delta = Y - m(X)
      cat(colSums(delta < 0)/n, "\n")
      L = sum(apply(delta < 0, 1, all))/n
      cat(L, "\n")

      #evolution lambda
      S = sapply(seq_len(p), function(j){
        sum(apply(delta[, -j] < 0, 1, all))/n
      })
      S = S/sum(S)
      lambda = lambda + (tau - L) * S

      #avoid learn too high repartition = avoid estimation pb
      lambda[lambda + tau > plafond] = plafond - tau


      if (all(lambda == lambda_prec)){
        break
      }
      lambda_prec = lambda

    }


    list(model = m, alpha = L, tau = tau + lambda)
  }

  if (globale_tau){
    qrgr = quant_reg_global_risk(X = X_mat, Y, tau)
    m = qrgr$m
  }else{
    m = fit_model(tau)
    m(X)
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
  Y = Y * sapply(sens, function(x) ifelse(x == "min", -1, 1))
  res$Y = res$Y * sapply(sens, function(x) ifelse(x == "min", -1, 1))
  res$Y0 = Y
  colnames(res$Y) = colnames(Y)
  res$tau = tau

  res
}
