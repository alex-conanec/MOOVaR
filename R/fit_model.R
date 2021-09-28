#' fit_model
#'
#' @param X matrix or data.frame of n observation of the decision variable of size d.
#' @param Y matrix or data.frame of n observation of the p objectifs.
#' @param reg_method only "linear" available yet.
#' @param tau vector of size NCOL(Y) containing the risk in ]0,1[ of each objectif.
#' @param method either "quantile" or "expected"
#' @param penalty NULL by default

#' @param path_tracking path where to write the step of the running function.

#' @export
fit_model = function(X, Y, allowed_dependence, reg_method = "linear", tau = NULL,
                     method = "quantile", penalty = NULL,
                     path_tracking = NULL){


  #tracking stop
  is_it_stopped(path_tracking)

  names_X_mat = colnames(model.matrix(~., X[,-1,F])[,-1, drop = F])

  p = NCOL(Y)
  if (length(tau) == 1) tau = rep(tau, p)
  is_fac = !sapply(X[,-1], is.numeric)

  X_list = lapply(1:p, function(j){
    # print(j)

    mask_dupli = !duplicated(X[,c(TRUE, allowed_dependence[,j]), drop = FALSE])
    id = X[mask_dupli,1,T]
    X = X[mask_dupli,-1,F]
    X_quanti = scale(X[, !is_fac, drop=F])
    mu = attr(X_quanti, "scaled:center")
    sd = attr(X_quanti, "scaled:scale")



    #NaN introduit car tout a zero donc sd a zero donc div a zero. On remet zero le tout
    for (i in which(sd == 0)){
      X_quanti[,i]=mu[i]
      sd[i]=1 #sinon divise par zero dans fonction m!!
    }

    X[, !is_fac] = X_quanti
    X = X[,allowed_dependence[,j]]
    mask_na_yj = !is.na(Y[mask_dupli,j,T])
    id = id[mask_na_yj]
    X_mat = model.matrix(~., X[mask_na_yj,])[,-1, drop = F]
    # X_mat = X_mat[mask_na_yj,]

    # print(NROW(X_mat))
    list(
      X_mat = X_mat,
      mu = mu,
      sd = sd,
      mask_dupli = mask_dupli,
      mask_na_yj = mask_na_yj,
      id = id
    )

  })


  Y_list = lapply(1:p, function(j){
    y = scale(Y[X_list[[j]]$mask_dupli , j, T][X_list[[j]]$mask_na_yj])
    list(
      y = y[,1],
      mu = attr(y, "scaled:center"),
      sd = attr(y, "scaled:scale")
    )
  })
  names(Y_list) = colnames(Y)

  X = X[,-1]
  allowed_dependence = lapply(1:NCOL(X), function(i=1){
    if (is.numeric(X[,i,TRUE])){
      matrix(allowed_dependence[i,], nrow = 1,
             ncol = NCOL(allowed_dependence), byrow = TRUE) %>%
        as.data.frame()
    }else{
      matrix(allowed_dependence[i,], nrow = nlevels(as.factor(X[,i, drop = TRUE]))-1,
             ncol = NCOL(allowed_dependence), byrow = TRUE) %>%
        as.data.frame()
    }
  }) %>% bind_rows()

  if (method == "quantile"){
    if (reg_method == "linear"){
      beta = sapply(seq_len(p), function(j){
        vide = capture.output(cv_rq <- hqreg::cv.hqreg(
          X = X_list[[j]]$X_mat,
          y = Y_list[[j]]$y,
          method = "quantile",
          tau = tau[j],
          alpha = 0,
          ncores = 1, nfolds = 10, seed = 123))

        beta = cv_rq$fit$beta[,which(cv_rq$lambda == cv_rq$lambda.1se)]

        beta_res = rep(0, NROW(allowed_dependence) + 1)
        beta_res[c(TRUE, allowed_dependence[,j])] = beta
        beta_res
      })

      rownames(beta) = c("intercept", names_X_mat)
      colnames(beta) = colnames(Y)

      m = function(X, X_list, Y_list, beta){

        # X_list = formals(m)$X_list
        # Y_list = formals(m)$Y_list
        # beta = formals(m)$beta

        is_fac = !sapply(X, is.numeric)
        p = NCOL(beta)

        X[,!is_fac] = scale_x(X[,!is_fac], mu = X_list[[1]]$mu,
                              s = X_list[[1]]$sd)

        X_mat = model.matrix(~., as.data.frame(X))

        scale_x(X_mat %*% beta,
                mu = sapply(Y_list, function(y) y$mu),
                s = sapply(Y_list, function(y) y$sd),
                reverse = TRUE)
      }

    }else if (reg_method == "neural_network"){
      #pass
    }
  }else if (method == "expected"){
    if (reg_method == "linear"){
      # X_mat = cbind(1, X_mat)
      allowed_dependence = rbind(TRUE, allowed_dependence)
      beta = sapply(seq_len(p), function(j=1){
        cv = glmnet::cv.glmnet(
          x = cbind(1, X_list[[j]]$X_mat),
          y = Y_list[[j]]$y,
          alpha = 0,
          intercept = TRUE)

        lambda = cv$lambda.1se

        model = glmnet::glmnet(
          x = X_list[[j]]$X_mat,
          y = Y_list[[j]]$y,
          lambda = lambda,
          alpha = 0,
          intercept = TRUE)

        beta = c(model$a0, model$beta[,1])

        beta_res = rep(0, NROW(allowed_dependence))
        beta_res[allowed_dependence[,j]] = beta
        beta_res
      })
      rownames(beta) = c("intercept", names_X_mat)
      colnames(beta) = colnames(Y)

      m = function(X, X_list, Y_list, beta){

        # X_list = formals(m_R2)$X_list
        # Y_list = formals(m_R2)$Y_list
        # beta = formals(m_R2)$beta

        is_fac = !sapply(X, is.numeric)
        p = NCOL(beta)

        X[,!is_fac] = scale_x(X[,!is_fac], mu = X_list[[1]]$mu,
                              s = X_list[[1]]$sd)

        X_mat = model.matrix(~., as.data.frame(X))

        scale_x(X_mat %*% beta,
                mu = sapply(Y_list, function(y) y$mu),
                s = sapply(Y_list, function(y) y$sd),
                reverse = TRUE)
      }

    }else if (reg_method == "neural_network"){
      # pass
    }
  }
  formals(m)$X_list = X_list
  formals(m)$Y_list = Y_list
  formals(m)$beta = beta
  m
}

# penalty_cv = function(X, y, tau){
#
#   prop_train = 0.8
#   n = NROW(X)
#   train_idx = sample(seq_len(n), n*prop_train)
#   X = model.matrix(~., as.data.frame(X))[,-1]
#   test_X = scale(X[-train_idx,])
#   test_y = scale(y[-train_idx])
#   train_X = scale(X[train_idx,])
#   train_y = scale(y[train_idx])
#
#   rho = function(u){
#     tau*u*as.numeric(u>=0) - (1-tau) * u * as.numeric(u<0)
#   }
#
#   loss = function(penalty){
#     vide = capture.output(m <- qrnn::qrnn2.fit(x = as.matrix(train_X),
#                                                y = as.matrix(train_y),
#                                                n.hidden=2, n.hidden2=2,
#                                                tau = tau, penalty=penalty))
#
#     # sum((qrnn::qrnn2.predict(as.matrix(test_X), m) - test_y)^2)
#     sum(rho(qrnn::qrnn2.predict(as.matrix(test_X), m) - test_y))
#   }
#
#   # pena = seq(0, 0.35, 0.05)
#   # res_pena = parallel::mclapply(pena, loss, mc.cores = 4)
#   # plot(pena, res_pena, type = 'b')
#
#   optimize(f = loss, lower = 0, upper = .4)$minimum
#
# }

#' @export
scale_x = function(x, mu, s, reverse = FALSE){

  x = as.matrix(x)
  if (reverse){
    as.matrix(sweep(sweep(x, 2, s, "*"), 2, mu, "+"))
  }else{
    as.matrix(sweep(sweep(x, 2, mu, "-"), 2, s, "/"))
  }

}
