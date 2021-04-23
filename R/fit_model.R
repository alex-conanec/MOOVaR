#' @export
fit_model = function(X, Y, reg_method, tau = NULL, method = "quantile", penalty = NULL){
  p = NCOL(Y)
  X_mat = model.matrix(~., X)[,-1]
  if (length(tau) == 1 ) tau = rep(tau, p)

  #normalisation
  X_mat = scale(X_mat)
  X_mat_scale = list(mu = attr(X_mat, "scaled:center"),
                     sd = attr(X_mat, "scaled:scale"))
  Y_mat = scale(Y)
  Y_mat_scale = list(mu = attr(Y_mat, "scaled:center"),
                     sd = attr(Y_mat, "scaled:scale"))

  if (method == "quantile"){
    if (reg_method == "linear"){
      beta_rq = sapply(seq_len(p), function(j){
        vide = capture.output(cv_rq <- hqreg::cv.hqreg(X_mat, Y_mat[,j],
                                                       method = "quantile",
                                                       tau = tau[j],
                                                       alpha = 0,
                                                       ncores = 1, nfolds = 10, seed = 123))
        # cv_rq$fit$beta[,which(cv_rq$lambda == cv_rq$lambda.1se)]
        hqreg::hqreg(X_mat, Y_mat[,j], method = "quantile", tau = tau[j],
                     alpha = 0, lambda = 0)$beta
      })

      m = function(X, X_mat_scale, Y_mat_scale){
        X_mat = model.matrix(~., as.data.frame(X))[,-1]
        X_mat = scale_x(X_mat, mu = X_mat_scale$mu, s = X_mat_scale$sd)

        y_pred = cbind(1, X_mat) %*% beta_rq
        scale_x(y_pred, mu = Y_mat_scale$mu, s = Y_mat_scale$sd, reverse = TRUE)
      }

    }else if (reg_method == "neural_network"){

      if (is.null(penalty)){
        penalty = rep(NA, p)
      }else if (length(penalty) < p){
        stop("The penalty must be a vector of size p")
      }


      #penalty tuning
      core_max_lim = 2
      if (.Platform$OS.type == "unix"){
        penalty = parallel::mclapply(seq_len(p), function(j){
          if (is.na(penalty[j])){
            cat("cross-validate determining ridge penalty of the neural network",
                j, "\n")
            penalty_cv(X = X_mat, y = Y_mat[,j], tau = tau[j])
          }else{
            penalty[j]
          }
        }, mc.cores = min(core_max_lim, p)) %>% unlist()
      }else{
        for (j in seq_len(p)){
          if (is.na(penalty[j])){
            cat("cross-validate determining ridge penalty of the neural network",
                j, "\n")
            penalty[j] = penalty_cv(X = X_mat, y = Y_mat[,j], tau = tau[j])
          }
        }
      }

      cat("train neural network \n")
      m_j = lapply(seq_len(p), function(j){
        vide = capture.output(m <- qrnn::qrnn2.fit(
          x = X_mat, y = as.matrix(Y_mat[,j, drop=FALSE]),
          n.hidden=2, n.hidden2=2, tau = tau[j], penalty = penalty[j]))
        m
      })

      m = function(X, X_mat_scale, Y_mat_scale, m_j){

        X_mat = model.matrix(~., as.data.frame(X))[,-1]
        X_mat = scale_x(X_mat, mu = X_mat_scale$mu, s = X_mat_scale$sd)


        #be sure that NROW(X) > 1 sinon ca fou le bordel
        y_pred = lapply(m_j, function(m){
          qrnn::qrnn2.predict(X_mat, m)
        }) %>% as.data.frame()
        colnames(y_pred) = colnames(Y)
        scale_x(y_pred, mu = Y_mat_scale$mu, s = Y_mat_scale$sd, reverse = TRUE)
      }
      formals(m)$m_j = m_j

    }
  }else if (method == "expected"){
    if (reg_method == "linear"){
      X_mat = cbind(1, X_mat)
      beta_re = sapply(seq_len(p), function(j=1){
        # cv = glmnet::cv.glmnet(x = X_mat, y = Y_mat[,j])
        # lambda = cv$lambda.1se
        lambda = 0
        glmnet::glmnet(x = X_mat, y = Y_mat[,j], lambda = lambda)$beta[,1]
        # solve(t(X_mat)%*%X_mat)%*%t(X_mat)%*%Y_mat[,j]
      })
      m = function(X, X_mat_scale, Y_mat_scale){
        X_mat = model.matrix(~., as.data.frame(X))[,-1]
        X_mat = scale_x(X_mat, mu = X_mat_scale$mu, s = X_mat_scale$sd)
        X_mat = cbind(1, X_mat)
        y_pred = X_mat %*% beta_re
        scale_x(y_pred, mu = Y_mat_scale$mu, s = Y_mat_scale$sd, reverse = TRUE)
      }

    }else if (reg_method == "neural_network"){

      ANN_model = function(){

        library(keras)
        inputs <- keras::layer_input(shape = d, name = 'aux_input')

        outputs = inputs %>%
          keras::layer_dense(units = 4, activation = 'relu') %>%
          layer_dropout(rate = 0.5) %>%
          keras::layer_dense(units = 4, activation = "relu") %>%
          layer_dropout(rate = 0.5) %>%
          keras::layer_dense(units = 1)

        model <- keras::keras_model(inputs = inputs, outputs = outputs)

        model %>% keras::compile(
          optimizer = keras::optimizer_rmsprop(),
          loss = "mse",
          metrics = c("mean_absolute_error")
        )

        cat("ANN_built \n")
        model

      }

      cat("train neural network \n")
      d = NCOL(X_mat)
      m_j = lapply(seq_len(p), function(j){
        model = ANN_model()

        model %>% keras::fit(x = list(X_mat), y = list(Y_mat[,j]),
                             batch_size  = n/10, epochs = 500,
                             validation_split = 0.2, shuffle = T, verbose = 0)
        model
      })

      m = function(X, X_mat_scale, Y_mat_scale, m_j){
        X_mat = model.matrix(~., as.data.frame(X))[,-1]
        X_mat = scale_x(X_mat, mu = X_mat_scale$mu, s = X_mat_scale$sd)

        y_pred = lapply(m_j, function(model){
          predict(model, X_mat)
        }) %>% as.data.frame()
        colnames(y_pred) = colnames(Y)
        scale_x(y_pred, mu = Y_mat_scale$mu, s = Y_mat_scale$sd, reverse = TRUE)
      }
      formals(m)$m_j = m_j
    }
  }
  formals(m)$X_mat_scale = X_mat_scale
  formals(m)$Y_mat_scale = Y_mat_scale
  m
}

#' @export
penalty_cv = function(X, y, tau){

  prop_train = 0.8
  n = NROW(X)
  train_idx = sample(seq_len(n), n*prop_train)
  X = model.matrix(~., as.data.frame(X))[,-1]
  test_X = scale(X[-train_idx,])
  test_y = scale(y[-train_idx])
  train_X = scale(X[train_idx,])
  train_y = scale(y[train_idx])

  rho = function(u){
    tau*u*as.numeric(u>=0) - (1-tau) * u * as.numeric(u<0)
  }

  loss = function(penalty){
    vide = capture.output(m <- qrnn::qrnn2.fit(x = as.matrix(train_X),
                                               y = as.matrix(train_y),
                                               n.hidden=2, n.hidden2=2,
                                               tau = tau, penalty=penalty))

    # sum((qrnn::qrnn2.predict(as.matrix(test_X), m) - test_y)^2)
    sum(rho(qrnn::qrnn2.predict(as.matrix(test_X), m) - test_y))
  }

  # pena = seq(0, 0.35, 0.05)
  # res_pena = parallel::mclapply(pena, loss, mc.cores = 4)
  # plot(pena, res_pena, type = 'b')

  optimize(f = loss, lower = 0, upper = .4)$minimum

}

#' @export
scale_x = function(x, mu, s, reverse = FALSE){

  x = as.matrix(x)
  if (reverse){
    as.matrix(sweep(sweep(x, 2, s, "*"), 2, mu, "+"))
  }else{
    as.matrix(sweep(sweep(x, 2, mu, "-"), 2, s, "/"))
  }

}
