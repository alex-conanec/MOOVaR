#' @export
R2 = function(X, Y, B = 30, seed= NULL,
              allowed_dependence = matrix(TRUE, nrow = NCOL(X),
                                          ncol = NCOL(Y))){

  if (!is.null(seed)) set.seed(seed)

  p = NCOL(Y)
  lapply(1:B, function(i){
    print(i)
    lapply(1:p, function(j=1){
      mask = !is.na(Y[,j,T])
      XX=X[mask,]
      YY=Y[mask,j,F]
      n = NROW(XX)
      train_idx = sample(1:n, size = n*0.7)

      X_train = XX[train_idx,]
      Y_train = YY[train_idx,,F]
      X_test = XX[-train_idx,-1,F]
      Y_test = YY[-train_idx,,F]

      m_R2 = fit_model(X = X_train, Y = Y_train,
                       allowed_dependence = allowed_dependence[,j,drop=F],
                       method = "expected")

      X_train = X_train[,-1,F][,allowed_dependence[,j,drop=T],F]
      is_fac = !sapply(X_train, is.numeric)
      cstr_X_train_R2 = def_cstr_X_space(X_train[,is_fac,F])$g

      feasible = cstr_X_train_R2(X_test[,allowed_dependence[,j,drop=T],F][,is_fac,F])
      X_test = X_test[feasible,,F]
      Y_test = Y_test[feasible,,F]

      for (i in which(!allowed_dependence[,j,drop=T])){
        a = XX[,-1,F][,i,T]
        X_test[is.na(X_test[,i]),i] = a[!is.na(a)][1]
      }

      SCE_residual = colSums( (m_R2(X=X_test) - Y_test )^2)
      SCE_tot = colSums((matrix(apply(Y_test, 2, mean),
                                nrow = NROW(Y_test),
                                ncol = p, byrow = TRUE) -
                           Y_test)^2)
      as.data.frame(t(1 - SCE_residual/SCE_tot))
    }) %>% bind_cols()
  }) %>%
    bind_rows() %>% apply(2, mean)

}

