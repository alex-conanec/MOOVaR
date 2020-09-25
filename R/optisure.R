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
optisure <- function(X, Y, sens = rep("min", NCOL(Y)), alpha = 0.5, N = NULL,
                     B = NULL, seed = NULL){

  cat("Recherche de la/les fenetre(s) mobile(s) optimale(s)\n")

  #anayse des factors dans X
  is_fac = sapply(X, is.factor)

  if (any(is_fac)){
    X_fac = X[,is_fac]
    combinaisons = unique(X_fac)
    # a checker quand combinaisons est nul!
    lapply(seq_len(NROW(combinaisons)), function(i){
      combi = combinaisons[i,]

      mask = which(apply(data.frame(
        X[, is_fac] == matrix(as.character(combi), ncol = 2, nrow = NROW(X), byrow = TRUE)
      ), 1, all))

      list(X_num = X[mask, !is_fac], X_fac = X[mask, is_fac], Y = Y[mask,],
           combi = combi)
    }) -> train_datas


    #calcul du h_n (des h_n si on a des var quali)
    lapply(seq_along(Y), function(j){
      cat(colnames(Y)[j], "\n")
      lapply(train_datas, function(data){
        cat("Combinaison",
            paste(colnames(X)[is_fac], data$combi, collapse = ' -- ', sep = ": "),
            "\n")
        hopt(X = data$X_num, Y = data$Y[,j])$h
      })
    }) -> all_h


    # les h tres importants s'explique par une tres faible correlation entre le X en question et le Y
    # combiné à une bonne correlation d'un des autres X.

    # library(FactoMineR)
    # PCA(X = cbind(train_datas[[3]]$X_num, train_datas[[3]]$Y))
    # cor(cbind(train_datas[[2]]$X_num, train_datas[[2]]$Y))

    #initiation de la fonction avec quantile
    fn = lapply(seq_len(p), function(j){
      function(X){
        lapply(seq_along(train_datas), function(q){

          #filtre par combinaison
          mask = which(apply(data.frame(
            X[, is_fac] == matrix(as.character(train_datas[[q]]$combi), ncol = 2, nrow = NROW(X), byrow = TRUE)
          ), 1, all))
          X_combi = X[mask,]

          #calcul le quantile pour chaque combinaison
          if (NROW(X_combi) > 0){
            data.frame(
              y = conditionnal_quantile(X = train_datas[[q]]$X_num,
                                        Y = train_datas[[q]]$Y[,j],
                                        x = X_combi[, !is_fac],
                                        alpha,
                                        h_n = all_h[[p]][[q]])$y,
              id = mask
            )
          }else NULL
        }) %>% bind_rows() %>% arrange(id) %>% pull(y)
      }
    })

  }else{ #si que des quanti

    h_n = lapply(seq_along(Y), function(j) hopt(X = X, Y = Y[,j])$h )
    fn = lapply(seq_len(p), function(j){
      function(x){
        conditionnal_quantile(X = X,
                              Y = Y[,j],
                              x = x,
                              alpha,
                              h_n = h_n[[j]])$y
      }
    })

  }

  names(fn) = paste0("Y", seq_along(Y))

  #tunage nsga parametres
  if (is.null(N)) N = 50
  if (is.null(B)) B = 20

  #appelle de NSGA ou optim si NCOL(Y)==1
  res = NSGA(X = X, fn, N = N, seed = seed, B = B, verbose = TRUE)

  #mise en forme des resultats
  res

}
