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
optisure <- function(X, Y, fn, sens = rep("min", length(Y)),
                     g = NULL, X_space_csrt = TRUE, parametric = FALSE,
                     alpha = 0.5, deg = 3, N = NULL, B = NULL, seed = NULL){

  fn_to_opt = fn
  d = NCOL(X)
  p = NCOL(Y)
  col_names_Y = colnames(Y)

  alpha_j = sapply(sens, function(x) if (x == "min") alpha else 1 - alpha)

  X = as.data.frame(X)
  Y = as.data.frame(Y)

  if (!parametric){
    cat("Recherche de la/les fenetre(s) mobile(s) optimale(s)\n")
  }else{
    cat("fit beta\n")
  }


  #anayse des factors dans X
  is_fac = !sapply(X, is.numeric)

  if (any(is_fac)){
    X_fac = X[,is_fac]
    for (k in seq_along(X_fac)) X_fac[,k] = as.character(X_fac[,k])
    combinaisons = unique(X_fac)

    # a checker quand combinaisons est nul!
    lapply(seq_len(NROW(combinaisons)), function(i){
      combi = combinaisons[i,]
      mask = which(apply(data.frame(
        X[, is_fac] == as.data.frame(matrix(combi, ncol = 2, nrow = NROW(X), byrow = TRUE))
      ), 1, all))

      list(X_num = X[mask, !is_fac],
           # X_fac = X[mask, is_fac],
           Y = Y[mask,],
           combi = combi)
    }) -> train_datas


    if (!parametric){
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
    }else{
      #calcul du h_n (des h_n si on a des var quali)
      lapply(seq_along(Y), function(j){
        cat(colnames(Y)[j], "\n")
        lapply(train_datas, function(data){
          cat("Combinaison",
              paste(colnames(X)[is_fac], data$combi, collapse = ' -- ', sep = ": "),
              "\n")
          fit_beta(X = data$X_num, Y = data$Y[,j], deg = deg, alpha = alpha_j[j])
        })
      }) -> all_beta
    }

    # function to optimize in NSGA
    fn = function(x){

      if (NCOL(x) < d){
        x = matrix(x, ncol = d)
      }
      x = data.frame(x)
      n = NROW(x)

      #on cherche les combinaison de x
      combi_x = unique(x[,is_fac])

      #pour toutes les combinaison de x
      Y = lapply(seq_len(NROW(combi_x)), function(k){
        mask = sapply(seq_len(n), function(i){
          all(x[i, is_fac, drop = FALSE] == combi_x[k, ,drop = FALSE])
        }) %>% which()

        train_datas_idx_k = sapply(seq_along(train_datas), function(i){
          if (all(train_datas[[i]]$combi == combi_x[k, , drop = FALSE])){
            i
          }else{
            NULL
          }
        }) %>% unlist()

        if (!parametric){
          #pour toutes les variables de Y
          res = lapply(seq_len(p), function(j){
            conditionnal_quantile(X = train_datas[[train_datas_idx_k]]$X_num,
                                  Y = train_datas[[train_datas_idx_k]]$Y[,j],
                                  x = x[mask, !is_fac],
                                  alpha = alpha_j[j],
                                  h_n = all_h[[j]][[train_datas_idx_k]])$y
          }) %>% as.data.frame() %>% mutate(id = mask)
        }else{
          res = lapply(seq_len(p), function(j){
            conditionnal_quantile(X = train_datas[[train_datas_idx_k]]$X_num,
                                  Y = train_datas[[train_datas_idx_k]]$Y[,j],
                                  x = x[mask, !is_fac],
                                  alpha = alpha_j[j],
                                  parametric = TRUE,
                                  beta = all_beta[[j]][[train_datas_idx_k]],
                                  deg = deg)
          }) %>% as.data.frame() %>% mutate(id = mask)
        }


        colnames(res) = c(col_names_Y, "id")
        res
      }) %>% bind_rows() %>% arrange(id) %>% select(-id)

      lapply(fn_to_opt, function(f){
        f(x, Y)
      }) %>% as.data.frame()

    }



  }else{ #si que des quanti

    if (!parametric){
      h_n = lapply(seq_along(Y), function(j) hopt(X = X, Y = Y[,j])$h )
      fn = function(x){

        Y = lapply(seq_len(p), function(j){
          conditionnal_quantile(X = X,
                                Y = Y[,j],
                                x = x,
                                alpha = alpha_j[j],
                                h_n = h_n[[j]])$y
        }) %>% as.data.frame()

        colnames(Y) = col_names_Y

        lapply(fn_to_opt, function(f){
          f(x, Y)
        }) %>% as.data.frame()

      }
    }else{
      beta = lapply(seq_along(Y), function(j) fit_beta(X = X, Y = Y[,j], alpha = alpha_j[j], deg = deg) )
      fn = function(x){

        Y = lapply(seq_len(p), function(j){
          conditionnal_quantile(X = X,
                                Y = Y[,j],
                                x = x,
                                alpha = alpha_j[j],
                                parametric = TRUE,
                                beta = beta[[j]])
        }) %>% as.data.frame()

        colnames(Y) = col_names_Y

        lapply(fn_to_opt, function(f){
          f(x, Y)
        }) %>% as.data.frame()

      }
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
  res = NSGA(X = X, fn, n_objective = length(fn_to_opt), g=g,
             sens = sens,
             N = N, seed = seed,
             B = B, verbose = TRUE)

  #mise en forme des resultats
  res$fn = function(X, Y){
    as.data.frame(sapply(fn_to_opt, function(f) f(X, Y)))
  }

  formals(res$fn)$Y = Y
  res$alpha = alpha

  res
}
