#' def_cstr_X_space
#'
#' Find the space where points where observed. Send back a function to check if
#' a/several point(s) belong to the space where the points where observed
#'
#'  @param X matrix of dimension n x d observed points
#'  @param m a vector of size d to indicate how many reference points create for each dimension.
#'
#' @examples
#' library(tidyverse)
#' library(optisure)
#'
#' #parameters constraint
#' alpha = 0.05
#'
#' #parameters data
#' n = 1000
#' d1 = 2
#' mini = rep(0, d1)
#' maxi = rep(10, d1)
#' b = matrix(c(0, 0, 5, 5), byrow = T, ncol = d1)
#' B = matrix(c(4, 4, 10, 10), byrow = T, ncol = d1)
#' p = NROW(B) #nombre de "foyer"
#'
#'
#' #generate data
#' X_U = sapply(seq_len(d1), function(i){
#'   sapply(seq_len(p), function(k){
#'     runif(n/p, b[k,i], B[k,i])
#'   })
#' }) %>% as.data.frame()
#'
#' #Monte Carlo simulation of 5000 points to check the constraint
#' n_MC = 5000
#' X_MC = sapply(seq_len(d1), function(i){
#'   runif(n_MC, mini[i], maxi[i])
#' })
#'
#' #visualisation
#' res = def_cstr_X_space(X=X_U, alpha = alpha)
#' feasible_X_MC = res$g(x = X_MC)
#'
#' cols <- c("FALSE" = "red", "training" = "blue", "TRUE" = "green")
#' ggplot(data.frame(X_MC, feasible = feasible_X_MC)) +
#'   geom_point(aes(x = X1, y = X2, colour = feasible)) +
#'   geom_point(data = X_U, aes(x = V1, y = V2, colour = "training"),
#'              shape = 17, size = 1) +
#'   scale_colour_manual(values = cols)
#'
#'
#' # verification avec les aires de MC
#' # aire theorique de 95% de la surface
#' A_true_uni = (1-alpha) * sum(apply(B-b, 1, prod))
#' A_true_uni
#'
#' # Aire estime par Monte Carlo
#' A_MC = prod(maxi - mini)
#' A_estime = A_MC*sum(feasible_X_MC)/n_MC
#' A_estime
#' (A_estime - A_true_uni)/A_true_uni #ratio
#'
#' @import solitude
#' @export



def_cstr_X_space <- function(X, alpha = 0.15){ #/!\ pb qd n d'une combi < a p+1, non ?

  X = as.data.frame(X)
  is_fac = !sapply(data.frame(X), is.numeric)

  for (k in which(is_fac)) X[,k] = as.factor(X[,k, drop=T])

  d1 = sum(!is_fac)
  d2 = sum(is_fac)
  n = NROW(X)

  novelty_detection = function(X, alpha){

    # cst_col = sapply(X, function(x) length(unique(x)) == 1)
    # if (any(cst_col)){
    #   values_cst = unique(X[,cst_col])
    #   X = X[,!cst_col]
    # }else{
    #   values_cst = NULL
    # }

    #isolation forest
    iforest <- isolationForest$new(sample_size = NROW(X),
                                   num_trees = 1000)
    iforest$fit(dataset = X)
    pred_train = iforest$predict(X)
    threshold = quantile(pred_train$anomaly_score, 1 - alpha)


    list(
      # cst_col = cst_col,
      # values_cst = values_cst,
      iforest = iforest,
      threshold = threshold
    )
  }

  if (all(is_fac)){
    res = unique(X)
  }else if (any(is_fac)){
    combinaison = data.frame(X) %>% dplyr::group_by_if(is.factor) %>%
      dplyr::summarise(n = n()) %>% filter(n > (sum(!is_fac) + 1)) %>% select(-n) #comment il gere ceux qui ont moins de (sum(!is_fac) + 1) ind ?

    K = NROW(combinaison)
    lapply(seq_len(K), function(k){
      mask = sapply(seq_len(n), function(i) all(X[i, is_fac, drop = FALSE] == combinaison[k,]))
      X_num = X[mask, !is_fac]
      res = novelty_detection(X=X_num, alpha = alpha) #/K

      res$combi = combinaison[k,]
      res$mask = mask
      res
    }) -> res

  }else{
    res = novelty_detection(X, alpha)
  }

  g = function(x, res, d){

    if (NCOL(x) < d){
      x = matrix(x, ncol = d)
    }
    x = as.data.frame(x)
    n = NROW(x)
    id = seq_len(n)

    if (class(res) == "list"){
      if (class(res[[1]]) == "list"){ #melange quanti / quali
        is_fac = !sapply(x, is.numeric)
        for (k in which(is_fac)) x[,k] = as.factor(x[,k, drop=TRUE])
        combinaison = unique(x[,is_fac, drop = FALSE])

        lapply(seq_len(NROW(combinaison)), function(k){
          mask = sapply(seq_len(n), function(i) all(x[i, is_fac, drop = FALSE] == combinaison[k, ,drop = FALSE]))

          res_idx_k = sapply(seq_along(res), function(i){
            if (all(res[[i]]$combi == combinaison[k, , drop = FALSE])){
              i
            }else{
              NULL
            }
          }) %>% unlist()

          if (!is.null(res_idx_k)){ #si le n ind par comb etait suffisant
            res_k = res[[res_idx_k]]

            # #check cst
            # n_k = length(which(mask))
            # res_cst = rep(TRUE, n_k)
            # if (any(res_k$cst_col)){
            #   res_cst = x[mask,!is_fac][,res_k$cst_col, drop=FALSE] == #a checker quand il ya plus d'une cst
            #     matrix(res_k$values_cst, nrow = n_k,
            #            ncol = length(res_k$values_cst), byrow = T)
            #   x = x[mask,!is_fac][,!res_k$cst_col]
            # }else{
              x = x[mask, !is_fac]
            # }

            feasible = res_k$iforest$predict(x)$anomaly_score <= res_k$threshold

            data.frame(
              id = id[mask],
              feasible = feasible)
              # feasible = apply(cbind(feasible, res_cst), 1, all))
          }else{ #sinon False
            data.frame(id = id[mask], feasible = FALSE)
          }


        }) %>% bind_rows() %>% arrange(id) %>% pull(feasible)

      }else{ #que des quanti

        res$iforest$predict(x)$anomaly_score <= res$threshold

      }

    }else{ #que des quali
      combinaison = unique(x)
      id = seq_len(n)
      lapply(seq_len(NROW(combinaison)), function(k){

        mask = sapply(seq_len(n), function(i) all(x[i, , drop = FALSE] == combinaison[k, ,drop = FALSE]))
        res_idx_k = sapply(seq_len(NROW(res)), function(i){
          if (all(res[i, , drop = FALSE] == combinaison[k, , drop = FALSE])){
            i
          }else{
            NULL
          }
        }) %>% unlist()

        data.frame(id = id[mask], feasible = !is.null(res_idx_k))

      }) %>% bind_rows() %>% arrange(id) %>% pull(feasible)

    }


  }

  formals(g)$res = res
  formals(g)$d = NCOL(X)

  res = list(
    g = g,
    X = X
  )

  class(res) = "combi_space"
  res

}







#' @export
plot.combi_space <- function(res, filtre = NULL){


  d = NCOL(res$X)
  n = NROW(res$X)
  is_num = sapply(res$X, is.numeric)

  mini = floor(apply(res$X[, is_num], 2, min))
  maxi = ceiling(apply(res$X[, is_num], 2, max))

  n_MC = 5000
  X_MC = sapply(seq_len(d), function(i){
    runif(n_MC, mini[i], maxi[i])
  }) %>% as.data.frame()

  combinaison = unique(res$X[,!is_num])
  X_MC = cbind(X_MC, combinaison[sample(seq_len(NROW(combinaison)),
                                        size = n_MC, replace = TRUE), ])



  if (!is.null(colnames(res$X))){
    colnames(X_MC) = colnames(res$X)
  }else{
    colnames(X_MC) = paste0("X_", seq_along(X_MC))
  }


  if (!is.null(filtre)){

    df_res = X_MC
    for (x_name in names(filtre)){

      if (is.numeric(res$X[,x_name])){
        df_res = df_res %>% filter(!!rlang::sym(x_name) > filtre[[x_name]][1] &
                                         !!rlang::sym(x_name) < filtre[[x_name]][2])
      }else{
        df_res = df_res %>% filter(!!rlang::sym(x_name) %in% filtre[[x_name]])
      }


    }

  }


  df_res_num = df_res[,is_num] %>% mutate(id = seq_len(NROW(.)),
                             feasible = res$g(df_res)) %>%
        gather(key = X_j, value = coord, -id, - feasible) %>% select(-id) %>%
    mutate(x = as.numeric(as.factor(X_j)))

  df_res_fac = df_res[,!is_num] %>% mutate(id = seq_len(NROW(.)),
                                          feasible = res$g(df_res)) %>%
    gather(key = X_j, value = coord, -id, - feasible) %>% select(-id)

  ggplot(df_res_num, aes(x = x, y = coord)) +
    facet_wrap(~X_j, "free") +
    geom_point(aes(colour = feasible ))

  ggplot(df_res_fac, aes(x = X_j, y = coord)) +
    geom_point(aes(colour = feasible )) +
    facet_wrap(~X_j) #, "free")
}





