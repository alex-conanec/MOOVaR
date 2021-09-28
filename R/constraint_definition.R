#' def_cstr_X_space
#'
#' Find the space where points where observed. Send back a function to check if
#' a/several point(s) belong to the space where the points where observed
#'
#'  @param X matrix of dimension n x d observed points
#'  @param alpha probability to deal with the trade-off false positive false
#'  negatif. Alpha closed to 1 mean low false positive and vice versa
#'  @param path_tracking path where to write the step of the running function.
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



def_cstr_X_space <- function(X, alpha = 0.15,
                             allowed_dependence = matrix(TRUE, nrow = NCOL(X), ncol = 1),
                             path_tracking = NULL){ #/!\ pb qd n d'une combi < a p+1, non ?

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


    #reduction dimensionnelle
    # eig = eigen(cor(X))
    # P = eig$vectors


    #isolation forest
    theta = 0.5
    train_idx = sample(x = 1:NROW(X), size = theta*NROW(X))
    X_train1 = X[train_idx,]
    X_train2 = X[-train_idx,]

    iforest <- isolationForest$new(sample_size = NROW(X_train1),
                                   num_trees = 1000)
    iforest$fit(dataset = X_train1)
    pred_train = iforest$predict(X_train2)
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
    dep = apply(allowed_dependence, 1, any)
    row_to_del = NULL
    col_to_del = NULL
    for (i in 1:NCOL(res)){
      is_na = is.na(res[,i])
      if (dep[i] & any(is_na)){
        row_to_del = unique(c(row_to_del, which(is_na)))
      }else if (!dep[i]){
        col_to_del = c(col_to_del, i)
      }
    }

    if (!is.null(col_to_del)){
      res = res[,-col_to_del]
    }
    if (!is.null(row_to_del)){
      res = res[-row_to_del,]
    }
    res = lapply(res, as.factor) %>% data.frame()

  }else if (any(is_fac)){
    combinaison = data.frame(X) %>% dplyr::group_by_if(is.factor) %>%
      dplyr::summarise(n = n()) %>% filter(n > (sum(!is_fac) + 1)) %>% select(-n) #comment il gere ceux qui ont moins de (sum(!is_fac) + 1) ind ?

    K = NROW(combinaison)
    lapply(seq_len(K), function(k){

      is_it_stopped(path_tracking)

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

  g = function(x, res, d, path_tracking){
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
          is_it_stopped(path_tracking)
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

      check_combi = function(X_test, X_ref){
        XX = array(as.matrix(X_ref), dim = c(dim(X_ref), NROW(X_test)))
        YY = array(as.matrix(X_test), dim = c(dim(X_test), NROW(X_ref))) %>%
          apply(c(3, 2, 1), function(x) x)
        apply(XX - YY == 0, c(1,3), all) %>% apply(2, any)
      }

      lev_fac = lapply(res, function(yy){
        data.frame(lev = levels(yy), num = 1:nlevels(yy))
      })

      xx = lapply(1:NCOL(res), function(i=1){
        df = data.frame(lev = x[,i,T]) %>% left_join(lev_fac[[i]], by = "lev") %>%
          select(num)
        is_na = is.na(res)
        if (any(is_na)) df[is_na] = 0.5
        df
      }) %>% as.data.frame()
      colnames(xx) = colnames(res)

      yy = lapply(res, as.numeric) %>% as.data.frame()

      check_combi(X_test = xx, X_ref = yy)

    }


  }

  formals(g)$res = res
  formals(g)$d = NCOL(X)
  formals(g)$path_tracking = path_tracking

  res = list(
    g = g,
    X = X
  )

  class(res) = "combi_space"
  res

}







#' @import ggplot2
#' @importFrom plotly ggplotly highlight subplot highlight_key
#'
#' @export
plot.combi_space <- function(res, as_list_plot = FALSE){


  d = NCOL(res$X)
  n = NROW(res$X)
  is_num = sapply(res$X, is.numeric)

  mini = floor(apply(res$X[, is_num], 2, min))
  maxi = ceiling(apply(res$X[, is_num], 2, max))

  combinaison = unique(res$X[,!is_num])
  n_MC = 10 ^ sum(is_num)
  X_quanti_MC = sapply(seq_len(d)[is_num], function(i){
    runif(n_MC, mini[i], maxi[i])
  }) %>% as.data.frame()

  if (any(is_num)){
    df = lapply(1:NROW(combinaison), function(i){
      cbind(combinaison[i,], X_quanti_MC)
    }) %>% bind_rows()
    df = df[res$g(df),]
  }else{
    df = combinaison
  }

  df$id = 1:NROW(df)

  dd = highlight_key(df, ~id)
  lapply(colnames(df)[-NCOL(df)], function(X_i){
    if (is.numeric(df[, X_i]) | is.integer(df[, X_i])){
      p = ggplot(dd) +
        geom_point(aes_string(x = 1, y = X_i)) +
        xlab(X_i) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
    }else{
      p = ggplot(dd) + geom_text(
        aes(x = 1, y = as.numeric(!! rlang::sym(X_i)),
            label = !! rlang::sym(X_i))) +
        ylim(0.5, nlevels(df[, X_i]) + 0.5) +
        xlab(X_i) +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title.y = element_blank())
    }
    p + theme_bw()
  }) -> p_x

  if (!as_list_plot){
    ggplotly(subplot(p_x, titleX = TRUE)) %>%
      highlight(on = "plotly_selected", off= "plotly_deselect",
                color = "red")
  }else{
    p_x
  }

}





