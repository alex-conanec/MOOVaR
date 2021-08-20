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



def_cstr_X_space <- function(X, alpha = 0.15, path_tracking = NULL){ #/!\ pb qd n d'une combi < a p+1, non ?

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
      combinaison = unique(x)
      id = seq_len(n)
      lapply(seq_len(NROW(combinaison)), function(k){

        is_it_stopped(path_tracking)

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
  formals(g)$path_tracking = path_tracking

  res = list(
    g = g,
    X = X
  )

  class(res) = "combi_space"
  res

}







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

  dd = plotly::highlight_key(df, ~id)
  lapply(colnames(df)[-NCOL(df)], function(X_i){
    if (is.numeric(df[, X_i]) | is.integer(df[, X_i])){
      p = ggplot2::ggplot(dd) +
        ggplot2::geom_point(ggplot2::aes_string(x = 1, y = X_i)) +
        ggplot2::xlab(X_i) +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank())
    }else{
      p = ggplot2::ggplot(dd) + ggplot2::geom_text(
        ggplot2::aes(x = 1, y = as.numeric(!! rlang::sym(X_i)),
                     label = !! rlang::sym(X_i))) +
        ggplot2::ylim(0.5, nlevels(df[, X_i]) + 0.5) +
        ggplot2::xlab(X_i) +
        ggplot2::theme(axis.text = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank())
    }
    p
  }) -> p_x

  if (!as_list_plot){
    plotly::ggplotly(plotly::subplot(p_x, titleX = TRUE)) %>%
      plotly::highlight(on = "plotly_selected", off= "plotly_deselect",
                        color = "red")
  }else{
    p_x
  }

}





