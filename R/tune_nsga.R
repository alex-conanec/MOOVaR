# library(tidyverse)
#
# #fonction de volume
#' @export
HI = function(Y_front, Y_MC, A_MC, N_MC, sens = rep("max", NCOL(Y_MC))){

  Y_front = sweep(Y_front, 2, sapply(sens, function(x) ifelse(x == "min", -1, 1)), "*")
  Y_MC = sweep(Y_MC, 2, sapply(sens, function(x) ifelse(x == "min", -1, 1)), "*")

  A = array(as.matrix(Y_front), dim = c(dim(Y_front), N_MC) )
  A = apply(A, c(3, 2, 1), function(x) x)

  B = array(Y_MC, dim = c(dim(Y_MC), dim(A)[3]))

  res=apply(A-B < 0, c(1, 3), any) %>%
    apply(1, all)

  (sum(!res)/N_MC)*A_MC

}

#' @export
evaluate = function(X, param_fixed, mc_cores = 1){

    N = NROW(X)
    X0 = lapply(1:N, function(i){
      X_quanti = matrix(runif(n = X$N[i]*param_fixed$d1),
                        ncol = param_fixed$d1)
      colnames(X_quanti) = paste0("X", 1:NCOL(X_quanti))

      idx = sample(1:NROW(param_fixed$quali_design), X$N[i], replace = TRUE)
      X_quali = param_fixed$quali_design[idx,]
      colnames(X_quali) = paste0("X", 1:NCOL(X_quali) + param_fixed$d1)
      cbind(X_quanti, X_quali)
    })

  res_all = parallel::mclapply(X = 1:N, mc.cores = 2, FUN = function(i){
      res_intra = parallel::mclapply(X = 1:param_fixed$B, mc.cores = mc_cores/2,
                                     FUN = function(j){

        seed = sample(1:10^6, 1)
        t = system.time({
          res_nsga = purrr::quietly(NSGA)(
            X = X0[[i]],
            fn = param_fixed$f,
            n_objective = param_fixed$p,
            sens = param_fixed$sens,
            type_var = param_fixed$type_var,

            #k tournois
            k_tournament = X[i,]$k_tournament,

            #crossing over
            crossing_over_method = X[i,]$crossing_over_method,
            n_echange = X[i,]$n_echange,
            n_bloc = X[i,]$n_bloc,

            #mutation
            mutation_method = X[i,]$mutation_method,
            freq_m=X[i,]$freq_m,

            #budget
            N = X[i,]$N,
            TT = X[i,]$TT,

            #verbose
            verbose = TRUE,
            front_tracking = TRUE,
            seed = seed
          )$result
        })[3]

        accuracy = HI(Y_front = res_nsga$Y,
                      Y_MC = param_fixed$Y_MC,
                      # A_MC = param_fixed$A_MC,
                      A_MC = 1,
                      N_MC = param_fixed$N_MC,
                      sens = param_fixed$sens)
        list(t = t, accuracy = accuracy,
             # res_nsga = res_nsga,
             seed = seed)
      })

      list(
        perf = data.frame(
          accuracy = mean(sapply(res_intra, function(x) x$accuracy)),
          t = mean(sapply(res_intra, function(x) x$t)),
          N_eval = X$TT[i] * X$N[i]
        ),
        seed = sapply(res_intra, function(x) x$seed)
      )
    })

  list(
    perf = bind_rows(lapply(res_all, function(x) x$perf)),
    memory = lapply(1:N, function(j) list(seed = res_all[[j]]$seed, X0 = X0[[j]]))
    # seed = lapply(res_all, function(x) x$seed)
  )

}

#' @export
tune_nsga = function(param,
                     param_fixed = list(p = 2,
                                        n_lev = 2:4,
                                        d1 = sample(9:12, 1),
                                        d2 = sample(2:4, 1),
                                        B = 10),
                     N = 30,
                     TT = 10,
                     base_N_MC = 30,
                     mc_cores = 1,
                     seed=1234){

  set.seed(seed)

  # param_fixed$lev = sample(param_fixed$n_lev, param_fixed$d2, replace = TRUE)
  # if (param_fixed$d2 > 2){
  #   param_fixed$lev = sample(param_fixed$n_lev, param_fixed$d2, replace = TRUE)
  # }else{
  #   param_fixed$lev = c(2, 5)
  # }


  param_fixed$quali_design = purrr::cross(
    lapply(param_fixed$n_lev, seq, from = 1)) %>%
    lapply(function(x){
      df = as.data.frame(x)
      colnames(df) = paste0("X", 1:NCOL(df))
      df
    }) %>% bind_rows() %>% mutate_if(is.integer, as.factor)

  beta = matrix(runif(n = param_fixed$p*(param_fixed$d1 + sum(param_fixed$n_lev)),
                      -1, 1), ncol = param_fixed$p)


  f = function(X, beta){
    is_fac = sapply(X, is.factor)
    if (any(is_fac)){
      X_dummy = lapply(which(is_fac), function(j){
        as.data.frame(model.matrix(~. + 0, data = X[,j, drop = FALSE]))
      }) %>% bind_cols()

      X = purrr::quietly(bind_cols)(X[!is_fac], X_dummy)$result
    }

    as.matrix(X) %*% beta

  }
  formals(f)$beta = beta

  param_fixed$f = f

  #set a X0 for cor
  nX0 = 100
  X_quanti = matrix(runif(n = nX0 * param_fixed$d1),
                    ncol = param_fixed$d1)
  colnames(X_quanti) = paste0("X", 1:NCOL(X_quanti))
  idx = sample(1:NROW(param_fixed$quali_design), nX0, replace = TRUE)
  X_quali = param_fixed$quali_design[idx,]
  X0 = cbind(X_quanti, X_quali)
  Y0 = param_fixed$f(X0)
  Y_MC_coord = data.frame(objectif = 1:NCOL(Y0),
                          min = floor(apply(Y0, 2, mean) - 3*apply(Y0, 2, sd)),
                          max = ceiling(apply(Y0, 2, mean) + 3*apply(Y0, 2, sd)))

  param_fixed$A_MC = apply(Y_MC_coord[,-1], 1, function(x) x[2] - x[1]) %>% prod()
  param_fixed$N_MC = base_N_MC^param_fixed$p
  param_fixed$Y_MC = sapply(1:param_fixed$p, function(j){
    runif(param_fixed$N_MC, min = unlist(Y_MC_coord[j, 2]),
          max = unlist(Y_MC_coord[j, 3]))
  })
  param_fixed$sens = rep("max", param_fixed$p)


  #nb de combinaison
  prod(sapply(param, length))

  X = lapply(param, function(x){
    sample(x, N, replace = TRUE)
  }) %>% bind_cols()

  #check type
  sapply(X, class)

  formals(evaluate)$param_fixed = param_fixed
  formals(evaluate)$mc_cores = mc_cores

  #### NSGA adapted
  all_front = NULL
  res_memory = NULL
  type_var = sapply(X, class)
  sens = c("max", "min", "min")
  k_tournament = 4
  n_echange = 3
  freq_m = 0.3
  X0 = X

  distri_Xi = lapply(param, function(x){
    if (class(x) == "factor"){
      list(levels = levels(x))
    }else{
      list(min = min(x), max = max(x))
    }
  })

  param_mut = lapply(1:N, function(i){
    list(alpha = matrix(0, nrow = sum(type_var != "factor"),
                        ncol = sum(type_var != "factor")),
         sigma = sapply(X[,which(type_var!="factor")], sd)*1)
  })


  # evaluate the function at X
  res_eval = evaluate(X)
  res_memory = res_eval$memory

  Y <- res_eval$perf %>% as.data.frame() %>% mutate(id = 1:NROW(.))
  Y$rank = dominance_ranking(Y[,-NCOL(Y)], sens)
  Y = crowding_distance(Y)

  time_deb = Sys.time()
  for (t in 1:TT){
    parents_idx <- tournament_selection(Y, k = k_tournament,
                                        N = ifelse(N%%2==0, N, N+1))

    res_cross = crossing_over_uniforme(
      parents = X[parents_idx,],
      n_echange = n_echange,
      param_mut = lapply(parents_idx, function(i) param_mut[[i]]))

    X_C = res_cross$X_C
    param_mut_C = res_cross$param_mut

    res_mut = mutation_mixte(X = X_C, freq_m = freq_m,
                             param_mut = param_mut_C,
                             distri_Xi = distri_Xi,
                             type_var = type_var)
    X_C = res_mut$X
    param_mut_C = res_mut$param_mut

    rownames(X_C)=1:NROW(X_C)

    X_C = X_C[seq_len(N),]
    param_mut = c(param_mut, param_mut_C[1:N])

    #evaluate the children
    res_eval = evaluate(X_C)
    res_memory = c(res_memory, res_eval$memory)

    Y_C = res_eval$perf %>% as.data.frame() %>% mutate(id = 1:NROW(.) + N)
    Y = bind_rows(Y %>% select(-rank, -crowding_distance) %>%
                    mutate(id = 1:N),
                  Y_C)

    Y$rank = dominance_ranking(Y[,-NCOL(Y)], sens)
    Y = crowding_distance(Y)

    # Select the best
    Y = Y %>% arrange(rank, desc(crowding_distance)) %>%
      filter(crowding_distance != 0) %>%
      head(N)

    # 5) Gather the Pt and Qt and proceed to the next iteration
    X = bind_rows(X, X_C)[Y$id,]
    param_mut = param_mut[Y$id]
    res_memory = res_memory[Y$id]
    Y = Y %>% mutate(id = 1:N)

    all_front[[t]] = data.frame(t=t, Y[, 1:3], rank = Y$rank)

    cat("Population ", t, '\n')
    print(round((Sys.time() - time_deb), 1))
    cat("\n")

  }

  X = as.data.frame(X)

  #if crowding distance parameter thing
  no_domi = Y$rank == 1
  Y <- Y[no_domi, 1:3]
  X = X[no_domi,]
  res_memory = res_memory[which(no_domi)]

  #tri
  Y = Y %>% mutate(id = as.factor(1:NROW(Y))) %>% arrange_all()
  X = X[as.numeric(Y$id),]
  res_memory = res_memory[as.numeric(Y$id)]
  Y = Y %>% select(-id)

  rownames(X) <- seq_len(NROW(X))
  res = list(X = X, Y = Y, time = Sys.time() - time_deb, fn = evaluate,
             mc_cores = mc_cores, X0 = X0, all_front = all_front,
             sens = sens, res_memory = res_memory)
  class(res) = "tune_nsga"
  res

}


#' @export
plot.tune_nsga = function(res, choice = "PF", id = NULL, mc_cores = 1,
                          as_list_plot = FALSE){

  if (choice == "PF"){
    tmp = res
    class(tmp) = "nsga"
    plot(tmp, choice = "PF", as_list_plot = as_list_plot)
  }else if (choice == "evol"){
    tmp = res
    res=tmp
    class(tmp) = "nsga"
    plot(tmp, choice = "evol")
  }else if (choice == "qlt"){
    tmp = res
    class(tmp) = "nsga"
    plot(tmp, choice = "qlt", mc_cores = mc_cores)
  }else if (choice == "comp"){

    df = lapply(1:length(res$res_nsga), function(i){

      parallel::mclapply(res$res_nsga[[i]], function(x){
        df = x$all_front %>% bind_rows() %>% mutate(t = as.factor(t)) %>%
          filter(rank == 1) %>% select(-rank)
        data.frame(
          id = id[i],
          t = 1:max(as.numeric(df$t)),
          Aire = sapply(1:max(as.numeric(df$t)), function(tt){
            HI(Y_front = df %>% filter(t == tt) %>% select(-t),
               Y_MC = formals(res$fn)$param_fixed$Y_MC,
               # A_MC = formals(res$fn)$param_fixed$A_MC,
               A_MC = 1,
               N_MC = formals(res$fn)$param_fixed$N_MC,
               sens = formals(res$fn)$param_fixed$sens)

          })
        )
      }, mc.cores = mc_cores) %>% bind_rows()
    }) %>% bind_rows()


    df %>% group_by(id, t) %>% summarise(moy = mean(Aire),
                                         q_.1 = quantile(Aire, 0.1),
                                         q_.9 = quantile(Aire, 0.9)) %>%
      mutate(id = as.factor(id)) %>%

      ggplot(aes(x = t)) +
      geom_line(aes(y = moy, colour = id))+
      geom_ribbon(aes(ymin= q_.1, ymax= q_.9, fill = id), alpha = 0.3) +
      ylab("HI")


  }else if (choice == "comp_front"){
    print("to dev")
  }

}


#' @export
res_nsga_run = function(res, id, mc_cores=1){

  lapply(id, function(i){
    seeds = res$res_memory[[i]]$seed
    X = res$X
    X0 = res$res_memory[[i]]$X0
    parallel::mclapply(X = 1:length(seeds), mc.cores = mc_cores, FUN = function(j){
      purrr::quietly(NSGA)(
        X = X0,
        fn = formals(res$fn)$param_fixed$f,
        n_objective = formals(res$fn)$param_fixed$p,
        sens = formals(res$fn)$param_fixed$sens,
        type_var = formals(res$fn)$param_fixed$type_var,

        #k tournois
        k_tournament = X[i,]$k_tournament,

        #crossing over
        crossing_over_method = X[i,]$crossing_over_method,
        n_echange = X[i,]$n_echange,
        n_bloc = X[i,]$n_bloc,

        #mutation
        mutation_method = X[i,]$mutation_method,
        freq_m=X[i,]$freq_m,

        #budget
        N = X[i,]$N,
        TT = X[i,]$TT,

        #verbose
        verbose = TRUE,
        front_tracking = TRUE,
        seed = seeds[j]
      )$result
    })
  })
}
