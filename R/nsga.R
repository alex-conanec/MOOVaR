#' NSGA II algorithm
#'
#' The non sorting genetic algorithm is finding the Pareto front
#'
#' @param X a matrix/data.frame of an initial population
#'
#' @param ff a list of list containing the function or predicting model and
#' if the objectif function should be maximise "max" or minimize "min"
#'
#' @param N the size of the population P before the generation of individus.
#' By default N=NROW(X)/2
#'
#' @param crossing_over_size number of variable to switch when doing the
#' crossing over
#'
#' @param freq_mutation a vector of numerics value varying from 0 to 1 which
#' represent the frequency of mutation of each variable
#'
#' @return a vecteur of rank associated with the individus of X
#'
#' @examples
#' library(tidyverse)
#' set.seed(1234)
#' n <- 300
#' q = 3
#' mini <- rep(x = -2, times = q)
#' maxi <- rep(x = 2, times = q)
#'
#' X <- lapply(seq_len(q), function(k){
#'   res <- data.frame(X=runif(n = n, min = mini[k], max = maxi[k]))
#'   names(res) <- paste0("X", k)
#'   res
#' }) %>% bind_cols()
#'
#' fn <- list(
#'   Y1 = function(X) unlist(X[1] + 2*X[2] + X[3]),
#'   Y2 = function(X) unlist(-X[1] - X[2] - 3*X[3])
#' )
#'
#' res = NSGA(X = X, fn, N=50)
#' plot(res$Y)
#'
#' res = NSGA(X = X, fn, N=50, sens = rep("max", length(fn)))
#' plot(res$Y)
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @export

NSGA <- function(X, fn, n_objective, sens = rep("min", n_objective),
                 N = NROW(X), g = NULL,
                 k_tournament, n_echange=1, n_bloc = 2,
                 crossing_over_method = "uniforme",
                 mutation_method = "simple",
                 freq_m=0.2,
                 type_var = NULL,
                 distri_Xi = NULL,
                 seed = NULL, TT=50, verbose = TRUE, front_tracking = TRUE,
                 updateProgress = NULL,
                 path_tracking = NULL){

    #tracking
    if (is.function(updateProgress)) {
        updateProgress(detail = "NSGA init")
    }
    tracking_msg(path_tracking, msg = "NSGA init")

    all_front = NULL

    if (!is.null(seed)) set.seed(seed)
    X = as.data.frame(X)
    X0 = X

    if (is.null(distri_Xi)){
        distri_Xi <- lapply(X, function(x){
            if (class(x) == "numeric"){
                list(min = floor(min(x)), max = ceiling(max(x)),
                     mean = mean(x), sd = sd(x))
            }else if (class(x) == "factor"){
                list(levels = levels(x))
            }
        })
    }


    if (is.null(type_var)) type_var = sapply(X, class)

    #check constraint for the initial pop
    K = length(g)
    if (K > 0){
        feasible = sapply(g, function(gg){
            gg(X)
        })
        # if (K > 1){
        feasible = apply(feasible, 1, all)
        # }
        X = X[feasible,]

        threshold = 5
        if (NROW(X) < threshold){
            stop(paste0("must provide at least", threshold, "feasible solution at the begenning"))
        }
    }

    if (mutation_method == "mixte"){
        param_mut = lapply(1:N, function(i){
            list(alpha = matrix(0, nrow = sum(type_var != "factor"),
                                ncol = sum(type_var != "factor")),
                 sigma = sapply(X[,which(type_var!="factor")], sd)*1)
        })
    }else{
        param_mut = lapply(1:N, function(i) NA)
    }

    # X = X[1:min(2*N, NROW(X)), ]
    X = head(X, N)
    time_deb <- Sys.time()

    # 1) evaluate the function at X
    Y <- fn(X) %>% as.data.frame() %>% mutate(id = 1:NROW(.))
    Y$rank = dominance_ranking(Y[,-NCOL(Y)], sens)
    Y = crowding_distance(Y)


    for (t in 1:TT){

        if (is.function(updateProgress)) {
            updateProgress(detail = paste0("NSGA, pop: ", t))
        }
        tracking_msg(path_tracking, msg = paste0("NSGA, pop: ", t))


        # while not generate new feasible possible solution
        res_X_C = NULL
        res_param_mut_C = NULL
        counter = 1

        is_it_stopped(path_tracking)

        while (NROW(res_X_C) < N){

            # #tracking stop
            # if (!is.null(path_tracking)){
            #     tracking = readRDS(path_tracking)
            #     if (tracking[[length(tracking)]] == "stop"){
            #         stop()
            #     }
            # }

            # 3) Pick a couples to be reproduce after a tournement selection
            parents_idx <- tournament_selection(Y, k = k_tournament,
                                                N = ifelse(N%%2==0, N, N+1))

            # #pas cool, il selectionne que les plus petit index... surement pas les nvx
            # lapply(1:10, function(s){
            #     tournament_selection(Y, k = k_tournament,
            #                          N = ifelse(N%%2==0, N, N+1))
            # }) %>% unlist() %>% hist()


            # 4) Generate new individu (Qt) with genetic operator apply on each
            if (crossing_over_method == "uniforme"){
                res_cross = crossing_over_uniforme(
                    parents = X[parents_idx,],
                    n_echange = n_echange,
                    param_mut = param_mut[parents_idx])
                X_C = res_cross$X_C
                param_mut_C = res_cross$param_mut
            }else if (crossing_over_method == "bloc"){
                res_cross = crossing_over_bloc(
                    parents = X[parents_idx,],
                    n_bloc = n_bloc,
                    param_mut = param_mut[parents_idx])
                X_C = res_cross$X_C
                param_mut_C = res_cross$param_mut
            }

            if (mutation_method == "simple"){
                X_C = mutation_simple(X_C, freq_m = freq_m, distri_Xi = distri_Xi,
                                      type_var = type_var)

            }else if (mutation_method == "mixte"){
                res_mut = mutation_mixte(X = X_C, freq_m = freq_m,
                                         param_mut = param_mut_C,
                                         distri_Xi = distri_Xi,
                                         type_var = type_var)
                X_C = res_mut$X
                param_mut_C = res_mut$param_mut
            }

            rownames(X_C)=1:NROW(X_C)

            #check that the new individu respect the g constraint
            if (K > 0){
                cat("check constraint:", counter, "\n")
                feasible = sapply(g, function(gg){
                    gg(x=X_C)
                })

                # if (K > 1){
                    feasible = apply(feasible, 1, all)
                # }
                res_X_C = rbind(res_X_C, X_C[feasible,])
                res_param_mut_C = c(res_param_mut_C, param_mut_C[which(feasible)])
                counter = counter + 1
                print(NROW(res_X_C))
            }else{
                res_X_C = X_C
                res_param_mut_C = param_mut_C
            }

        }

        res_X_C = res_X_C[seq_len(N),]
        param_mut = c(param_mut, res_param_mut_C[1:N])

        #evaluate the children
        Y_C = fn(res_X_C) %>% as.data.frame() %>% mutate(id = 1:NROW(.) + NROW(Y))
        Y = bind_rows(Y %>% select(-rank, -crowding_distance) %>%
                          mutate(id = 1:NROW(Y)),
                      Y_C)

        Y$rank = dominance_ranking(Y[,-NCOL(Y)], sens)
        Y = crowding_distance(Y)

        # Select the best
        Y = Y %>% arrange(rank, desc(crowding_distance)) %>%
            filter(crowding_distance != 0) %>%
            head(N)

        # 5) Gather the Pt and Qt and proceed to the next iteration
        X = bind_rows(X, res_X_C)[Y$id,]
        param_mut = param_mut[Y$id]
        Y = Y %>% mutate(id = 1:NROW(Y))

        if (front_tracking){
            all_front[[t]] = data.frame(t=t, Y[, 1:n_objective], rank = Y$rank)
        }

        if (verbose) {
            cat("Population ", t, '\n')
            print(round((Sys.time() - time_deb), 1))
            cat("\n")
        }
    }


    #if crowding distance parameter thing
    no_domi = Y$rank == 1
    Y <- Y[no_domi, seq_len(n_objective)]
    X = X[no_domi,]

    #tri
    Y = Y %>% mutate(id = as.factor(1:NROW(Y))) %>% arrange_all()
    X = X[Y$id,]
    Y = Y %>% select(-id)

    # Y <- Y[, seq_along(fn)] #a remettre si crowding distance thing est supprime finalement.
    # colnames(Y) <- names(fn)
    rownames(X) <- seq_len(NROW(X))
    res = list(X = X, Y = Y, time = Sys.time() - time_deb, fn = fn, g = g,
               X0 = X0, all_front = all_front, sens = sens)
    class(res) = "nsga"
    res
}


#' plot nsga
#'
#' @export

plot.nsga <- function(res, choice = "PF", alpha_ellipse = NULL, mc_cores = 1,
                      as_list_plot = FALSE){

    library(magrittr)
    if (choice == "PF"){
        # colnames(res$Y0) = colnames(res$Y)
        is_fac = !sapply(res$X, is.numeric)
        for (k in which(is_fac)) res$X[,k] = as.factor(res$X[,k])

        df = dplyr::bind_cols(res$X, res$Y) %>% dplyr::mutate(id = 1:NROW(.))
        dd = plotly::highlight_key(df, ~id)

        combi_y = combn(colnames(res$Y), 2)

        lapply(seq_len(NCOL(combi_y)), function(i){

            p = ggplot2::ggplot(dd) +
                ggplot2::geom_point(ggplot2::aes_string(x = combi_y[1, i],
                                                        y = combi_y[2, i]),
                                    colour = "blue", size = 2) +
                ggplot2::xlab(combi_y[1, i]) + ggplot2::ylab(combi_y[2, i]) +
                ggplot2::geom_point(data = res$Y0,
                                    ggplot2::aes_string(x = combi_y[1, i],
                                                        y = combi_y[2, i]),
                                    alpha = 0.3)

            if (NCOL(res$Y) < 3){
                p + ggplot2::geom_line(data = df,
                                       ggplot2::aes_string(x = combi_y[1, i],
                                                           y = combi_y[2, i]),
                                       colour = "blue")
            }else{
                p
            }

        }) -> p_y

        lapply(colnames(res$X), function(X_i){
            if (is.numeric(res$X[, X_i]) | is.integer(res$X[, X_i])){
                p = ggplot2::ggplot(dd) +
                    ggplot2::geom_point(ggplot2::aes_string(x = 1, y = X_i)) +
                    ggplot2::xlab(X_i) +
                    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                   axis.ticks.x = ggplot2::element_blank())
            }else{
                p = ggplot2::ggplot(dd) + ggplot2::geom_text(
                    ggplot2::aes(x = 1, y = as.numeric(!! rlang::sym(X_i)),
                                 label = !! rlang::sym(X_i))) +
                    ggplot2::ylim(0.5, nlevels(res$X[, X_i])+0.5) +
                    ggplot2::xlab(X_i) +
                    ggplot2::theme(axis.text = ggplot2::element_blank(),
                                   axis.ticks = ggplot2::element_blank(),
                                   axis.title.y = ggplot2::element_blank())
            }

            # if (NCOL(res$X) > 5){
            #     p + ggplot2::theme(axis.title.x = ggplot2::element_text(angle = 90))
            # }else{
            p
            # }
        }) -> p_x


        if (!as_list_plot){
            n_y = length(p_y)
            n_x = length(p_x)

            s_y = plotly::subplot(p_y, titleY = TRUE, titleX = TRUE) #, nrows = ceiling(n_y/3))
            s_x = plotly::subplot(p_x, titleX = TRUE) #, nrows = ceiling(n_x/10))

            if (is.null(res$tau)){
                pareto_title = "Pareto front"
            }else{
                pareto_title = paste0("Pareto front (alpha=", res$tau, ")")
            }


            plotly::ggplotly(plotly::subplot(s_y, s_x, nrows = 2, margin = 0.05,
                                             titleY = TRUE, titleX = TRUE)) %>%
                plotly::highlight(on = "plotly_selected", off= "plotly_deselect",
                                  color = "red") %>%
                plotly::layout(annotations = list(
                    list(x = 0 , y = 1.05, text = pareto_title,
                         showarrow = F, xref='paper', yref='paper'),
                    list(x = 0 , y = 0.48, text = "Decision solution",
                         showarrow = F, xref='paper', yref='paper'))
                )
        }else{
            list(p_x = p_x, p_y = p_y, combi_y = combi_y)
        }

    }else if (choice == "qlt"){


        df = res$all_front %>% bind_rows() %>% select(-rank, -t)
        Y_MC_coord = data.frame(min = floor(apply(df, 2, min)),
                                max = ceiling(apply(df, 2, max)) )

        p = NROW(Y_MC_coord)

        N_MC = 30^p
        # A_MC = apply(Y_MC_coord, 1, function(x) x[2] - x[1]) %>% prod()
        A_MC = 1
        Y_MC = sapply(1:p, function(j){
            runif(N_MC, min = unlist(Y_MC_coord[j, 1]),
                  max = unlist(Y_MC_coord[j, 2]))
        })

        df = res$all_front %>% bind_rows() %>% mutate(t = as.factor(t)) %>%
            filter(rank == 1)
        df = df[order(df[,2]),]

        qlt <- data.frame(
                t = 1:max(as.numeric(df$t)),
                aire = parallel::mclapply(1:max(as.numeric(df$t)), function(tt){
                    HI(Y_front = df %>% filter(t == tt) %>% select(-t, -rank),
                       Y_MC, A_MC, N_MC, sens = res$sens)
                }, mc.cores = mc_cores) %>% unlist()
            )

        ggplot(qlt, aes(x=as.numeric(t), y=aire)) +
            geom_line() +
            xlab('t')


    }else if (choice == "evol"){
        df = res$all_front %>% bind_rows() %>% mutate(t = as.factor(t)) %>%
            filter(rank == 1)

        df = df[order(df[,2]),]

        ggplot(df, aes_string(x = colnames(df)[2], y = colnames(df)[3])) +
            geom_point(aes(colour = t)) +
            geom_line(aes(group = t, colour = t))
    }

}
