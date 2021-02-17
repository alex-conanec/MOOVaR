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
                 N = round(NROW(X) / 2, 0), g = NULL,
                 crossing_over_size = round(NCOL(X)/2),
                 freq_mutation=rep(1/NCOL(X), NCOL(X)),
                 seed = NULL, B=50, verbose = TRUE){

    if (!is.null(seed)) set.seed(seed)
    X = as.data.frame(X)
    X0 = X
    distri_Xi <- lapply(X, function(x){
        if (class(x) == "numeric"){
            list(min = floor(min(x)), max = round(max(x)),
                 mean = mean(x), sd = sd(x))
        }else if (class(x) == "factor"){
            list(levels = levels(x))
        }
    })

    #check constraint for the initial pop
    K = length(g)
    if (K > 0){
        feasible = sapply(g, function(gg){
            gg(X)
        })
        if (K > 1){
            feasible = apply(feasible, 1, all)
        }
        X = X[feasible,]

        threshold = 5
        if (NROW(X) < threshold){
            stop(paste0("must provide at least", threshold, "feasible solution at the begenning"))
        }
    }

    X = X[1:min(2*N, NROW(X)), ]

    time_deb <- Sys.time()
    for (b in seq_len(B)){

        # 1) evaluate the function at X
        Y <- fn(X) %>% as.data.frame() %>% mutate(id = 1:NROW(.))


        # 2) Select the best indivdu base (Pt) on the domination notion
        Y$rank = dominance_ranking(Y[,-NCOL(Y)], sens)
        Y = crowding_distance(Y)
        Y = Y %>% arrange(rank, desc(crowding_distance)) %>%
          filter(crowding_distance != 0) %>%
          head(N)

        Pt <- X[Y$id,]


        # while not generate new feasible possible solution
        res_Qt = NULL
        counter = 1
        while (NROW(res_Qt) < N){

            # 3) Pick a couples to be reproduce after a tournement selection
            parents <- tournament_selection(Y, N) # ici il  ya un pb si NROW(X)>N

            # 4) Generate new individu (Qt) with genetic operator apply on each
            Qt <- crossing_over(X, parents = parents, crossing_over_size = crossing_over_size)
            Qt <- mutation(Qt, freq = freq_mutation, distri_Xi = distri_Xi) %>%
                as.data.frame()

            rownames(Qt)=1:NROW(Qt)

            #check that the new individu respect the g constraint
            if (K > 0){
                cat("check constraint:", counter, "\n")
                feasible = sapply(g, function(gg){
                    gg(Qt)
                })
                if (K > 1){
                    feasible = apply(feasible, 1, all)
                }
                res_Qt = rbind(res_Qt, Qt[feasible,])
                counter = counter + 1
            }else{
                res_Qt = Qt
            }

        }
        res_Qt = res_Qt[seq_len(N),]


        # 5) Gather the Pt and Qt and proceed to the next iteration
        X <- rbind(Pt, res_Qt)


        if (verbose) {
            cat("Population ", b, '\n')
            print(round((Sys.time() - time_deb), 1))
            cat("\n")
        }
    }


    #if crowding distance parameter thing
    no_domi = Y$rank == 1
    Y <- Y[no_domi, seq_len(n_objective)]
    Pt = Pt[no_domi,]

    #tri
    Y = Y %>% mutate(id = as.factor(1:NROW(Y))) %>% arrange_all()
    Pt = Pt[Y$id,]
    Y = Y %>% select(-id)

    # Y <- Y[, seq_along(fn)] #a remettre si crowding distance thing est supprime finalement.
    # colnames(Y) <- names(fn)
    rownames(Pt) <- seq_len(NROW(Pt))
    res = list(X = Pt, Y = Y, time = Sys.time() - time_deb, fn = fn, g = g,
               X0 = X0)
    class(res) = "nsga"
    res
}


#' plot nsga
#'
#' @export

plot.nsga <- function(res, alpha_ellipse = NULL){

    library(magrittr)
    cost_fn_0 = res$fn(res$X0)

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
            ggplot2::geom_point(data = cost_fn_0,
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
        if (is.numeric(res$X[, X_i])){
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


    n_y = length(p_y)
    n_x = length(p_x)

    s_y = plotly::subplot(p_y, titleY = TRUE, titleX = TRUE) #, nrows = ceiling(n_y/3))
    s_x = plotly::subplot(p_x, titleX = TRUE) #, nrows = ceiling(n_x/10))

    if (is.null(res$alpha)){
        pareto_title = "Pareto front"
    }else{
        pareto_title = paste0("Pareto front (alpha=", res$alpha, ")")
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

}
