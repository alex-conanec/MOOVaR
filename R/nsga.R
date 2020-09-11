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
#' set.seed(123)
#' n <- 300
#' q = 3
#' mini <- rep(x = -2, times = q)
#' maxi <- rep(x = 2, times = q)
#'
#' X <- lapply(seq_len(q), function(k){
#'     res <- data.frame(X=runif(n = n, min = mini[k], max = maxi[k]))
#'     names(res) <- paste0("X", k)
#'     res
#' }) %>% bind_cols()
#'
#' ff <- list(
#'     list(f = function(X) unlist(X[1] + 2*X[2] + X[3]),
#'          sens = "min"),
#'     list(f = function(X) unlist(-X[1] - X[2] - 3*X[3]),
#'          sens = 'min')
#' )
#'
#' library(NSGA2)
#' NSGA(X=X[1:50,], ff, N = 50, crossing_over_size = 2, freq_mutation = rep(0.3, q),
#'     seed = 123, constraints = NULL, B=50, crossing_type = "hybrid",
#'     verbose = TRUE, time_out = 5)
#'
#' @export

NSGA <- function(X, ff, N, crossing_over_size, freq_mutation, seed,
                 B=500, verbose = TRUE){

    require(dplyr)
    require(magrittr)

    if (!missing(seed)){
        set.seed(seed)
    }

    if (missing(N)) N <- round(NROW(X) / 2, 0)

    sens <- sapply(ff, function(f){
        f$sens
    })

    distri_Xi <- lapply(X, function(x){
        if (class(x) == "numeric"){
            list(min = floor(min(x)), max = round(max(x)),
                 mean = mean(x), sd = sd(x))
        }else if (class(x) == "factor"){
            list(levels = levels(x))
        }
    })

    time_deb <- Sys.time()
    for (b in seq_len(B)){

        # 1) Evaluation of the performances of each individu for each objectif function
        Y <- sapply(1:length(ff), function(j){

            if (class(ff[[j]]$f) == "function"){
                ff[[j]]$f(X)
            }else{
                predict(ff[[j]]$f, newdata = X)
            }
        }) %>% as.data.frame() %>% mutate(id = 1:NROW(.))

        # 2) Select the best indivdu base (Pt) on the domination notion
        Y$rank = dominance_ranking(Y[,-NCOL(Y)], sens)
        Y = crowding_distance(Y)
        Y = Y %>% arrange(rank, desc(crowding_distance)) %>% head(N)

        Pt <- X[Y$id,]

        # 3) Pick a couples to be reproduce after a tournement selection
        parents <- tournament_selection(Y, N) # ici il  ya un pb si NROW(X)>N

        # 4) Generate new individu (Qt) with genetic operator apply on each
        Qt <- crossing_over(X, parents = parents, crossing_over_size = crossing_over_size)
        Qt <- mutation(Qt, freq = freq_mutation, distri_Xi = distri_Xi) %>%
            as.data.frame()


        # 5) Gather the Pt and Qt and proceed to the next iteration
        X <- rbind(Pt, Qt)


        if (verbose) {
            cat("Population ", b, '\n')
            print(round((Sys.time() - time_deb), 1))
            cat("\n")
        }
    }

    Y <- Y[,seq_along(ff)]
    colnames(Y) <- names(ff)
    rownames(Pt) <- seq_len(NROW(Pt))
    list(X = Pt, Y = Y)
}
