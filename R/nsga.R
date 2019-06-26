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
#' sum(1:10)

NSGA <- function(X, ff, N, crossing_over_size, freq_mutation,
                 B=500, verbose = TRUE){

    require(tidyverse)
    require(magrittr)
    require(rlist)

    if (missing(N)) N <- NROW(X) / 2
    sens <- sapply(ff, function(f){
        f$sens
    })

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
        Y %<>% mutate(rank = dominance_ranking(.[,-NCOL(Y)], sens))

        Y %<>% bind_cols(crowding_distance = sapply(Y$rank %>% unique(), function(rank){
            crowding_distance(S = Y[Y$rank==rank, -((NCOL(Y)-1):NCOL(Y))])
        }) %>% unlist() %>% as.numeric()) %>% arrange(rank, desc(crowding_distance)) %>%
            head(N)

        Pt <- X[Y$id,]

        # 3) Pick a couples to be reproduce after a tournement selection
        parents <- tournament_selection(Y, floor(NROW(X) - N)/2)

        # 4) Generate new individu (Qt) with genetic operator apply on each couple
        #crossing over
        Qt <-lapply(parents, function(couple){
            crossing_over(couple = X[couple,],
                          crossing_position = sample(1:4, crossing_over_size) %>% sort())
        }) %>% list.rbind()

        #mutation
        Qt <- mutation(Qt, freq = freq_mutation)

        # 5) Gather the Pt and Qt and proceed to the next iteration
        X <- rbind(Pt, Qt)


        if (verbose) cat("Population ", b, '\n')
    }
    list(Pt = Pt, Y = Y[,1:2])
}
