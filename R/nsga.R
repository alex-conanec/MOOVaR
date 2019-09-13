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
#' @export

NSGA <- function(X, ff, N, crossing_over_size, freq_mutation, seed,
                 constraints = NULL, B=500, crossing_type = "hybrid",
                 verbose = TRUE, time_out = 5){

    require(dplyr)
    require(magrittr)

    if (!missing(seed)){
        set.seed(seed)
    }

    if (! crossing_type %in% c("hybrid", "auto"))
        stop("crossing_type must be either 'hybrid' or 'auto'")

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
        Y <- Y %>% mutate(rank = dominance_ranking(.[,-NCOL(Y)], sens))

        if (NCOL(Y) > 3){
            Y <- Y %>% bind_cols(crowding_distance = sapply(Y$rank %>% unique(), function(rank){
                crowding_distance(S = Y[Y$rank==rank, -((NCOL(Y)-1):NCOL(Y))])
            }) %>% unlist() %>% as.numeric()) %>% arrange(rank, desc(crowding_distance)) %>%
                head(N)
        }else{
            Y <- head(Y, N) %>% mutate(crowding_distance = runif(NROW(.)))
        }


        Pt <- X[Y$id,]

        # 3) Pick a couples to be reproduce after a tournement selection
        parents <- tournament_selection(Y, floor(NROW(X) - N)/2)

        # 4) Generate new individu (Qt) with genetic operator apply on each
        # couple, in respect to the set constraints
        start_loop <- Sys.time()
        repeat{

            if (crossing_type == "hybrid"){
                #crossing over
                Qt <- lapply(parents, function(couple){
                    crossing_over(couple = X[couple,],
                                  crossing_position = sample(seq_along(X),
                                                             crossing_over_size) %>%
                                      sort())
                }) %>% bind_rows()
            }else{
                #auto breeding
                Qt <- lapply(seq_len(N), function(i){
                    crossing_position = sample(seq_along(X), crossing_over_size)
                    structure(cbind(Pt[i, crossing_position], Pt[i, -crossing_position]),
                              colnames = colnames(Pt))
                    }) %>% bind_rows()
            }

            #mutation
            Qt <- mutation(Qt, freq = freq_mutation, distri_Xi = distri_Xi) %>%
                as.data.frame()

            #check for constraints
            if (!is.null(constraints)){
                passed_constraints <- sapply(seq_len(NROW(Qt)), function(i){
                    sapply(constraints, function(constraint){
                        constraint(as.data.frame(Qt[i,]))
                    }) %>% all()
                })
            }else{
                passed_constraints <- rep(TRUE, N)
            }

            Qt <- Qt[passed_constraints,]

            # out of the loop if everything is all right
            if (NROW(Qt) == N){
                break
            }

            # stop the algorithme if time out
            if (difftime(time1 = Sys.time(),
                         time2 = start_loop,
                         units = "min") > time_out){
                stop(paste("Time out when generating new population. Check if",
                     "the constraints are not too strick (even impossible)",
                     "or increase time_out arg"))
            }
        }

        # 5) Gather the Pt and Qt and proceed to the next iteration
        X <- rbind(Pt, Qt)


        if (verbose) {
            cat("Population ", b, '\n')
            print(round((Sys.time() - time_deb), 1))
            cat("\n")
        }
    }

    if (length(ff) > 1){
        Y <- Y[,seq_along(ff)]
        colnames(Y) <- names(ff)
    }else{
        Y <- Y[,1]
    }

    list(Pt = Pt, Y = Y)
}
