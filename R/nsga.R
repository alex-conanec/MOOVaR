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

NSGA <- function(X, fn, sens = rep("min", length(fn)), N = round(NROW(X) / 2, 0),
                 crossing_over_size = round(NCOL(X)/2),
                 freq_mutation=rep(1/NCOL(X), NCOL(X)),
                 seed = NULL, B=50, verbose = TRUE){

    if (!is.null(seed)) set.seed(seed)

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

        Y <- sapply(fn, function(f){
          if (class(f) == "function"){
            f(X)
          }else{
            predict(f, newdata = X)
          }
        }) %>% as.data.frame() %>% mutate(id = 1:NROW(.))

        # 2) Select the best indivdu base (Pt) on the domination notion
        Y$rank = dominance_ranking(Y[,-NCOL(Y)], sens)
        Y = crowding_distance(Y)
        Y = Y %>% arrange(rank, desc(crowding_distance)) %>%
          filter(crowding_distance != 0) %>%
          head(N)

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

    #if crowding distance parameter thing
    no_domi = Y$rank == 1
    Y <- Y[no_domi, seq_along(fn)]
    Pt = Pt[no_domi,]

    #tri
    Y = Y %>% mutate(id = as.factor(1:NROW(Y))) %>% arrange_all()
    Pt = Pt[Y$id,]
    Y = Y %>% select(-id)

    # Y <- Y[, seq_along(fn)] #a remettre si crowding distance thing est supprime finalement.
    colnames(Y) <- names(fn)
    rownames(Pt) <- seq_len(NROW(Pt))
    list(X = Pt, Y = Y, time = Sys.time() - time_deb)
}
