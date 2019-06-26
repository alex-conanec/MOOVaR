#' Dominance ranking
#'
#' Use the pairwisecomparison function to attribute a rank of individu belonging
#' to the same front which mean that they are not comparable to each other.
#' The lowest the rank the best is the individu
#'
#' @param X a matrix/data.frame of all the selected individus evaluation
#'
#' @param sens the objectif fonctions goal which can be either "min" or "max"
#'
#' @return a vecteur of rank associated with the individus of X
#'
#' @examples
#' sum(1:10)

dominance_ranking <- function(X, sens){

    require(tidyverse)
    res <- rep(NA, NROW(X))
    rank <- 1
    PwCp <- pairwise_comparison(X, sens)
    last_points <- 1:length(PwCp)

    while (any(is.na(res))){
        temp <-
            sapply(last_points, function(i){
                if (PwCp[[i]]$dominated_count == 0) i
            }) %>% unlist()

        res[temp] <- rank

        for (j in temp){
            for (k in PwCp[[j]]$dominating_index){
                PwCp[[k]]$dominated_count <- PwCp[[k]]$dominated_count - 1
            }
        }

        rank <- rank + 1
        last_points <- last_points[!last_points %in% temp]
    }
    res
}


# dominance_ranking <- function(X, sens){
#
#     res <- rep(NA, NROW(X))
#     rank <- 1
#     last_points <- 1:length(PwCp)
#     PwCp <- pairwise_comparison(X, sens)
#
#     while (any(is.na(res))){
#         temp <- NULL
#         for (i in last_points){
#             if (PwCp[[i]]$dominated_count == 0){
#                 temp <- c(temp, i)
#             }
#         }
#
#         res[temp] <- rank
#
#         for (j in temp){
#             for (k in PwCp[[j]]$dominating_index){
#                 PwCp[[k]]$dominated_count <- PwCp[[k]]$dominated_count - 1
#             }
#         }
#
#         rank <- rank + 1
#         last_points <- last_points[!last_points %in% temp]
#     }
#     res
# }



