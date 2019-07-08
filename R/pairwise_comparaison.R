#' Pairwise comparison
#'
#' Compare one by one all the individu evaluations and give them a rank and
#' count how many individu dominate each individus and the individus that are
#' dominated by each individus
#'
#' @param X a matrix/data.frame of all the selected individus evaluation
#'
#' @param sens the objectif fonctions goal which can be either "min" or "max"
#'
#' @return a list of two element for each individu. First the dominated_count
#' which mean how many individu are better. Secondly the dominating_index
#' which gather the index of the dominted individu.
#'
#' @examples
#' sum(1:10)

pairwise_comparison <- function(X, sens){

    if (any(!sens %in% c("min", "max"))) stop("sens must be either 'min' or 'max'")
    dom_sens <- c('>', '<')[as.factor(sens=="min")]
    res <- list()

    for (k in 1:NROW(X)){
        res[[as.character(k)]]$dominated_count <- 0
        res[[as.character(k)]]$dominating_index <- NULL

        for (i in (1:NROW(X))[-k]){
            comparison <- sapply(1:NCOL(X), function(j){
                outer(X[k,j], X[i,j], dom_sens[j])
            })

            equals <- sapply(1:NCOL(X), function(j){
                outer(X[k,j], X[i,j], '==')
            })

            if (!(any(comparison) & any(!comparison)) & !any(equals)){ #comparable
                if (!all(comparison)){
                    res[[k]]$dominated_count <- res[[k]]$dominated_count + 1
                }else if (all(comparison)){
                    res[[k]]$dominating_index <- c(res[[k]]$dominating_index, i)
                }
            }
        }
    }
    res
}

# pairwise_comparison <- function(X, sens){
# require(tidyverse)
#     dom_sens <- c('<', '>')[as.factor(sens=="min")]
#
#     result <- lapply(1:NROW(X), function(k){
#         res <- list(dominated_count = 0, dominating_index = NULL)
#
#         comparisons <- lapply((1:NROW(X))[-k], function(i){
#             sapply(1:NCOL(X), function(j){
#                 outer(X[k,j], X[i,j], dom_sens[j])
#             })
#         })
#
#         res$dominated_count <-
#             sapply(comparisons, function(comparison){
#                 all(!comparison)
#             }) %>% which() %>% length()
#
#         res$dominating_index <-
#             sapply(1:length(comparisons), function(i){
#                 if (all(comparisons[[i]])) (1:NROW(X))[-k][i]
#             }) %>% unlist()
#
#         res
#     })
#     result
# }
