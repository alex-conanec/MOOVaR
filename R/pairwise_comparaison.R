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
#' @export

pairwise_comparison <- function(X, sens){

    if (any(!sens %in% c("min", "max"))) stop("sens must be either 'min' or 'max'")
    dom_sens <- sens
    dom_sens[which(sens == "min")] <- '<'
    dom_sens[which(sens == "max")] <- '>'
    res <- list()

    if (class(X) == "numeric") X <- data.frame(X)
    for (k in 1:NROW(X)){
        res[[as.character(k)]]$dominated_count <- 0
        res[[as.character(k)]]$dominating_index <- NULL

        for (i in (1:NROW(X))[-k]){

            comp <- comparaison(X1 = X[k,],
                                X2 = X[i,],
                                sens = sens)

            if (!is.null(comp)){
                if (comp){
                    res[[k]]$dominating_index <- c(res[[k]]$dominating_index, i)
                }else{
                    res[[k]]$dominated_count <- res[[k]]$dominated_count + 1
                }
            }
        }
    }
    res
}
