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

    if (class(X) == "numeric") X <- data.frame(X)

    b = lapply(seq_len(NCOL(X)), function(j){

        A = matrix(X[,j], nrow = NROW(X), ncol = NROW(X), byrow = F)
        outer(A - t(A), 0, dom_sens[j])

    })

    m = array(unlist(b), dim = c(NROW(X), NROW(X), NCOL(X)))

    domination = apply(m, 1:2, all)

    lapply(seq_len(NROW(X)), function(i){
        list(
            dominating_index = which(domination[i,]),
            dominated_count = length(which(domination[,i]))
            )
    })

}
