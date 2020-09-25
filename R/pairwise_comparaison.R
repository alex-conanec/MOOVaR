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

    A = array(as.matrix(X), dim = c(NROW(X), NCOL(X), NROW(X)))
    A_T = apply(A, c(3,1), t)
    B = apply(A, c(2,1,3), function(x) x)
    m = B - A_T


    min_dom_strict = apply((m < 0)[which(sens == "min"),, ,drop = FALSE], 2:3, any)
    min_dom_eq = apply((m <= 0)[which(sens == "min"),, ,drop = FALSE], 2:3, all)

    max_dom_strict = apply((m > 0)[which(sens == "max"),, ,drop = FALSE], 2:3, any)
    max_dom_eq = apply((m >= 0)[which(sens == "max"),, ,drop = FALSE], 2:3, all)

    if (all(sens == "min")){
      dom = min_dom_strict & min_dom_eq
    }else if (all(sens == "max")){
      dom = max_dom_strict & max_dom_eq
    }else{
      dom = (min_dom_strict | max_dom_strict) & max_dom_eq & min_dom_eq
    }

    lapply(seq_len(NROW(X)), function(i){
        list(
            dominating_index = which(dom[i,]),
            dominated_count = length(which(dom[,i]))
            )
    })

}
