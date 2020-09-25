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


#'no_dominated
#' @export
no_dominated <- function(front, X, sens){

  X <- as.matrix(X)
  if (any(!sens %in% c("min", "max"))) stop("sens must be either 'min' or 'max'")

  #intrapolation front
  nb_point_min = 1000
  front = unique(front)
  n = NROW(front)
  p = NCOL(front)
  pas = ( max(front[,1]) - min(front[,1]) ) /
    nb_point_min * n
  k = 0
  for (i in 2:n){
    sign_grad = sign(front[i+k, 1] - front[i-1+k, 1])
    Y1 = seq(front[i-1+k, 1], front[i+k, 1], sign_grad * pas)[-1]
    if (length(Y1) > 0){
      Y_autre_col =
        matrix( unlist(front[i+k-1, 2:(p)]), ncol = p-1, nrow = length(Y1), byrow = T) +
        matrix( unlist((Y1 - front[i+k-1, 1]) / (front[i+k, 1] - front[i+k-1, 1])),
                ncol = (p-1), nrow = length(Y1)) *
        matrix( unlist(front[i+k, 2:(p)] - front[i+k-1, 2:(p)]),
                ncol = p-1, nrow = length(Y1), byrow = T)
      colnames(Y_autre_col) = colnames(front)[2:p]
      front = rbind(front[1:(i+k-1),], cbind(Y1, Y_autre_col), front[(i+k):(n+k),])
      k = k + length(Y1)
    }
  }

  b = lapply(seq_len(NCOL(X)), function(j){

    A = matrix(X[,j], nrow = NROW(X), ncol = NROW(front), byrow = F)
    B = matrix(front[,j], nrow = NROW(X), ncol = NROW(front), byrow = T)
    A - B
  })

  delta = array(unlist(b), dim = c(NROW(X), NROW(front), NCOL(X)))
  delta = apply(delta, c(3,1,2), function(x) x)

  min_dom_strict = apply((delta < 0)[which(sens == "min"),, ,drop = FALSE], 2:3, any)
  max_dom_strict = apply((delta > 0)[which(sens == "max"),, ,drop = FALSE], 2:3, any)

  if (all(sens == "min")){
    apply(min_dom_strict, 1, all)
  }else if (all(sens == "max")){
    apply(max_dom_strict1, all)
  }else{
    apply(min_dom_strict | max_dom_strict, 1, all)
  }
}
