# function that mesure the aire between two lines where only few points of it
# are known

#' front_diversity
#'
#' Function to evaluate the diversity of a set of pareto front which aim to
#' reach the same optimal but which reach different fronts
#'
#' @param all_front a numerical matrix n x p, n being the number of point
#' belonging at all the differents pareto front.
#'
#' @param sens vector of size p with value in c("max", "min") giving the objectif
#' seeking for each component.
#'
#' @param n number of point sample to proceed to the Monte Carlos estimation
#'
#' @return return a numerical value, being a area for p=2, a volume if p=3, and
#' so on, corresponding to the diversity of the fronts
#'
#' @export

front_diversity <- function(all_front, sens, n = 10^4){
    all_front <- unique(all_front)
    anti_sens <- sens
    anti_sens[which(sens=="max")] <- "min"
    anti_sens[which(sens=="min")] <- "max"

    boundary_high <- pareto_front(X=all_front, sens = sens)
    boundary_low <- pareto_front(X=all_front, sens = anti_sens)

    area_btw_front(front_sup = boundary_high,
                   front_inf = boundary_low,
                   sens = sens,
                   n = n)
}

#' area_btw_front
#'
#' Estimate the area between two front of dimension p throught a Monte Carlos
#' algorithme.
#'
#' @param front_sup set of n point having p coordinates store in a matrix (n x p).
#' front_sup is better than front_inf as define by sens
#'
#' @param front_inf set of n point having p coordinates store in a matrix (n x p).
#' front_inf is worst than front_sup as define by sens
#'
#' @param sens vector of size p with value in c("max", "min") giving the objectif
#' seeking for each component.
#'
#' @param n number of point sample to proceed to the Monte Carlos estimation
#'
#' @return return a numerical value, being a area for p=2, a volume if p=3, and
#' so on
#'
#' @export

area_btw_front <- function(front_sup, front_inf, sens, n = 10^4){

    anti_sens <- sens
    anti_sens[which(sens=="max")] <- "min"
    anti_sens[which(sens=="min")] <- "max"

    min_point <- sapply(seq_len(NCOL(front_sup)), function(j){
        min(c(front_sup[,j], front_inf[,j]))
    })

    max_point <- sapply(seq_len(NCOL(front_sup)), function(j){
        max(c(front_sup[,j], front_inf[,j]))
    })

    spaceMC <- sapply(seq_along(min_point), function(j){
        runif(n=n, min = min_point[j], max = max_point[j])
    })

    above_sup <- count_under_front(front = front_sup,
                               X = spaceMC,
                               sens = sens)

    above_inf <- count_under_front(front = front_inf,
                               X = spaceMC,
                               sens = sens)

    area_covered <- prod(max_point - min_point)
    area_covered * (above_sup - above_inf) / n

    # below <- count_under_front(front = front_inf,
    #                            X = spaceMC,
    #                            sens = anti_sens)
    #
    # area_covered <- prod(max_point - min_point)
    # (n - above - below)/n * area_covered

}

#' pareto_front
#'
#' Find the pareto front inside a set of point
#'
#' @param X a numerical matrix n x p of the set of point
#'
#' @param sens a character vector of size p allowing only the two value in
#' c("min", "max"). A jth "min" value of sens means that the function
#' seeking the pareto front in order to minmize the jth component
#'
#' @return return a numerical matrix corresponding to the set of point belonging
#' to pareto front
#'
#' @export

pareto_front <- function(X, sens){

    X <- as.matrix(X)
    res <- list(X[1,])

    for (i in seq_len(NROW(X))[-1]){
        for (k in seq_along(res)){
            comp <- comparaison(X1 = X[i,], X2 = res[[k]], sens)
            if (!is.null(comp)){
                if (comp){
                    res[[k]] <- X[i,]
                    break
                }else{
                    break
                }
            }
        }
        if (is.null(comp)){
            res[[length(res) + 1]] <- X[i,]
        }
    }
    matrix(unlist(res), ncol = NCOL(X), byrow = TRUE)
}

#' count_under_front
#'
#' Function that count the number of point below/above the front given priorities
#'
#'
#' @param front a numerical matrix n x p, n being the number of point forming
#' the front.
#'
#' @param X a numerical matrix n2 x p, n2 being the number of position points
#' to test
#'
#' @param sens a character vector of size p allowing only the two value in
#' c("min", "max"). A "min" value at the j index of sens means that the function
#' will count the number of point below the front of the j component and
#' vice versa.
#'
#' @return return the number of point
#'
#' @export

count_under_front <- function(front, X, sens){
    X <- as.matrix(X)
    res <- 0
    for (i in seq_len(NROW(X))){
        for (k in seq_len(NROW(front))){
            comp <- comparaison(X1 = X[i,], X2 = front[k,], sens)
            if (!is.null(comp)){
                if (comp){
                    res <- res + 1
                    break
                }else{
                    break
                }
            }
        }
    }
    res
}


#' comparaison
#'
#' Function to compare to vector on all their component given a set of
#' priorities
#'
#' @param X1 a numeric vector of size p which will be compare to the second one
#' @param X2 a numeric vector of size p which will be compare to the first one
#' @param sens a character vector of size p allowing only the two value in
#' c("min", "max"). A value "min" of sens at the j position means to test if
#' X1 < X2.
#'
#' @return return TRUE if X1 is better than X2 on all the component given the
#' set of instruction sens, FALSE if it is the opposite and NULL if for all the
#' component none between X1 and X2 is better.
#'
#' @export

comparaison <- function(X1, X2, sens){

    if (any(!sens %in% c("min", "max"))) stop("sens must be either 'min' or 'max'")
    dom_sens <- sens
    dom_sens[which(sens == "min")] <- '<'
    dom_sens[which(sens == "max")] <- '>'

    X1 <- as.matrix(X1)
    comparison <- sapply(seq_along(X1), function(j){
        outer(X1[j], X2[j], dom_sens[j])
    })

    equals <- sapply(seq_along(X1), function(j){
        outer(X1[j], X2[j], "==")
    })

    if (!(any(comparison) & any(!comparison)) & !any(equals)){ #comparable
        if (all(comparison)){
            TRUE
        }else{
            FALSE
        }
    }else NULL
}
