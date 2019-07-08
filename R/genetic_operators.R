#' Crossing over operator
#'
#' This function switch components for the given variables position between two
#' individus (parent) which have been selected after the tounament
#'
#' @param couple a matrix/data.frame of two individus
#' @param crossing_position a sorted vector of variable position to switch
#'
#' @return a matrix/data.frame of the two childs from the crossing over of the
#' parents
#'
#' @examples
#' sum(1:10)

crossing_over <- function(couple, crossing_position){

    temp <- couple[1, crossing_position]
    couple[1, crossing_position] <- couple[2, crossing_position]
    couple[2, crossing_position] <- temp
    couple
}

#' Mutation operator
#'
#' This function assign random value to the individu caracteristique following
#' a frequency of occurence
#'
#' @param Qt a matrix/data.frame of all the individus generated after crossing
#' over
#' @param freq a vector of numerics value varying from 0 to 1 which represent
#' the frequency of mutation of each variable
#'
#' @param distribution choice of the distribution where for the mutation value
#'
#' @return a matrix/data.frame of the two childs from the crossing over of the
#' parents
#'
#' @examples
#' sum(1:10)

mutation <- function(Qt, freq, distribution = "uniform"){

    if (missing(Qt)) stop("Qt is missing")
    if (missing(freq)) stop("Qt is missing")
    if (length(freq) != NCOL(Qt)) stop("You have to assign a frequence of mutation for
                                each variable")

    sapply(seq_len(NCOL(Qt)), function(j){
        i <- sample(seq_len(NROW(Qt)), size = ceiling(NROW(Qt)*freq[j]))
        if (distribution=="uniform"){
            Qt[i, j] <- runif(length(i), min = min(Qt[,j]), max = max(Qt[,j]))
        }else if (distribution=="normal"){
            Qt[i, j] <- rnorm(length(i), mean = mean(Qt[,j]), sd = sd(Qt[,j]))
        }else{
            stop("Distribution must be either 'uniform' or 'normal'")
        }
        Qt[,j]
    })
}
