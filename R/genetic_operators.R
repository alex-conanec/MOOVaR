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

mutation <- function(Qt, freq, distribution = "uniform", distri_Xi){

    if (missing(Qt)) stop("Qt is missing")
    if (missing(freq)) stop("freq is missing")
    if (length(freq) != NCOL(Qt)) stop("You have to assign a frequence of mutation for
                                each variable")

    for (j in seq_len(NCOL(Qt))){
        i <- sample(seq_len(NROW(Qt)), size = ceiling(NROW(Qt)*freq[j]))
        var <- names(Qt)[j]
        if (class(Qt[1,j]) == "numeric"){
            if (distribution=="uniform"){
                min_j <- distri_Xi[[var]]$min
                max_j <- distri_Xi[[var]]$max
                Qt[i, j] <- runif(length(i), min = min_j, max = max_j)
            }else if (distribution=="normal"){
                mean_j <- distri_Xi[[var]]$mean
                sd_j <- distri_Xi[[var]]$sd
                Qt[i, j] <- rnorm(length(i), mean = mean_j, sd = sd_j)
            }else{
                stop("Distribution must be either 'uniform' or 'normal'")
            }
        }else if (class(Qt[1,j]) == "factor"){
            levels_j <- distri_Xi[[var]]$levels
            Qt[i, j] <- sample(x = levels_j,
                               size = length(i),
                               replace = TRUE)
        }
    }
    Qt
}
