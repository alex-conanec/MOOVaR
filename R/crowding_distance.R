# il y a apparement des amaliroation possible
# https://www.atlantis-press.com/journals/ijcis/125905769/view


#' Crowding distance
#'
#' Compute the crowding distance between all the point belonging to a same
#' front
#'
#' @param S a matrix/data.frame of evaluated points belonging to the same front
#'
#' @return A vector containing the crowding distances
#'
#' @examples
#' sum(1:10)

crowding_distance <- function(S){
    require(tidyverse)
    distance <- rep(0, NROW(S))
    names(distance) <- 1:NROW(S)
    for (j in 1:NCOL(S)){
        temp <- sort(S[,j])
        id <- data.frame(id = 1:NROW(S), value = S[,j]) %>%
            arrange(value) %>% pull(id)
        interval <- max(temp) - min(temp)
        if (interval == 0) interval <- 10^20
        distance[id[1]] <- 10^50
        distance[id[NROW(S)]] <- 10^50

        if (NROW(S) > 2){
            for (i in 2:(NROW(S)-1)){
                distance[id[i]] <- distance[id[i]] + (temp[i+1]-temp[i-1])/interval
            }
        }
    }
    distance
}
