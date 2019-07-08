#' Tournement selection
#'
#' Select a couple of individu (parents) to the crossing over
#'
#' @param X a matrix/data.frame of all the selected individus evaluation
#'
#' @param N size of the Qt population to generate
#'
#' @return a list of couple index which must be use for crossing over
#'
#' @examples
#' sum(1:10)

tournament_selection <- function(X, N){

    lapply(seq_len(N), function(i){
        temp <- sample(1:NROW(X), 4, replace = FALSE)
        sapply(c(1,3), function(j){
            if (X[temp[j],]$rank < X[temp[j+1],]$rank){
                temp[j]
            }else if (X[temp[j],]$rank > X[temp[j+1],]$rank){
                temp[j+1]
            }else if (X[temp[j],]$crowding_distance > X[temp[j+1],]$crowding_distance){
                temp[j]
            }else if (X[temp[j],]$crowding_distance < X[temp[j+1],]$crowding_distance){
                temp[j+1]
            }else{
                temp[j]
            }
        })
    })
}
