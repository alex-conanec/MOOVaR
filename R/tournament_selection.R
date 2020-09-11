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
#' @export

tournament_selection <- function(Y, N){

    a = matrix(sample(1:N, 2*N, replace = TRUE), ncol = 2, byrow = T, nrow = N)

    diff_rank = Y$rank[a[,1]] - Y$rank[a[,2]]

    equal = diff_rank == 0

    diff_crowding_distance = Y$crowding_distance[a[,1]] - Y$crowding_distance[a[,2]]

    matrix(c(a[diff_rank < 0, 1], a[diff_rank > 0, 2],
             a[diff_crowding_distance >= 0 & equal, 1], a[diff_crowding_distance < 0 & equal, 2]),
           ncol = 2, byrow = TRUE)
}
