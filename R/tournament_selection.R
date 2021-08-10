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
tournament_selection <- function(Y, k, N = NROW(Y)){

    # k = 4
    # N=30
    # p=2
    # Y = matrix(runif(N*p), ncol = p, nrow = N)
    # Y = as.data.frame(Y)
    # Y$rank = dominance_ranking(Y, sens = rep("min", 2))
    # Y$id = 1:N
    # Y = crowding_distance(Y)

    idx = matrix(sample(1:NROW(Y), k*N, replace = TRUE), ncol = k,
               byrow = TRUE, nrow = N)

    apply(idx, 1, function(x){
        A = matrix(Y[x,]$rank, ncol = k, nrow = k, byrow = TRUE)
        have_best_rank = apply(A - t(A) <= 0, 2, all)
        if (sum(have_best_rank) > 1){
            k_exeaquo = sum(have_best_rank)
            B = matrix(Y[x,][have_best_rank, ]$crowding_distance,
                       ncol = k_exeaquo, nrow = k_exeaquo, byrow = TRUE)
            have_best_CD = apply(B - t(B) >= 0, 2, all)
            Y[x,][have_best_rank, ][have_best_CD,]$id[1]
        }else{
            Y[x,][have_best_rank, ]$id
        }
    }) %>% unlist()

}
