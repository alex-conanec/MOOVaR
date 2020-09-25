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
#' @export

crowding_distance <- function(Y, threshold = 10^(-5)*sapply(Y[,-((NCOL(Y)-1):NCOL(Y))], sd)){

    lapply(unique(Y$rank), function(r){
        a = Y %>% filter(rank == r)

        if (NROW(a) < 2){
            cbind(a, crowding_distance = rep((NCOL(Y)-2)*10^5, NROW(a)))
        }else{
            lapply(seq_len(NCOL(Y)-2), function(j){
                b = a[order(a[,j]),]
                nn = NROW(b)
                ecart_tot = b[NROW(b),j] - b[1,j]
                if (ecart_tot < threshold[j]){
                  b[,j] = c(10^5, rep(0, nn-1))
                }else{
                  b[,j] = c(10^5,
                            (b[-(1:2), j] - b[-((nn-1):nn), j])/ecart_tot,
                            10^5)
                }
                b[order(b$id),j]
            }) %>% unlist() %>%  matrix(ncol = NCOL(Y)-2) %>% rowSums() -> distance

          cbind(a, crowding_distance = distance)
        }
    }) %>% bind_rows() %>% arrange(id)

}
