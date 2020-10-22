#' def_cstr_X_space
#'
#' Find the space where points where observed. Send back a function to check if
#' a/several point(s) belong to the space where the points where observed
#'
#'  @param X matrix of dimension n x d observed points
#'  @param m a vector of size d to indicate how many reference points create for each dimension.
#'
#' @examples
#' library(tidyverse)
#'
#' #generate data ---
#' #parameter of the data
#' n = 100
#' d = 2 #dim de X
#'
#' mini = rep(0, d)
#' maxi = rep(10, d)
#'
#' set.seed(123)
#'
#' #uniform
#' X = sapply(seq_len(d), function(i){
#'   runif(n, mini[i], maxi[i])
#' })
#' plot(X)
#'
#' #avec des trous
#' b = matrix(c(0, 0, 5, 5), byrow = T, ncol = d)
#' B = matrix(c(4, 4, 10, 10), byrow = T, ncol = d)
#' p = NROW(B) #nombre de "foyer"
#'
#' X_c = sapply(seq_len(d), function(i){
#'   sapply(seq_len(p), function(k){
#'     runif(n/p, b[k,i], B[k,i])
#'   })
#' })
#' plot(X_c)
#'
#' #use function
#' res = def_cstr_X_space#check several points
(X_c)
#'
#' #plot reference points with two colors depending or their feasability
#' plot(res$ref_points, pch=3, col = "orange")
#' points(res$feasible_ref_points, pch=3, col = "blue")
#'
#' #check only one point
#' i=3
#' points(X[i,1], X[i,2], col = c("red", "green")[as.numeric(res$g(X = X[i,])) + 1])
#'
#' #check several points
#' points(X, col = c("red", "green")[as.numeric(res$g(X)) + 1])
#'
#' @export

def_cstr_X_space <- function(X, m = rep(10, NCOL(X)) ){

  n = NROW(X)
  mini = floor(apply(X, 2, min))
  maxi = ceiling(apply(X, 2, max))
  h = (maxi - mini)/m

  #all the combination possible (depending of m)
  K = prod(m+1)
  ref_points = matrix(unlist(purrr::cross(lapply(seq_len(d), function(j){
    sapply(seq_len(m[j] + 1) - 1, function(k){
      mini[j] + k*h[j]
    })
  }))), ncol = d, byrow = T)

  #generate a XX and a AA to have a conformable dimension
  XX = array(X, dim = c(n, d, K))
  AA = array(ref_points, c(K, d, n)) %>%
    apply(c(3, 2, 1), FUN = function(x) x)

  #the points from which you have to be close to be in the feasible space
  feasible_ref_points = apply(XX - AA, c(1,3), function(x) x[1]^2 + x[2]^2) %>% #distance euclidienne
      apply(1, function(x) x == min(x)) %>%  #trouve point le plus proche
      apply(1, any) # combinaison de A ayant au moins un X proche


  #function that check if point(s) X is/are feasible
  g = function(X, ref_points, feasible_ref_points, K){

    d = NCOL(ref_points)
    X = matrix(X, ncol = d)
    n = NROW(X)

    XX = array(X, dim = c(n, d, K))
    AA = array(ref_points, c(K, d, n)) %>%
      apply(c(3, 2, 1), FUN = function(x) x)

    closest_ref_points = apply(XX - AA, c(1,3), function(x) x[1]^2 + x[2]^2) %>%
      apply(1, function(x) x == min(x))

    res = closest_ref_points[feasible_ref_points, ]

    if (length(dim(res)) > 1){
      apply(res, 2, any)
    }else{
      any(res)
    }
  }

  formals(g)$ref_points = ref_points
  formals(g)$feasible_ref_points = feasible_ref_points
  formals(g)$K = K

  list(
    ref_points = ref_points,
    feasible_ref_points = ref_points[feasible_ref_points,],
    g = g
  )

}
