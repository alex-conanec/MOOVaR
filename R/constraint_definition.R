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
#' res = def_cstr_X_space(X_c)#check several points
#'
#' #plot reference points with two colors depending or their feasability
#' plot(res$ref_points, pch=3, col = "orange")
#' points(res$feasible_ref_points, pch=3, col = "blue")
#' points(X_c)
#'
#' #check the feasibility
#' plot(res$ref_points, pch=3, col = "orange")
#' points(res$feasible_ref_points, pch=3, col = "blue")
#'
#' #of only one point
#' i=3
#' points(X[i,1], X[i,2], col = c("red", "green")[as.numeric(res$g(X = X[i,])) + 1])
#'
#' #of several points
#' points(X, col = c("red", "green")[as.numeric(res$g(X)) + 1])
#'
#' #visualize the combination of variable. Very useful when d > 2
#' plot(res)
#'
#' @export

def_cstr_X_space <- function(X, m = rep(10, NCOL(X)) ){

  n = NROW(X)
  d = NCOL(X)
  mini = floor(apply(X, 2, min))
  maxi = ceiling(apply(X, 2, max))
  h = (maxi - mini)/m

  #all the combination possible (depending of m)
  K = prod(m)
  ref_points = matrix(unlist(purrr::cross(lapply(seq_len(d), function(j){
    sapply(seq_len(m[j])-0.5, function(k){
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

  res = list(
    ref_points = ref_points,
    feasible_ref_points = ref_points[feasible_ref_points,],
    g = g,
    h = h
  )

  class(res) = "combi_space"
  res
}

#' @export
plot.combi_space <- function(X, m = rep(3, NCOL(X$feasible_ref_points))){

  df = lapply(seq_along(X$h), function(j){
    if (m[j] == 3){
      l = seq(-1, 1, length.out = m[j])
    }else{
      l = 0
    }
    lapply(l, function(k){
      data.frame(X = X$feasible_ref_points[,j] + k*X$h[j]/2,
                 id = 1:NROW(X$feasible_ref_points))
    }) %>% bind_rows()
  }) %>% purrr::reduce(full_join, by = "id") %>%
    gather(key = X_j, value = coord, -id) %>%
    mutate(id = as.factor(id),
           X_j = as.factor(X_j))

  dd <- plotly::highlight_key(df, ~id)
  p = ggplot(dd, aes(x = X_j, y = coord)) +
    geom_point(aes(group = id)) +
    facet_wrap(~ X_j)

  plotly::ggplotly(p, tooltip = "coord") %>%
    highlight(on = "plotly_selected", color = "red")

}






##### interpolation hypothese basé sur absence d'interations entre deux X pour une fonction donnee... mais sur un espace observé seulement...
# f1 = function(x){
#   x[1] + x[2]
# }
#
# f2 = function(x){
#   x[1] + x[2] + x[1] * x[2]
# }
#
# f3 = function(x){
#   if ((x[1] < B[2,1] & x[2] < B[2,2] & x[1] > b[2,1] & x[2] > b[2,2]) |
#       (x[1] > b[1,1] & x[2] > b[1,2] & x[1] < B[1,1] & x[2] < B[1,2])){
#     x[1] + x[2] + x[1] * x[2]
#   }else{
#     x[1] + x[2] - 3*x[1] * x[2]
#   }
# }
#
#
# Y = sapply(list(f1, f2, f3), function(f){
#   sapply(seq_len(NROW(X)), function(i){
#     f(X[i,])
#   })
# })
#
#
# Y_c = sapply(list(f1, f2, f3), function(f){
#   sapply(seq_len(NROW(X_c)), function(i){
#     f(X_c[i,])
#   })
# })
#
#
# #on regarde sur l'ensemble de R2 maintenant dans le cas ou on a f3
# Z = data.frame(X, Y = Y) %>% gather(key = f, value = Y, -X1, -X2) %>%
#   mutate(f = as.factor(f)) %>% filter(f == "Y.3")
#
#
# plot3D::scatter3D(x = Z$X1, y = Z$X2, z = Z$Y)
# plot3Drgl::plotrgl()
#
#
# fit1 = lm(Y_c[,1] ~ X_c[,1]*X_c[,2])
# summary(fit1)
#
# fit2 = lm(Y_c[,2] ~ X_c[,1]*X_c[,2])
# summary(fit2)
#
# df = data.frame(y=Y_c[,3], X_c)
# deg = 3
# formule = paste0("y ~ ", "(", paste(colnames(df)[-1], collapse = " + "), ")^", deg)
# fit3 = lm(formule, data = df)
#
# #on se base sur p.value associé a chaque beta pour le model globale
# summary(fit3)$coefficients[-1, "Pr(>|t|)"]
#
# #on se base sur le meilleur model --> selection varable
# best_fit = step(fit3)
# summary(best_fit)$coefficients[-1, "Pr(>|t|)"]


