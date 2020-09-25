#' hypervolume ellipse
#' @export
hypervolume_ellipse <- function(Y, alpha = 0.05){

  x = eigen(var(Y))$values * qnorm(1-alpha/2)
  n = length(x)
  pi^(n/2)*prod(x)/gamma(n/2 + 1)

}

#' ellipse_var
#' @export
ellipse_var = function(Y1, Y2, alpha = 0.01, N=100){

  Y = data.frame(Y1, Y2)
  scale_ellipse = sapply(Y, sd)
  center_ellipse = sapply(Y, mean)

  dd=eigen(var(Y))

  a = sqrt(dd$values[1]) * qnorm(1-alpha/2)
  b = sqrt(dd$values[2]) * qnorm(1-alpha/2)

  x = seq(-a, a, length.out = N)
  y = sqrt(b^2*(1-x^2/a^2))

  ellipse = as.matrix(data.frame(x = c(x, x[N:1]), y = c(y, (-y)[N:1])))

  ellipse = ellipse %*% solve(dd$vectors)
  ellipse[,1] = (ellipse[,1] + center_ellipse[1])
  ellipse[,2] = (ellipse[,2] + center_ellipse[2])

  ellipse

}


#' sensibilite_nsga
#' @export
sensibilite_nsga = function(X, fn, sens = rep("min", length(fn)), N_MC = 10^4,
                            B = 50, N = 50){

  # pour paraleliser un peu
  numCores <- parallel::detectCores()

  #definition de la grille de recherche des parametres
  crossing_over_size_choice = seq_len(round(NCOL(X)/2))
  freq_mutation_choice = seq(0, 1, 0.1)

  combi = matrix(unlist(purrr::cross2(crossing_over_size_choice,
                                      freq_mutation_choice)),
                 ncol = 2, byrow = T) %>% as.data.frame()

  colnames(combi) = c("crossing_over", "mutation")

  #parametre de NSGA a optimiser aussi peut etre
  X0 = X

  #recherche du front avec les parametres choisis
  parallel::mclapply(seq_len(NROW(combi)), function(i){

    res = NSGA(X = X0, fn, N = N, crossing_over_size = combi[i,1],
               freq_mutation = rep(combi[i,2], NCOL(X)),
               seed = 123, B = B, verbose = TRUE)

    data.frame(id = i, combi[i,], res$Y)

  }, mc.cores = numCores) %>% bind_rows() -> res

  for (k in 1:3) res[,k]  = as.factor(res[,k])

  #creation du bloc de MC
  #recherche des coordonnees des coins du cuboid. La regle est que le cuboid doit etre
  #entierement traverse par chaque front

  min_max = bind_rows(
    res %>% group_by(id) %>% summarise_if(is.numeric, min) %>% mutate(type = "min"),
    res %>% group_by(id) %>% summarise_if(is.numeric, max) %>% mutate(type = "max")
  )

  min_front = min_max %>% dplyr::filter(type == "min") %>% summarise_if(is.numeric, max) %>% unlist()
  max_front = min_max %>% dplyr::filter(type == "max") %>% summarise_if(is.numeric, min) %>% unlist()

  #creation du MC avec une repartition uniforme
  MC = sapply(seq_along(max_front), function(j){
    runif(n = N_MC, min = min_front[j], max = max_front[j])
  })

  #aire du cube de MC
  A = prod(max_front - min_front)

  #calcul du critere (air non dominée par le front)
  parallel::mclapply(seq_len(NROW(combi)), function(i){
    no_dom = no_dominated(front = res %>% dplyr::filter(id == i) %>% select_if(is.numeric),
                          X = MC, sens = sens)

    data.frame(res %>% select(id, crossing_over, mutation) %>%
                 distinct() %>%
                 filter(id == i), aire = A*sum(no_dom)/NROW(MC))
  }, mc.cores = round(numCores/2)) %>% bind_rows() -> aa

  #calcul du volume total
  Y = lapply(fn, function(f) f(X)) %>% bind_cols()
  vol_tot = hypervolume_ellipse(Y, alpha = 0.01)

  #plot des ecart relatifs a la raille du nuage
  A = matrix(aa$aire, ncol = length(aa$aire), nrow = length(aa$aire), byrow = F)
  hist(100*((A - t(A)))[upper.tri(A)]/vol_tot, main = "", xlab = "delta aire front / aire espace objectif")

  #critere en fonction des combinaisons
  print(
    ggplot(aa, aes(x = crossing_over, y = mutation)) +
    geom_point(aes(size = aire, colour = aire)) )

  if (length(fn) < 3){
    #front
    print(
      ggplot(res, aes(x = Y1, y = Y2)) +
      geom_point(aes(colour = id)) +
      geom_line(aes(colour = id)) +
      geom_point(data = Y, aes(x=Y1, y=Y2))
    )
  }
}

#' nsga_tunage
#' @export
nsga_tunage = function(X, fn, sens, N_MC = 10^4){

  # pour paraleliser un peu
  numCores <- parallel::detectCores()

  #definition de la grille de recherche des parametres
  grille = c(
    list(N = rep(seq(10, 50, 10), 2), B = rep(seq(10, 50, 10), 2),
         crossing_over_size = rep( seq_len(round(NCOL(X)/2)), 5)),
         freq_mutation = lapply(seq_along(X), function(i) seq(0, 0.9, 0.1))
  ) %>% as.data.frame()

  for (k in 1:3) grille[,k] = as.factor(grille[,k])

  ff = list(
    performance = function(grille){

      parallel::mclapply(seq_len(NROW(grille)), function(i){
        data.frame(
          NSGA(X = X, fn = fn, N = as.numeric(as.character(grille$N[i])),
               B = as.numeric(as.character(grille$B[i])),
               crossing_over_size = as.numeric(as.character(grille$crossing_over_size[i])),
               freq_mutation = grille[i, 4:NCOL(grille)],
               verbose = F)$Y,
          id = as.factor(i)
        )
      }, mc.cores = numCores) %>% bind_rows() -> all_fronts

      min_max = bind_rows(
        all_fronts %>% group_by(id) %>% summarise_if(is.numeric, min) %>% mutate(type = "min"),
        all_fronts %>% group_by(id) %>% summarise_if(is.numeric, max) %>% mutate(type = "max")
      )

      min_front = min_max %>% dplyr::filter(type == "min") %>% summarise_if(is.numeric, max) %>% unlist()
      max_front = min_max %>% dplyr::filter(type == "max") %>% summarise_if(is.numeric, min) %>% unlist()

      #creation du MC avec une repartition uniforme
      MC = sapply(seq_along(max_front), function(j){
        runif(n = N_MC, min = min_front[j], max = max_front[j])
      })

      #aire du cube de MC
      A = prod(max_front - min_front)

      #calcul du critere (air non dominée par le front)
      parallel::mclapply(seq_len(NROW(grille)), function(i){
        no_dominated(front = all_fronts %>% dplyr::filter(id == i) %>%
                       select_if(is.numeric),
                     X = MC, sens = sens) %>%
          sum() * A / NROW(MC)
      }, mc.cores = round(numCores/2)) %>% unlist()

    },
    time = function(grille){
      parallel::mclapply(seq_len(NROW(grille)), function(i){
        NSGA(X = X, fn = fn, N = as.numeric(as.character(grille$N[i])),
             B = as.numeric(as.character(grille$B[i])),
             crossing_over_size = as.numeric(as.character(grille$crossing_over_size[i])),
             freq_mutation = grille[i, 4:NCOL(grille)],
             verbose = F)$time
      }) %>% unlist()
    }
  )

  NSGA(X = grille, fn = ff, sens = c("min", "min"), N = 10, B = 10, verbose = TRUE)

}
