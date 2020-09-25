context("pairewise_comparaison")
library(dplyr)

test_that("la longueur du resultat est egale au nombre de point", {
    p <- 5
    X <- matrix(runif(sample(100, size = 1) * p), ncol = p)
    expect_equal(pairwise_comparison(X, sens = rep("min", p)) %>% length(),
                 NROW(X))
})

test_that("la somme des domines est egale a la somme des dominants",{
    p <- 5
    X <- matrix(runif(sample(100, size = 1) * p), ncol = p)

    comp <- pairwise_comparison(X, sens = rep("min", p))

    sum_dominated <- sapply(comp, function(point){
        point$dominated_count
    }) %>% sum()

    sum_dominating <- sapply(comp, function(point){
        point$dominating_index %>% length()
    }) %>% sum()

    expect_equal(sum_dominated, sum_dominating)
})

test_that("un point ne peut etre domine par un point qu'il domine", {
    p <- 5
    X <- matrix(runif(sample(100, size = 1) * p), ncol = p)

    comp <- pairwise_comparison(X, sens = rep("min", p))

    for (i in seq_along(comp)){
      for (index in comp[[i]]$dominating_index){
        expect_false(i %in% comp[[index]]$dominating_index)
      }
    }

})

test_that("pairwise_comparaison behave well with two type examples of 5 points", {

  #double max
  p <- 2
  X <- matrix(c(0, 0,
                1, 1,
                1, 2,
                0, 2,
                -1, 0.1), ncol = p, byrow = T)

  comp <- pairwise_comparison(X, sens = rep("max", p))

  true_comp = list(
    list(
      dominating_index = integer(0),
      dominated_count = 3 #c(2, 3, 4)
    ),
    list(
      dominating_index = c(1, 5),
      dominated_count = 1 #3
    ),
    list(
      dominating_index = c(1, 2, 4, 5),
      dominated_count = 0
    ),
    list(
      dominating_index = c(1, 5),
      dominated_count = 1 #3
    ),
    list(
      dominating_index = integer(0),
      dominated_count = 3 #c(2, 3, 4)
    )
  )

  for (i in length(true_comp)){
    expect_equal(true_comp[[i]]$dominated_count, comp[[i]]$dominated_count)
    expect_equal(true_comp[[i]]$dominated_index, comp[[i]]$dominated_index)
  }

  #mixte max/min
  comp <- pairwise_comparison(X, sens = c("max", "min"))

  true_comp = list(
    list(
      dominating_index = 4:5,
      dominated_count = 0
    ),
    list(
      dominating_index = c(3, 4),
      dominated_count = 0
    ),
    list(
      dominating_index = 4,
      dominated_count = 2 #c(1, 2)
    ),
    list(
      dominating_index = integer(0),
      dominated_count = 3 #1:3
    ),
    list(
      dominating_index = integer(0),
      dominated_count = 1 #1
    )
  )

  for (i in length(true_comp)){
    expect_equal(true_comp[[i]]$dominated_count, comp[[i]]$dominated_count)
    expect_equal(true_comp[[i]]$dominated_index, comp[[i]]$dominated_index)
  }

})
