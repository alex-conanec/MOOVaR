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

    lapply(seq_along(comp), function(i){
        lapply(comp[[i]]$dominating_index, function(index){
            expect_false(i %in% comp[[index]]$dominating_index)
        })
    })
})
