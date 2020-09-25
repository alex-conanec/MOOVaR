test_that("pairwise_comparaison behave well with two type examples of 5 points", {

  p <- 2
  X <- matrix(c(0, 0,
                1, 1,
                1, 2,
                0, 2,
                -1, 0.1), ncol = p, byrow = T)

  #double max
  dom = dominance_ranking(X, sens = c("max", "max"))
  true_dom = c(3,2,1,2,3)
  expect_equal(dom, true_dom)

  #mixte max/min
  dom = dominance_ranking(X, sens = c("max", "min"))
  true_dom = c(1,1,2,3,2)
  expect_equal(dom, true_dom)
})
