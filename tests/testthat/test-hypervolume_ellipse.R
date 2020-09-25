test_that("test that volume calculated by the hypervolume formulas is equal to the MC method", {

# p = 2
n = 300
Y1 = rnorm(n)
Y2 = Y1 + rnorm(n)
Y = data.frame(Y1, Y2)

ellipse = ellipse_var(Y1=Y$Y1, Y2=Y$Y2, alpha = 0.01)
# plot(ellipse, type='l')
# points(Y, col = "grey")

B = 10^4
MC = data.frame(runif(B, -3, 3), runif(B, -4, 4))
A = 6*8

all_sens = matrix(unlist(purrr::cross2(c("max", "min"), c("max", "min"))),
                  ncol = p, byrow = T)

lapply(seq_len(NROW(all_sens)), function(i){
  front_ellipse = ellipse[dominance_ranking(ellipse, sens = all_sens[i,]) == 1, ]
  !no_dominated(front = front_ellipse, X = MC, sens = all_sens[i,])
}) %>% bind_cols() %>% apply(1, all) -> res


plot(MC, col = c("red", "black")[as.factor(res)])
points(ellipse, type = 'l', col = 'blue', lwd = 3)


MC_method = A*sum(res)/NROW(MC)
hypervolume_method = hypervolume_ellipse(Y, alpha = 0.01)

expect_that(abs((hypervolume_method - MC_method))/hypervolume_method,
            is_less_than(0.1))

expect_that(abs((hypervolume_method - MC_method))/MC_method,
            is_less_than(0.1))

})
