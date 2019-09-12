context("area_btw_front")
# library(dplyr)
# library(testthat)

test_that("test que l'aire d'un parallelograme egale 3", {

    # parallelogramme
    x <- seq(0, 3, 0.01)
    front1 <- cbind(x = x, y = x*1/3)
    front2 <- cbind(x = x, y = 1 + x*1/3)

    area_esti <- area_btw_front(front_sup = front1,
                                front_inf = front2,
                                sens = c("max", "min"),
                                n = 10^4)

    expect_true(area_esti >= 3*0.985 & area_esti <= 3*1.015) #estimation +/- 1.5%
})

test_that("test que le volume d'un cube d'arrete 2 egale 8", {

    x <- runif(200, 0, 2)
    y <- runif(200, 0, 2)
    z <- runif(200, 0, 2)

    front1 <- rbind(cbind(x = x, y = y, z=0),
                    cbind(x = x, y = 0, z=z),
                    cbind(x = 2, y = y, z=z)
    )

    front2 <- rbind(cbind(x = x, y = y, z=2),
                    cbind(x = x, y = 2, z=z),
                    cbind(x = 0, y = y, z=z)
    )

    front1_rot <- rotation_3D(X = front1, phi = pi/20, psi = pi/20, teta = pi/20)
    front2_rot <- rotation_3D(X = front2, phi = pi/20, psi = pi/20, teta = pi/20)

    # Show the cubes in a plot space
    # plot3D::scatter3D(x = front1_rot[,1], y = front1_rot[,2], z = front1_rot[,3])
    # plot3Drgl::plotrgl()
    # plot3D::scatter3D(x = front2_rot[,1], y = front2_rot[,2], z = front2_rot[,3])
    # plot3Drgl::plotrgl()

    area_esti <- area_btw_front(front_sup = front1_rot,
                                         front_inf = front2_rot,
                                         sens = c("max", "max", "min"),
                                         n = 10^4)

    expect_true(area_esti >= 8*0.90 & area_esti <= 8*1.1) #estimation +/- 10%

})


