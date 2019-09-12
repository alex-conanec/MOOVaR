rotation_2D <- function(X, teta){
    tt <- matrix(c(cos(teta), -sin(teta),
                   sin(teta), cos(teta)),
                 byrow = TRUE, ncol = 2)
    X %*% t(tt)
}


rotation_3D <- function(X, phi = 0, teta = 0, psi = 0,
                        rot_order = c("phi", "teta", "psi")){

    rot_order[which(rot_order=="phi")] <- "Rx"
    rot_order[which(rot_order=="teta")] <- "Ry"
    rot_order[which(rot_order=="psi")] <- "Rz"

    rotation <- list(
        Rx = matrix(c(1,     0,       0,
                      0, cos(phi), - sin(phi),
                      0, sin(phi),   cos(phi)),
                    byrow = TRUE, ncol = 3),

        Ry = matrix(c(  cos(teta), 0, sin(teta),
                        0,     1,      0,
                        - sin(teta), 0, cos(teta)),
                    byrow = TRUE, ncol = 3),

        Rz = matrix(c(cos(psi), - sin(psi), 0,
                      sin(psi),   cos(psi), 0,
                      0,           0,   1),
                    byrow = TRUE, ncol = 3)
    )



    X %*% t(rotation[[rot_order[1]]]) %*%
        t(rotation[[rot_order[2]]]) %*%
        t(rotation[[rot_order[3]]])
}

