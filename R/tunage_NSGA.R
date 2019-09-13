#' @export
tunage_NSGA <- function(boundary_param, boundary_space, NSGA_param, nMC = 10^3,
                        n = 32, B = 20){

    #on cherche le min/max de chaque function a optimiser
    x0 = unlist(NSGA_param$X[1,])
    boundary_MC_space <- lapply(seq_along(NSGA_param$ff), function(j){
        list(
            min = optim(par = x0, fn = NSGA_param$ff[[j]]$f,  #marchera pas si c'est un model
                        lower = boundary_space$lower_bound,
                        upper = boundary_space$upper_bound,
                        method = "L-BFGS-B",
                        control = list(fnscale = 1))$value,
            max = optim(par = x0, fn = NSGA_param$ff[[j]]$f,
                        lower = boundary_space$lower_bound,
                        upper = boundary_space$upper_bound,
                        method = "L-BFGS-B",
                        control = list(fnscale = -1))$value
        )
    })

    #on creer un cube de MC pour tout le monde
    spaceMC <- sapply(seq_along(boundary_MC_space), function(j){
        runif(n=nMC, min = boundary_MC_space[[j]]$min,
              max = boundary_MC_space[[j]]$max)
    })

    #on creer une function qui retourne l'air sous la courbe dans le cube MC
    accuracy <- function(x, NSGA_param, spaceMC){

        if (class(x) %in% c("data.frame", "matrix")){
            numCores <- parallel::detectCores()
            parallel::mclapply(seq_len(NROW(x)), function(i){
                if (is.null(NSGA_param$N)) N <- round(NROW(NSGA_param$X) / 2, 0)
                NSGA_param$B <- x[i,"B"]
                NSGA_param$crossing_over_size = x[i,"crossing_over_size"]
                NSGA_param$freq_mutation = rep(x[i,"freq_mutation"],
                                               NCOL(NSGA_param$X))

                front <- do.call(NSGA, NSGA_param)

                sens <- sapply(NSGA_param$ff, function(f){
                    f$sens
                })

                count_under_front(front = front$Y,
                                  X = spaceMC,
                                  sens = sens)
            }, mc.cores = numCores) %>% unlist()
        }else if (class(x) %in% c("factor", "numeric")){
            if (is.null(NSGA_param$N)) N <- round(NROW(NSGA_param$X) / 2, 0)
            NSGA_param$B <- x[,"B"]
            NSGA_param$crossing_over_size = x[,"crossing_over_size"]
            NSGA_param$freq_mutation = rep(x[,"freq_mutation"], NCOL(NSGA_param$X))

            front <- do.call(NSGA, NSGA_param)

            sens <- sapply(NSGA_param$ff, function(f){
                f$sens
            })

            count_under_front(front = front$Y,
                              X = spaceMC,
                              sens = sens)
        }
    }
    formals(accuracy)$NSGA_param <- NSGA_param
    formals(accuracy)$spaceMC <- spaceMC


    time_running <- function(x, NSGA_param){

        if (class(x) %in% c("data.frame", "matrix")){
            numCores <- parallel::detectCores()
            parallel::mclapply(seq_len(NROW(x)), function(i){
                if (is.null(NSGA_param$N)) N <- round(NROW(NSGA_param$X) / 2, 0)
                NSGA_param$B <- x[i,"B"]
                NSGA_param$crossing_over_size = x[i,"crossing_over_size"]
                NSGA_param$freq_mutation = rep(x[i,"freq_mutation"],
                                               NCOL(NSGA_param$X))

                system.time(do.call(NSGA, NSGA_param))[3]

            }, mc.cores = numCores) %>% unlist()
        }else if (class(x) %in% c("factor", "numeric")){
            if (is.null(NSGA_param$N)) N <- round(NROW(NSGA_param$X) / 2, 0)
            NSGA_param$B <- x[,"B"]
            NSGA_param$crossing_over_size = x[,"crossing_over_size"]
            NSGA_param$freq_mutation = rep(x[,"freq_mutation"], NCOL(NSGA_param$X))

            system.time(do.call(NSGA, NSGA_param))[3]
        }
    }
    formals(time_running)$NSGA_param <- NSGA_param


    #minimisation of to_opt
    B_bound <- c(boundary_param$lower_bound$B,
                 boundary_param$upper_bound$B)

    crossOvSize <- c(boundary_param$lower_bound$crossing_over_size,
                     boundary_param$upper_bound$crossing_over_size)

    param0 <- data.frame(
        B = factor(sample(
            x = seq(B_bound[1], B_bound[2]),
            size = n, replace = TRUE),
            levels = seq(B_bound[1], B_bound[2])
        ),
        crossing_over_size = factor(sample(
            x = seq(crossOvSize[1], crossOvSize[2]),
            size = n, replace = TRUE),
            levels = seq(crossOvSize[1], crossOvSize[2])
        ),
        freq_mutation = seq(
            boundary_param$lower_bound$freq_mutation,
            boundary_param$upper_bound$freq_mutation,
            by = (
                boundary_param$upper_bound$freq_mutation - boundary_param$lower_bound$freq_mutation
            ) / (n - 1)
        )
    )

    NSGA(X = param0,
         ff = list(accuracy = list(f = accuracy, sens = "min"),
                   time_running = list(f = time_running, sens = "min")),
         crossing_over_size = 2,
         freq_mutation = c(0.5, 0.5, 0.5),
         B = B)

}


# #' @export
# tunage_NSGA <- function(X, ff, N, seed = 123, constraints = NULL, n_MC = 10^4,
#                         n_comp = 32, B_range = 10:1000, crossing_type = "hybrid",
#                         verbose = TRUE, time_out = 5){
#
#
#     #amelioration possible, X big et tunage de N
#     if (!missing(seed)){
#         set.seed(seed)
#     }
#
#     if (! crossing_type %in% c("hybrid", "auto"))
#         stop("crossing_type must be either 'hybrid' or 'auto'")
#
#     if (missing(N)) N <- round(NROW(X) / 2, 0)
#
#     sens <- sapply(ff, function(f){
#         f$sens
#     })
#
#     numCores <- parallel::detectCores()
#
#     set.seed(seed)
#     candidates <- parallel::mclapply(seq_len(n_comp), function(i){
#         param <- list(crossing_over_size = sample(x = seq_along(X), size = 1,
#                                                   replace = TRUE),
#                     freq_mutation = rep(x = runif(1), times = NCOL(X)),
#                     B = sample(x = B_range, size = 1, replace = TRUE))
#
#         list(front =
#                  NSGA(X = X, ff = ff, N = N,
#                       crossing_over_size = param$crossing_over_size,
#                       freq_mutation = param$freq_mutation,
#                       B = param$B,
#                       seed = seed,
#                       constraints = constraints, crossing_type = crossing_type,
#                       verbose = verbose, time_out = time_out)$Y,
#              param = param)
#     }, mc.cores = numCores)
#
#     res <- candidates[[1]]
#     gain <- NULL
#     for (i in seq_along(candidates)[-1]){
#
#         diff <- area_btw_front(front_sup = res$front,
#                                front_inf = candidates[[i]]$front,
#                                sens = sens,
#                                n = n_MC)
#
#         if (diff >= 0){
#             if (FALSE){
#                 break
#             }
#         }else{
#             res <- candidates[[i]]
#             gain <- c(gain, diff)
#             names(gain)[length(gain)] <- i
#         }
#     }
#     candidate_param <- lapply(candidates, function(candidate){
#         candidate$param
#     })
#     list(res = res$param, gain = gain, candidates = candidate_param)
# }
