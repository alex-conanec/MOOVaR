#' @export
multi_front <- function(NSGA_param, model, type = "bootstrap", B = 100){

    require(parallel)
    require(magrittr)
    numCores <- detectCores()
    n <- NROW(model$X_train)

    if (type == "bootstrap"){
        model_list <- lapply(seq_len(B), function(i){
            train <- sample(x = seq_len(n), size = n, replace = TRUE)
            lapply(seq_along(model$param), function(j){
                generate_model(type = model$type,
                               X_train = model$X_train[train, ],
                               Y_train = model$Y_train[train,j],
                               param = model$param[[j]])
            })
        })
    }else if (type == "var_param"){

        formules <- sapply(model$param, function(param) param$formule)
        model_list <- var_param(formules,
                                X_train = model$X_train,
                                Y_train = model$Y_train,
                                B = B)
    }


    mclapply(model_list, function(model){
        for (j in seq_along(model)){
            NSGA_param$ff[[j]]$f <- model[[j]]
        }
        do.call(NSGA, NSGA_param)
    },  mc.cores = numCores)
}


#' @export
generate_model <- function(type = model$type, X_train = X_train,
                           Y_train = Y_train, param = model$param[[j]]){
    switch (type,
        "lm" = {
            lm(param$formule, data = cbind(Y = Y_train, X_train))
        },
        "rf" = {
            randomForest::randomForest(x = X_train, y = Y_train,
                                       mtry = param$mtry)
        }
    )
}

#' @export
var_param <- function(formules, X_train, Y_train, B = 50){

    reg_param <- lapply(seq_along(formules), function(j){
        reg <- lm(formules[[j]], data = cbind(Y = Y_train[,j], X_train))

        # deal with formulas and the matrix
        var_reg <- formules[[j]] %>%
            stringr::str_extract(pattern = "(?<=~).*") %>%
            strsplit(split = "+", fixed = T) %>% unlist() %>%
            strsplit(split = ":")

        var_reg <- lapply(var_reg, function(var){
            if (length(var) == 1){
                without_space <- strsplit(var, split = " ")[[1]]
                without_space[! without_space %in% ""]
            }else{
                sapply(var, function(v){
                    without_space <- strsplit(v, split = " ")[[1]]
                    without_space[! without_space %in% ""]
                })
            }
        })

        X_reg <- sapply(var_reg, function(var){
            if (length(var) == 1){
                X_train[, var]
            }else{
                sapply(var, function(v){
                    X_train[, v]
                }) %>% apply(1, prod)
            }
        })

        X_reg <- cbind(1, X_reg)


        #sce & mse
        sce <- (Y_train[,j] - reg$fitted.values)^2 %>% sum()
        p <- NCOL(X_reg); n <- NROW(X_reg)
        mse <- (1/(n - p - 1)) * sce

        # loi MVN
        beta <- mvtnorm::rmvnorm(n = B,
                         mean = reg$coefficients,
                         sigma = mse * solve(t(X_reg)%*%X_reg))

        list(beta = beta, var_reg = var_reg)
    })

    lapply(seq_len(B), function(i){
        lapply(seq_along(formules), function(j){
            model <- function(X, beta, var_reg){

                X_reg <- sapply(var_reg, function(var){
                    if (length(var) == 1){
                        X[, var]
                    }else{
                        sapply(var, function(v){
                            X[, v]
                        }) %>% apply(1, prod)
                    }
                })

                X_reg <- cbind(1, X_reg)

                as.matrix(X_reg) %*% beta
            }
            formals(model)$beta <- reg_param[[j]]$beta[i,]
            formals(model)$var_reg <- reg_param[[j]]$var_reg
            model
        })
    })
}




