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

        betas <- lapply(seq_along(model$param), function(j){
            reg <- lm(model$param[[j]]$formule,
                      data = cbind(Y = model$Y_train[,j], model$X_train))

            #sce & mse
            sce <- (model$Y_train[,j] - reg$fitted.values)^2 %>% sum()
            p <- NCOL(model$X_train)
            mse <- (1/(n - p - 1)) * sce

            # deal with formulas and the matrix
            var_reg <- model$param[[j]]$formule %>%
                stringr::str_extract(pattern = "(?<=~).*") %>%
                strsplit(split = "+", fixed = T) %>% unlist() %>%
                strsplit(split = ":")

            X_reg <- sapply(var_reg, function(var){
                if (length(var) == 1){
                    without_space <- strsplit(var, split = " ")[[1]]
                    X[, without_space[! without_space %in% ""]]
                }else{
                    sapply(var, function(v){
                        without_space <- strsplit(v, split = " ")[[1]]
                        X[, without_space[! without_space %in% ""]]
                    }) %>% apply(1, prod)
                }
            })

            X_reg <- cbind(1, X_reg)

            # loi MVN
            mvtnorm::rmvnorm(n = B,
                             mean = reg$coefficients,
                             sigma = mse * solve(t(X_reg)%*%X_reg))

        })

        model_list <- lapply(seq_len(B), function(i){
            lapply(seq_along(model$param), function(j){
                model <- function(X, beta){
                    as.matrix(X) %*% beta
                }
                formals(model)$beta <- betas[[j]][i,]
                model
            })
        })
    }

    model=model_list[[1]]
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
