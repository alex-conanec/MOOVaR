#' Crossing over operator
#'
#' This function switch components for the given variables position between two
#' individus (parent) which have been selected after the tounament
#'
#' @param couple a matrix/data.frame of two individus
#' @param crossing_position a sorted vector of variable position to switch
#'
#' @return a matrix/data.frame of the two childs from the crossing over of the
#' parents
#'
#' @examples
#' sum(1:10)
#' @export

# crossing_over <- function(X, parents, crossing_over_size){
crossing_over_uniforme <- function(parents, n_echange, param_mut){

    d = NCOL(parents)
    n = NROW(parents)

    if (n%%2 == 1){
      n = n-1
      warning("the number of parent is odd, therefore only n-1 children were created")
    }

    res = lapply(seq(1, n, 2), function(i){
      c1 = sample(seq_len(d), size = n_echange, replace = FALSE)
      c2 = seq_len(NCOL(parents))[-c1]

      list(
        X_C = bind_rows(
          data.frame(parents[i, c1, drop = F], parents[i+1, c2, drop = F]),
          data.frame(parents[i, c2, drop = F], parents[i+1, c1, drop = F])
        ),
        param_mut = list(
          ifelse(length(c1) >= length(c2), param_mut[i], param_mut[i+1]),
          ifelse(length(c1) >= length(c2), param_mut[i+1], param_mut[i])
        )
      )

    })


    param_mut_C = lapply(res, function(x) c(x$param_mut[[1]], x$param_mut[[2]]))
    param_mut_C = lapply(0:(n-1), function(i){
      param_mut_C[[(i%/%2)+1]][[i%%2+1]]
    })

    list(
      X_C = lapply(res, function(x) x$X_C) %>% bind_rows() %>%
        select(!!colnames(parents)),
      param_mut_C = param_mut_C
    )


}

#' @export
crossing_over_bloc <- function(parents, param_mut, n_bloc = 2){

  d = NCOL(parents)
  n = NROW(parents)

  if (n%%2 == 1){
    n = n-1
    warning("the number of parent is odd, therefore only n-1 children were created")
  }

  res = lapply(seq(1, n, 2), function(i){
    var_idx = c(0, sort(sample(2:(d-1), n_bloc)), d)
    c1 = sapply(seq(1, n_bloc+1, 2), function(j=5){
      (var_idx[j]+1):var_idx[j+1]
    }) %>% unlist()

    c2 = seq_len(d)[-c1]

    list(
      X_C = bind_rows(
        data.frame(parents[i, c1, drop = F], parents[i+1, c2, drop = F]),
        data.frame(parents[i, c2, drop = F], parents[i+1, c1, drop = F])
      ),
      param_mut = list(
        ifelse(length(c1) >= length(c2), param_mut[i], param_mut[i+1]),
        ifelse(length(c1) >= length(c2), param_mut[i+1], param_mut[i])
      )
    )

  })

  param_mut_C = lapply(res, function(x) c(x$param_mut[[1]], x$param_mut[[2]]))
  param_mut_C = lapply(0:(n-1), function(i){
    param_mut_C[[(i%/%2)+1]][[i%%2+1]]
  })

  list(
    X_C = lapply(res, function(x) x$X_C) %>% bind_rows() %>%
      select(!!colnames(parents)),
    param_mut_C = param_mut_C
  )

}


#' Mutation operator
#'
#' This function assign random value to the individu caracteristique following
#' a frequency of occurence
#'
#' @param Qt a matrix/data.frame of all the individus generated after crossing
#' over
#' @param freq_m a vector of numerics value varying from 0 to 1 which represent
#' the frequency of mutation of each variable
#'
#' @param distribution choice of the distribution where for the mutation value
#'
#' @return a matrix/data.frame of the two childs from the crossing over of the
#' parents
#'
#' @examples
#' sum(1:10)
#' @export
mutation_simple <- function(X, freq_m, distri_Xi, type_var){

  n = NROW(X)
  d = NCOL(X)

  mut_decision = matrix(sample(c(FALSE, TRUE), n*d, prob = c(1-freq_m, freq_m),
                               replace = TRUE), ncol = d)


  X_replacement = lapply(1:d, function(i){
    if (length(distri_Xi[[i]]) > 1){
      if (type_var[i] == "numeric"){
        runif(n, distri_Xi[[i]]$min, distri_Xi[[i]]$max)
      }else{
        round(runif(n, distri_Xi[[i]]$min, distri_Xi[[i]]$max))
      }
    }else{
      as.factor(sample(distri_Xi[[i]]$levels, n, replace = TRUE))
    }
  }) %>% as.data.frame()

  for (j in 1:d){
    X[mut_decision[,j], j] = X_replacement[mut_decision[,j],j]
  }

  X
}


#' @export
mutation_mixte <- function(X, freq_m, param_mut, distri_Xi,
                           type_var,
                           beta = 5,
                           epsilon0 = 10^(-4),
                           tau_1 = 1/sqrt(2*sqrt(sum(type_var != "factor"))),
                           tau_2 = 1/sqrt(2*sum(type_var != "factor")) ){

  n = NROW(X)
  is_fac = type_var == "factor"
  d1 = sum(!is_fac)
  d2 = sum(is_fac)

  # borne_quanti = lapply(which(!is_fac), function(i) distri_Xi[[i]])
  levels_quali = lapply(which(is_fac), function(i) distri_Xi[[i]]$levels)

  #quanti
  if (any(c("numeric", "integer") %in% type_var)){
    res_quanti = lapply(param_mut, function(x){
      C = matrix(0, nrow = d1, ncol = d1)
      diag(C) = x$sigma
      for (i in 1:d1){
        for (j in (1:d1)[-i]){
          if (i < j){
            C[i,j] = abs(1/2*(x$sigma[j]^2 - x$sigma[i]^2))*tan(2*x$alpha[i,j])
          }else{
            C[i,j] = C[j,i]
          }

        }
      }

      list(
        delta = purrr::quietly(mvtnorm::rmvnorm)(n=1, mean=rep(0, d1),
                                                 sigma = C, method="chol")$result %>%
          as.data.frame(),

        param_mut = list(
          sigma = x$sigma*exp(rnorm(d1, 0, tau_1) + rnorm(d1, 0, tau_2)) %>%
            sapply(function(x) max(epsilon0, x)),
          alpha = x$alpha + matrix(beta*rnorm(d1^2), d1, d1)
        )

      )

    })

    delta = lapply(res_quanti, function(x) x$delta) %>% bind_rows()
    param_mut = lapply(res_quanti, function(x) x$param_mut)

    X[,!is_fac] = X[,!is_fac] + delta
    X[,type_var=="integer"] = round(X[,type_var=="integer"])

    for (j in which(!is_fac)){
      X[X[,j] > distri_Xi[[j]]$max ,j] = distri_Xi[[j]]$max
      X[X[,j] < distri_Xi[[j]]$min ,j] = distri_Xi[[j]]$min
    }

  }

  #quali
  if (any(type_var == "factor")){
    mut_decision = matrix(sample(c(FALSE, TRUE), n*d2, prob = c(1-freq_m, freq_m),
                                 replace = TRUE), ncol = d2)
    X_replacement = sapply(levels_quali, function(levs){
      sample(levs, n, replace = TRUE)
    }) %>% as.data.frame()

    for (j in 1:d2){
      X[mut_decision[,j], which(is_fac)[j]] = X_replacement[mut_decision[,j],j]
    }
  }

  list(X = X, param_mut = param_mut)

}


# n = 30
# d1 = 3
# d2 = 2
# levs = sample(2:4, d2, replace = TRUE)
# X_quanti = as.data.frame(matrix(runif(n, -10, 10), nrow = n, ncol = d1))
# X_quanti[,1] = round(X_quanti[,1])
# colnames(X_quanti) = paste0("X", 1:d1)
#
# X_quali = lapply(1:d2, function(j){
#   factor(sample(1:levs[j], n, replace = TRUE))
# }) %>% bind_cols()
# colnames(X_quali) = paste0("X", (d1+1):(d1+d2))
#
# X = bind_cols(X_quanti, X_quali)
#
# levels_quali = lapply(levs, function(lev){
#   factor(1:lev)
# })
# mutation_mixte(X, freq_m = 0.2,
#                sigma = sapply(X[,1:d1], sd)*1,
#                levels_quali = levels_quali,
#                borne_quanti = lapply(1:d1, function(x) list(min = -10, max = 10)),
#                type_var = c("integer", "numeric", "numeric", "factor", "factor"))

