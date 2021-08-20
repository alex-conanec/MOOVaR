#' @export
is_it_stopped = function(path_tracking){
  if (!is.null(path_tracking)){
    tracking = readRDS(path_tracking)
    if (tracking[[length(tracking)]] == "stop"){
      stop("stop algo")
    }
  }
}

#' @export
tracking_msg = function(path_tracking, msg = ""){
  if (!is.null(path_tracking)){
    tracking = readRDS(path_tracking)
    tracking[[length(tracking)+1]] = msg
    saveRDS(tracking, path_tracking)
  }
}
