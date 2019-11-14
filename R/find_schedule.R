#' Restore the optical screening time from the DP result
#' @return \itemize{
#'  \item idx: indeces of the optimal solution
#'  \item val: solution value
#' }
#' @export
find_schedule <- function(dp_res, K){

  if (missing(K)) K = nrow(dp_res$opt)

  sol_idx <- rep(NA, K)

  sol_idx[K] <- which.max(dp_res$opt[K,])
  if (K > 1)
    for (k in (K-1):1)
      sol_idx[k]  <-  dp_res$opt_idx[k, sol_idx[k+1]]

  return(list(idx = sol_idx,
              val = dp_res$time[sol_idx]))
}
