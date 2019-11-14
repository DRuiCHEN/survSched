#' Dynamic Programming procedure
#' Algorithm 1 in the paper
#'
#' @param tp_est List of transition probability estimation.
#' @param maxK Integer.
#' @param output_eta Logical.
#'
#' @return List of \itemize{
#'  \item time, loss_const_term
#'  \item opt: C matrix in the paper
#'  \item opt_idx: argmin_{l} { C[k-1, t[l]] + eta[t[l], t[i]] }
#'  \item const_term: the constant term in the expected loss
#'  \item eta_0, eta (optional): eta function as defined in the paper, which is the partial loss.
#' }
#'
#' @export
#'
dp_proc <- function(tp_est, maxK, output_eta = FALSE) {

  if (is.null(tp_est$int_p01_0)) tp_est <- compute_integral_on_tp(tp_est)

  t <- tp_est$time
  L <- length(t)

  # compute eta: eta_0[i] = eta(0, t[i]), eta[i, j] = eta(t[i], t[j])
  eta_0 <- with(tp_est, p01_0 * int_p11)
  eta <- with(tp_est, outer(p00_0, int_p11) * p01)

  # opt[k, i] is the optimal loss involving the first k variables when a_k takes t[i]
  opt <- matrix(0, maxK, L)
  # opt_idx[k-1, i] = l means when a_{k}=t[i] the minimizer is a_{k-1} = t[l]
  opt_idx <- matrix(NA, maxK-1, L)
  # initialize opt[1,]
  opt[1,] <- eta_0
  # compute opt[k,] (k >= 2) recursively
  if (maxK > 1)
    for (k in 2:maxK)
      for (i in k:L){
        max_idx <- which.max(opt[k-1, (k-1):(i-1)] + eta[(k-1):(i-1), i]) + (k-2)
        opt_idx[k-1, i] <- max_idx
        opt[k, i] <- opt[k-1, max_idx] + eta[max_idx, i]
      }

  res <- list(opt = opt,
              opt_idx = opt_idx,
              time = t,
              const_term = tp_est$int_p01_0)
  if (output_eta) {
    res$eta_0 <- eta_0
    res$eta <- eta
  }
  return(res)
}
