#' Evaluating the (empirical) loss for a given schedule (in idx)
#'
#' @return loss_est: loss evaluated result.
#' @export
eval_loss <- function(sched_idx, tp_est, beta = 1, J = Inf) {
  t <- tp_est$time
  L <- length(t)

  sched_idx <- sort(sched_idx)
  K <- length(sched_idx)

  # compute eta/zeta on sched_idx and evaluate loss.
  if (beta == 1 || J == 0) {
    # eta_a: eta vector on sched_idx. eta_a[i] = eta(a_{i-1}, a_i).
    eta_a <- rep(NA, K)
    eta_a[1] <- with(tp_est,
                     p01_0[sched_idx[1]] * int_p11[sched_idx[1]])
    eta_a[-1] <- with(tp_est,
                      p00_0[sched_idx[-K]] * p01[cbind(sched_idx[-K], sched_idx[-1])] * int_p11[sched_idx[-1]])
    # evaluate loss.
    loss_est <- tp_est$int_p01_0 - beta * sum(eta_a)
  } else {
    # zeta_a: zeta matrix on sched_idx. zeta_a[i, j] = zeta(a_{j-1}, a_j, a_i).
    zeta_a <- matrix(NA, K, K)
    zeta_a[,1] <- with(tp_est,
                       p01_0[sched_idx[1]] * p11[sched_idx[1], sched_idx] * int_p11[sched_idx])
    zeta_a[-1, -1] <- with(tp_est,
                           outer(int_p11[sched_idx[-1]],
                                 p00_0[sched_idx[-K]] * p01[cbind(sched_idx[-K], sched_idx[-1])]) *
                             t(p11[sched_idx[-1], sched_idx[-1]]))
    # evaluate loss.
    loss_est <- tp_est$int_p01_0
    for (j in 0:min(J, K-1)) {
      loss_est <- loss_est - beta * (1 - beta)^j * sum(zeta_a[cbind((1+j):K, 1:(K-j))])
    }
  }

  return(loss_est)
}
