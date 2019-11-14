#' Compute integrals of TP
#'
#' @param trans_prob A list of transition probabilities.
#' @return trans_prob, with int_p01_0 and int_p11 added.
compute_integral_on_tp <- function(trans_prob) {
  t <- trans_prob$time
  L <- length(t)

  t_diff <- t - c(0, t[1:(L-1)])
  # integrate p01(0, t) from 0 to inf
  trans_prob$int_p01_0 <- sum(trans_prob$p01_0 * t_diff)
  # integrate p11: int_p11[i] = int_t[i]^inf p11(t[i], t) dt
  trans_prob$int_p11 <- rep(0, L)
  for (i in 1:(L-1))
    trans_prob$int_p11[i] <- sum(trans_prob$p11[i, (i+1):L] * t_diff[(i+1):L])

  # # compute the probability of T being less than U
  # p_tri <- rep(NA, L) #p_tri[i] = Prob(t[i]-t_gap < T <= t[i], T < D)
  # p_tri[1] <- p01_0[1]
  # for (i in 2:L)
  #   p_tri[i] <- p01_0[i] - p01_0[i-1] * p11[i-1, i]
  # prob_TlessU = sum(p_tri)

  return(trans_prob)
}
