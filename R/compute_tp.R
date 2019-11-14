#' Compute transition probability
#'
#' @param survdat A survTP object
#'
#' @return A list \itemize{
#'  \item p00_0, p01_0: \code{p00_0[i] = p00_0(0, t[i]), p01_0[i] = p01(0, t[i])}
#'  \item p01 and p11: Markov transition probabilites.
#'  \item time: time grid for searching the optimal schedule.
#'  \item int_p01_0, int_p11: integral of p01(0, s)ds and p11(t[i], s)ds
#' }
#'
#' @export
compute_tp <- function(survdat,
                       t_gap = 1/12,
                       method = "AJ",
                       end_time, t) {

  trans_fn <- switch(method,
                     "AJ" = TPmsm::transAJ,
                     "PAJ" = TPmsm::transPAJ
  )

  if (missing(t)) {
    if (missing(end_time)) end_time <- max(survdat[[1]]$Stime)
    t <- seq(0, end_time, by = t_gap)[-1]
    t <- t[t <= max(survdat[[1]]$Stime)]
    if (utils::tail(t, 1) < end_time) t <- c(t, end_time)
  } else {
    end_time <- max(t, survdat[[1]]$Stime)
    t <- t[t <= max(survdat[[1]]$Stime)]
    if (utils::tail(t, 1) < end_time) t <- c(t, end_time)
  }
  L <- length(t)

  # p00_0, p01_0, int_p01_0
  p00_0 <- p01_0 <- rep(NA, L)
  p_from0 <- trans_fn(object=survdat,
                      s = 0,
                      conf = FALSE,
                      state.names = c('0','1','2'))
  t_diff <- with(p_from0, c(time[-1], end_time) - time)
  int_p01_0 <- sum(p_from0$est[, '0 1'] * t_diff)
  i <- k <- 1
  while (i <= L) {
    if (k <= length(p_from0$time) && p_from0$time[k] <= t[i])
      k <- k + 1
    else {
      p00_0[i] <- p_from0$est[k-1, '0 0']
      p01_0[i] <- p_from0$est[k-1, '0 1']
      i <- i + 1
    }
  }

  # p01, p11, int_p11
  p01 <- p11 <- matrix(NA, L, L)
  int_p11 <- rep(0, L)
  diag(p01) <- 0
  diag(p11) <- 1
  for (i in 1:(L-1)){
    p_fromi <- trans_fn(object=survdat,
                        s = t[i],
                        conf = FALSE,
                        state.names = c('0','1','2'))
    t_diff <- with(p_fromi, c(time[-1], end_time) - time)
    int_p11[i] <- sum(p_fromi$est[, '1 1'] * t_diff)
    j <- i + 1
    k <- 1
    while (j <= L) {
      if (k <= length(p_fromi$time) && p_fromi$time[k] <= t[j])
        k <- k + 1
      else {
        p01[i, j] <- p_fromi$est[k-1, '0 1']
        p11[i, j] <- p_fromi$est[k-1, '1 1']
        j <- j + 1
      }
    }
  }

  return(list(time = t,
              p00_0 = p00_0,
              p01_0 = p01_0,
              p01 = p01,
              p11 = p11,
              int_p01_0 = int_p01_0,
              int_p11 = int_p11))
}
