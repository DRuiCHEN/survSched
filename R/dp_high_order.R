#' Dynamic Programming with J = 1 or 2 (under Markov assumption). dy_proc + find_schedule.
#'
#' @return \itemize{
#'  \item idx: indeces of the optimal solution
#'  \item val: solution value
#'  \item loss_est: estimated expected loss (with corresponding J)
#' }
#' @export
dp_high_order <- function(tp_est, maxK, beta = 1, J = 1,
                          echo = TRUE, only_maxK = TRUE) {

  if (is.null(tp_est$int_p01_0)) tp_est <- compute_integral_on_tp(tp_est)

  if (maxK == 1 || beta == 1 || J == 0)
    return(find_schedule(dp_proc(tp_est, maxK)))

  t <- tp_est$time
  L <- length(t)

  ## compute zeta ##
  # zeta_0[i, j] = zeta(0, t[i], t[j]), zeta[i, j, l] = zeta(t[i], t[j], t[l]).
  # (eta(v, w) = zeta(v, w, w) is also stored in the zeta array)
  zeta_0 <- with(tp_est, outer(p01_0, int_p11) * p11)
  zeta <- array(NA, dim = c(L, L, L))
  for (j in 1:L) zeta[,j,] <- with(tp_est, outer(p00_0 * p01[,j], p11[j,] * int_p11))

  if (echo) cat("Computation of zeta finished.\n")

  ## dp_proc ##
  if (J == 1) {
    if (maxK < 2) stop("maxK should be at least 2 for J = 1.")
    # opt[k, i, j] is the optimal loss involving the first k visits when a_{k-1}=t[i] and a_k=t[j], i < j
    opt <- array(0, dim = c(maxK, L, L))
    # opt_idx[k-2, i, j] = z means when when a_{k-1}=t[i] and a_k=t[j] the minimizer is a_{k-2} = t[z]
    opt_idx <- array(NA, dim = c(maxK-2, L, L))

    # skip opt[1,,], initialize opt[2,,]
    for (j in 2:L)
      for (i in 1:(j-1))
        opt[2, i, j] <- zeta_0[i, i] + zeta[i, j, j] + (1 - beta) * zeta_0[i, j]

    # compute opt[k,,] (k >= 3) recursively
    if (maxK > 2)
      for (k in 3:maxK) {
        if (echo) cat("Computing opt: k = ", k, "\n")
        for (j in k:L)
          for (i in (k-1):(j-1)) {
            max_idx <- which.max(
              opt[k-1, (k-2):(i-1), i] + (1 - beta) * zeta[(k-2):(i-1), i, j]
            )[1] + (k-3)
            opt_idx[k-2, i, j] <- max_idx
            opt[k, i, j] <- opt[k-1, max_idx, i] +
              (1 - beta) * zeta[max_idx, i, j] +
              zeta[i, j, j]
          }
      }

  } else if (J == 2) {
    if (maxK < 3) stop("maxK should be at least 3 for J = 2.")
    # opt[k, i, j, l] is the optimal loss involving the first k visits when a_{k-1}=t[i] and a_k=t[j], i < j < l
    opt <- array(0, dim = c(maxK, L, L, L))
    # opt_idx[k-2, i, j] = l means when when a_{k-2}t[i], a_{k-1}=t[j] and a_k=t[l] the minimizer is a_{k-3} = t[l]
    opt_idx <- array(NA, dim = c(maxK-3, L, L, L))

    # skip opt[1,,,] and opt[2,,,], initialize opt[3,,,]
    for (l in 3:L)
      for (j in 2:(l-1))
        for (i in 1:(j-1)) {
          opt[3, i, j, l] <- zeta_0[i, i] + zeta[i, j, j] + zeta[j, l, l] +
            (1 - beta) * zeta_0[i, j] + (1 - beta) * zeta[i, j, l] +
            (1 - beta)^2 * zeta_0[i, l]
        }

    # compute opt[k,,,] (k >= 3) recursively
    if (maxK > 3)
      for (k in 4:maxK) {
        if (echo) cat("Computing opt: k = ", k, "\n")
        for (l in k:L)
          for (j in (k-1):(l-1))
            for (i in (k-2):(j-1)) {
              max_idx <- which.max(
                opt[k-1, (k-3):(i-1), i, j] + (1 - beta)^2 * zeta[(k-3):(i-1), i, l]
              )[1] + (k-4)
              opt_idx[k-3, i, j, l] <- max_idx
              opt[k, i, j, l] <- opt[k-1, max_idx, i, j] +
                (1 - beta)^2 * zeta[max_idx, i, l] +
                (1 - beta) * zeta[i, j, l] +
                zeta[j, l, l]
            }
      }

  } else stop("Only support J = 0, 1 or 2!")

  if(echo) cat("Computation of opt finished.\n")

  ## find_schedule ##
  res <- list()
  if (only_maxK) Ks <- maxK else Ks <- (J+1):maxK

  for (K in Ks) {
    sol_idx <- rep(NA, K)
    if (echo) cat("Find optimal schedule with", K, "visits.\n")
    if (J == 1) {
      # trace back opt to find the optimal schedule
      sol_idx[(K-1):K] <- which(opt[K,,] == max(opt[K,,]), arr.ind = TRUE)[1,]
      if (K > 2) {
        for (k in (K-2):1) sol_idx[k]  <-  opt_idx[k, sol_idx[k+1], sol_idx[k+2]]
      }
      # evaluate loss
      loss_est <- tp_est$int_p01_0 -
        beta * zeta_0[sol_idx[1], sol_idx[1]] -
        beta*(1-beta) * zeta_0[sol_idx[1], sol_idx[2]] -
        beta * sum(zeta[cbind(sol_idx[-K], sol_idx[-1], sol_idx[-1])]) -
        beta*(1-beta) * ifelse(K == 2, 0, sum(zeta[cbind(sol_idx[1:(K-2)], sol_idx[2:(K-1)], sol_idx[3:K])]))

    } else if (J == 2) {
      # trace back opt to find the optimal schedule
      sol_idx[(K-2):K] <- which(opt[K,,,] == max(opt[K,,,]), arr.ind = TRUE)[1,]
      if (K > 3) {
        for (k in (K-3):1) sol_idx[k] <- opt_idx[k, sol_idx[k+1], sol_idx[k+2], sol_idx[k+3]]
      }
      # evaluate loss
      loss_est <- tp_est$int_p01_0 -
        beta * zeta_0[sol_idx[1], sol_idx[1]] -
        beta*(1-beta) * zeta_0[sol_idx[1], sol_idx[2]] -
        beta*(1-beta)^2 * zeta_0[sol_idx[1], sol_idx[3]] -
        beta * sum(zeta[cbind(sol_idx[-K], sol_idx[-1], sol_idx[-1])]) -
        beta*(1-beta) * sum(zeta[cbind(sol_idx[1:(K-2)], sol_idx[2:(K-1)], sol_idx[3:K])]) -
        beta*(1-beta)^2 * ifelse(K == 3, 0, sum(zeta[cbind(sol_idx[1:(K-3)], sol_idx[2:(K-2)], sol_idx[4:K])]))
    }

    resK <- list(idx = sol_idx,
                 val = t[sol_idx],
                 loss_est = loss_est)
    if (only_maxK) return(resK) else res[[K]] <- resK
  }

  return(res)
}
