library(TPmsm)
library(dplyr)

# Generate training set as a survTP
make_data <- function(n, setting = 1,
                      censored = TRUE, cen_setting = 1,
                      returnTU = FALSE, ...) {
  TU <- .generate_TU(n, setting)
  if(returnTU) return(TU)
  if (censored) C <- .generate_C(n, cen_setting) else C <- Inf
  .TU_to_train(TU, C, ...)
}

.hazfn <- function(t, x, betas) {
  h0 <- if (t < 2) 2 * t else
    if (t < 5) 6 - t else
      if (t < 13) 1.5 - .1 * t else .2
  .15 * h0
}

.generate_TU <- function(n, setting = 1) {
  if (setting == 0) {
    t <- sample(c(2, 3, 4, 6, 8, 10), n, replace = TRUE)
    u1 <- rweibull(n, 2, 15)
    u <- pmin(u1, t + rexp(n, .1))
  }
  if (setting == 1) {
    t <- simsurv::simsurv(hazard = .hazfn, x = data.frame(rep(NA, n)))$eventtime
    u1 <- rweibull(n, 2, 15)
    u <- pmin(u1, t + rexp(n, .1))
  }
  if (setting == 2) {
    t <- simsurv::simsurv(hazard = .hazfn, x = data.frame(rep(NA, n)))$eventtime
    u1 <- rweibull(n, 2, 15)
    u <- pmin(u1, t + rexp(n, .2))
  }
  if (setting == 3) {
    t <- simsurv::simsurv(hazard = .hazfn, x = data.frame(rep(NA, n)))$eventtime
    u <- rep(NA, n)
    i <- 1
    while(i <= n) {
      u[i] <- rweibull(1, 2, 8)
      if (is.infinite(t[i]) ||
          (u[i] < t[i] && runif(1) < .2) ||
          (u[i] > t[i] && runif(1) < .8)) i <- i + 1
    }
  }
  t[t > 16] <- Inf

  cbind(t = t, u = u)
}

.generate_C <- function(n, setting = 1) {
  if (setting == 1) return(runif(n, 0, 16))
  if (setting == 2) return(16 - runif(n, 0, 4)^2)
  if (setting == 3) truncdist::rtrunc(n, 'weibull', b = 16, shape = 1.7, scale = 7)
}

# Convert (T, U, C) to survTP, optionally with disturbance
.TU_to_train <- function(TU, C,
                         delay_detect = 0, prop_undetect = 0) {
  # delay detection by delay_detect * (U - T)
  if (delay_detect > 0) {
    TU[,1] <- (1 - delay_detect) * TU[,1] + delay_detect * TU[,2]
  }
  # make undetected if T + pesudo_delay > U
  if (prop_undetect > 0) {
    rate0 <- 1 / mean((TU[,2] - TU[,1])[TU[,2] > TU[,1]])
    pesudo_delay <- rexp(nrow(TU), rate0 * (1 - prop_undetect) / prop_undetect)
    make_undetect <- (TU[,1] + pesudo_delay) > TU[,2]
    TU[,1] <- ifelse(make_undetect, Inf, TU[,1])
  }
  # convert to survTP format
  minTU <- pmin(TU[,1], TU[,2])
  event1 <- minTU < C
  time1 <- ifelse(event1, minTU, C)
  event <- TU[,2] < C
  Stime <- ifelse(event, TU[,2], C)

  TPmsm::survTP(time1, event1, Stime, event)
}

# Evaluate loss on a uncensored test set (containing only T and U)
eval_on_test <- function(sched_t, TU){
  t <- TU[,1]; u <- TU[,2]
  recur <- t < u
  a <- rep(NA, length(t))
  for (i in (1:length(t))[recur]) a[i] <- sched_t[sched_t >= t[i]][1]
  detect <- recur & !is.na(a) & (a <= u)
  loss <- ifelse(recur,
                 ifelse(detect, a - t, u - t),
                 0)

  s <- c("min", "q1", "median", "mean", "q3", "max")
  lag_summary <- setNames(summary((a-t)[detect]), paste0("lag_", s))
  loss_summary <- setNames(summary(loss), paste0("loss_", s))

  c(recur_rate = mean(recur),
    detect_rate = sum(detect) / sum(recur),
    lag_summary,
    loss_summary)
}

# Calculate distance between two schedules
lqdist_sched <- function(sched_val1, sched_val2, M, q = 1) {
  sched_val1 <- drop(sched_val1)
  sched_val2 <- drop(sched_val2)
  K1 <- length(sched_val1)
  K2 <- length(sched_val2)
  tibble(sched_val = c(sched_val1, sched_val2),
         jump = c(rep(1 / K1, K1), rep(-1 / K2, K2))) %>%
    arrange(sched_val) %>%
    mutate(diff_ecdf = abs(cumsum(jump)),
           len_intvl = lead(sched_val) - sched_val) %>%
    summarise(dist = (sum(diff_ecdf ^ q * len_intvl, na.rm = TRUE) / M) ^
                (1 / q)) %>%
    unlist() %>%
    setNames(NULL)
}

# Impute U
adjust_U <- function(survdat, haz_drop) {
  Stime_adj <- with(survdat[[1]],
                    ifelse(event1 == 1 & time1 < Stime,
                           time1 + (1 - haz_drop) * (Stime - time1),
                           Stime))
  with(survdat[[1]],
       survTP(time1, event1, Stime_adj, event))
}

# Gernerate T and U given transition probability
make_tu_with_tp <- function(n, tp) {
  t(sapply(1:n, function(x) .make_tu_with_tp(tp)))
}

.make_tu_with_tp <- function(tp) {
  min_tu <- .invF(tp$p00_0, tp$time)
  if (is.na(min_tu$idx)) return(c(Inf, Inf))

  p0 <- ifelse(min_tu$idx > 1,
               tp$p00_0[min_tu$idx] / tp$p00_0[min_tu$idx - 1],
               tp$p00_0[1])
  p1 <- ifelse(min_tu$idx > 1,
               tp$p01[min_tu$idx - 1, min_tu$idx],
               tp$p01_0[1])
  if (p0 + p1 > 1) {
    warning("Probability > 1!")
    return(c(0, 0))
  }
  if (runif(1) > p1 / (1 - p0)) {
    return(c(Inf, min_tu$val))
  }

  t <- min_tu$val
  u <- .invF(tp$p11[min_tu$idx, ], tp$time)$val
  c(t, u)
}

.invF <- function(p, t) {
  r <- runif(1)
  idx <- which(p <= r)[1]
  val <- ifelse(is.na(idx), Inf, t[idx])
  list(idx = idx, val = val)
}



