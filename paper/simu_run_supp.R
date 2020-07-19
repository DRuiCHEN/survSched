library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel()

library(survSched)
source('utils.R')



# Disturbance: run experiments --------------------------------------------
# Simulation parameters
set.seed(222)
n <- 400
n_simu <- 100

setting <- 1
Ks <- c(6, 10, 14)

# Delayed detection (about 150 s)
system.time({
  prop <- seq(0, .3, by = .05)
  res_delay <- foreach(i = 1:n_simu, .combine = 'bind_rows') %dopar% {
    TU <- .generate_TU(n)
    C <- .generate_C(n)
    foreach(p = prop, .combine = 'cbind') %do% {
      dpres <- .TU_to_train(TU, C, delay_detect = p) %>%
        compute_tp(end_time = 16) %>%
        dp_proc(max(Ks))
      foreach(k = Ks, .combine = 'c') %do% find_schedule(dpres, k)$val
    } %>%
      as_tibble() %>%
      setNames(prop) %>%
      mutate(K = rep(Ks, Ks),
             round = rep(i, length(K)))
  }
})
saveRDS(res_delay, "rds/simu_supp_delay.rds")


# Missed detection about 150 s)
system.time({
  prop <- seq(0, .3, by = .05)
  res_undetect <- foreach(i = 1:n_simu, .combine = 'bind_rows') %dopar% {
    TU <- .generate_TU(n)
    C <- .generate_C(n)
    foreach(p = prop, .combine = 'cbind') %do% {
      dpres <- .TU_to_train(TU, C, prop_undetect = p) %>%
        compute_tp(end_time = 16) %>%
        dp_proc(max(Ks))
      foreach(k = Ks, .combine = 'c') %do% find_schedule(dpres, k)$val
    } %>%
      as_tibble() %>%
      setNames(prop) %>%
      mutate(K = rep(Ks, Ks),
             round = rep(i, length(K)))
  }
})
saveRDS(res_undetect, "rds/simu_supp_undetect.rds")


# Imputing U (about 180 s)
system.time({
  haz_drop <- seq(0, .6, by = .1)
  res_adjustu <- foreach(i = 1:n_simu, .combine = 'bind_rows') %dopar% {
    dtrain <- make_data(n)
    foreach(p = haz_drop, .combine = 'cbind') %do% {
      dpres <- adjust_U(dtrain, p) %>%
        compute_tp(end_time = 16) %>%
        dp_proc(max(Ks))
      foreach(k = Ks, .combine = 'c') %do% find_schedule(dpres, k)$val
    } %>%
      as_tibble() %>%
      setNames(haz_drop) %>%
      mutate(K = rep(Ks, Ks),
             round = rep(i, length(K)))
  }
})
saveRDS(res_adjustu, "rds/simu_supp_adjustu.rds")
