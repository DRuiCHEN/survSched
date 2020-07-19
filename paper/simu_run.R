library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel()

library(survSched)
source('utils.R')


# Simulation parameters ---------------------------------------------------
set.seed(222)
n <- 400
maxK <- 20
n_test <- 5000
n_simu <- 200
settings <- 1:3


# Run simulation and evaluate results -------------------------------------

# single run
do_run <- function(n, n_test, setting, maxK, sched_current){
  # generate data training and test datasets
  dtrain <- make_data(n, setting)
  dtest <- make_data(n_test, setting, returnTU = TRUE)
  # DP
  dp_res <- dp_proc(compute_tp(dtrain, end_time = 16), maxK = maxK)
  # to record DP schedules as a long vector
  sched_vec <- numeric(0)

  # evaluate results for three methods with k varying from 1 to maxK
  eval_res <- foreach(method = c('unif90', 'unif95', 'DP'), .combine = 'bind_rows') %do% {
    tmp <- foreach(k = 1:maxK, .combine = 'rbind') %do% {
      sched <- switch(
        method,
        DP = find_schedule(dp_res, k)$val,
        unif90 = seq(0, quantile(dtrain[[1]]$Stime, .9), length.out=k+1)[-1],
        unif95 = seq(0, quantile(dtrain[[1]]$Stime, .95), length.out=k+1)[-1]
      )
      # record DP schedules
      if (method == 'DP') sched_vec <- c(sched_vec, sched)
      # evaluate schedules
      eval_on_test(sched, dtest)
    }
    tmp <- as_tibble(tmp) %>%
      mutate(K = 1:maxK, method = !!method)
    tmp
  }

  if (!missing(sched_current)) {
    eval_res <- do.call(function(...) add_row(eval_res, ...),
                        c(as.list(eval_on_test(sched_current, dtest)),
                          K = 0, method = "current"))
  }

  list(dp_sched_vec = sched_vec,
       eval_res = eval_res)
}


# takes about 700 seconds (in parallel with 4 cores)
system.time({
  res_all <- foreach(setting = settings) %do% {
    cat("Setting ", setting, "\n")

    res <- foreach(i = 1:n_simu) %dopar% {  # use dopar for parallel computing
      cat("\tRound ", i, "\n")
      do_run(n, n_test, setting, maxK,
             sched_current = c(seq(0, 5, by = .5)[-1], 6:15))
    }

    list(
      sched_tib = foreach(r = res, .combine = "cbind") %do% {r$dp_sched_vec} %>%
        as_tibble() %>%
        setNames(paste0('round', 1:n_simu)) %>%
        mutate(K = rep(1:maxK, 1:maxK)),
      eval_tib = foreach(r = res, .combine = 'bind_rows') %do% {r$eval_res} %>%
        mutate(setting = setting)
    )
  }
})

# store the result
saveRDS(res_all, "rds/simu_1to3_40020200.rds")



