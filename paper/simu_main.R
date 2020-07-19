library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel()

library(survSched)
source('paper/utils.R')


## Simulation parameters ----------------------------------------
set.seed(222)
n <- 400
maxK <- 20
n_test <- 5000
n_simu <- 200
settings <- 1:3



## Run simulation and evaluate results ----------------------------------------

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
saveRDS(res_all, "paper/rds/simu_1to3_40020200.rds")



## Plots ----------------------------------------
# load stored results
res_all <- readRDS("paper/rds/simu_1to3_40020200.rds")

# Plot schedules (Figure 2)
Ks <- c(6, 10, 14)
sched_tb <- res_all[[1]]$sched_tib %>%
  filter(K %in% Ks) %>%
  mutate(rownum = row_number()) %>%
  gather(round, val, -rownum, -K) %>%
  group_by(rownum) %>%
  summarise(cen = median(val),
            y = sapply(cen, .hazfn),
            K = head(K, 1),
            hi = quantile(val, .975),
            lo = quantile(val, .025),
            hiq = quantile(val, .75),
            loq = quantile(val, .25))
haz_tb <- tibble(x = c(0, 2, 5, 13, 16),
                 y = sapply(x, .hazfn))
ggplot(sched_tb) +
  geom_path(aes(y, x), data = haz_tb) +
  geom_crossbar(aes(y, cen, ymin = loq, ymax = hiq), width = .02,
                color = 'orange') +
  geom_errorbar(aes(y, ymin = lo, ymax = loq), width = .02,
                color = 'orange') +
  geom_errorbar(aes(y, ymin = hiq, ymax = hi), width = .02,
                color = 'orange') +
  theme_bw() +
  labs(x = 'Hazard of T', y = 'Time (year)') +
  scale_y_continuous(breaks = seq(0, 16, by = 4)) +
  facet_wrap(~K, ncol = 3,
             labeller = purrr::partial(label_both, sep = ' = ')) +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(hjust = 0, size = 10)) +
  coord_flip()

ggsave("paper/fig/simu-sched-CI.pdf", width = 8, height = 2.6)


# plot mean, with percentile as confidence band (Figure 3)
tidy_eval_res <- do.call(bind_rows, lapply(res_all, function(r) r$eval_tib)) %>%
  select(detect_rate, lag_mean, loss_mean, K:setting)
res_current <- filter(tidy_eval_res, K == 0) %>%
  gather(key, val, detect_rate:loss_mean) %>%
  group_by(key, setting) %>%
  summarise(val_mean = mean(val))

tidy_eval_res %>%
  filter(K >= 4, K <= 20) %>%
  gather(key, val, detect_rate, lag_mean, loss_mean) %>%
  mutate(key = factor(key, levels = c("loss_mean", "detect_rate", "lag_mean"))) %>%
  group_by(K, method, setting, key) %>%
  summarise(mid = mean(val),
            hi = quantile(val, .975),
            lo = quantile(val, .025)) %>%
  ggplot(aes(K, mid, ymin = lo, ymax = hi)) +
  geom_hline(aes(yintercept = val_mean), data = res_current,
             color = "red", alpha = .8) +
  geom_hline(aes(yintercept = val), alpha = 0,
             data = data.frame(key = rep(c("loss_mean", "detect_rate", "lag_mean"), c(1, 2, 1)),
                               val = c(0, 0.6, 1, 0),
                               setting = 1)) +
  geom_ribbon(aes(fill = method),
              alpha = .22) +
  geom_point(aes(color = method),
             size = .6) +
  geom_line(aes(color = method, linetype = method),
            alpha = .8) +
  scale_color_manual(values=c("unif95"="#D55E00", "unif90"="#E69F00", "DP"="#0088CC")) +
  scale_fill_manual(values=c("unif95"="#D55E00", "unif90"="#E69F00", "DP"="#0088CC")) +
  facet_grid(key ~ setting,
             scales = "free_y", switch = "y",
             labeller = as_labeller(
               c(loss_mean= "Expected Loss",
                 detect_rate = "Detection Rate",
                 lag_mean = "Mean Lag Time",
                 `1` = "Setting 1", `2` = "Setting 2", `3` = "Setting 3")
             )) +
  ylab(NULL) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_x_continuous(breaks = seq(4, 20, by = 4))

ggsave('paper/fig/simu-eval.pdf', width = 8, height = 8)
