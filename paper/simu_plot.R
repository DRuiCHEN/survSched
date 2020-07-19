library(tidyverse)
source("utils.R")


# load stored results
res_all <- readRDS("rds/simu_1to3_40020200.rds")

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

ggsave("fig/simu-sched-CI.pdf", width = 8, height = 2.6)


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

ggsave('fig/simu-eval.pdf', width = 8, height = 8)
