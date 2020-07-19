library(tidyverse)
source("utils.R")

# Evaluation on quantiles -------------------------------------------------

res_all <- read_rds("rds/simu_1to3_40020200.rds")

tidy_eval_res <- do.call(bind_rows, lapply(res_all, function(r) r$eval_tib)) %>%
  select(lag_q1, loss_q1, lag_q3, loss_q3, K:setting) %>%
  gather(key, val, lag_q1, loss_q1, lag_q3, loss_q3) %>%
  mutate(type = ifelse(startsWith(key, "lag"), "lag", "loss"),
         quantity = ifelse(endsWith(key, "q1"), "1st quantile", "3rd quantile"))

res_current <- filter(tidy_eval_res, K == 0) %>%
  group_by(key, setting, type, quantity) %>%
  summarise(val_mean = mean(val))

tidy_eval_res %>%
  filter(K >= 4, K <= 20) %>%
  group_by(K, method, setting, key, type, quantity) %>%
  summarise(mid = mean(val),
            hi = quantile(val, .975),
            lo = quantile(val, .025)) %>%
  ggplot(aes(K, mid, ymin = lo, ymax = hi)) +
  # value achieved current guideline: red horizontal line
  geom_hline(aes(yintercept = val_mean), data = res_current,
             color = "red", alpha = .8) +
  # set ylim
  geom_hline(aes(yintercept = val), alpha = 0,
             data = data.frame(type = c("loss", "lag"),
                               val = c(0, 0.6, 1, 0),
                               setting = 1)) +
  # confidence band, dots and curves
  geom_ribbon(aes(fill = method, group = interaction(key, method)),
              alpha = .22) +
  geom_point(aes(color = method, shape = quantity),
             size = 1.2, alpha = .8) +
  geom_line(aes(color = method, linetype = method, group = interaction(key, method)),
            alpha = .8) +
  scale_color_manual(values=c("unif95"="#D55E00", "unif90"="#E69F00", "DP"="#0088CC")) +
  scale_fill_manual(values=c("unif95"="#D55E00", "unif90"="#E69F00", "DP"="#0088CC")) +
  # other plotting configurations
  facet_grid(type ~ setting,
             scales = "free_y", switch = "y",
             labeller = as_labeller(
               c(loss= "Quantiles of loss (years)",
                 lag = "Quantiles of lag time (years)",
                 `1` = "Setting 1", `2` = "Setting 2", `3` = "Setting 3")
             )) +
  xlab("Number of visits") +
  ylab(NULL) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_x_continuous(breaks = seq(4, 20, by = 4))

ggsave("fig/simu-supp-eval-quantile.pdf", width = 8, height = 6)





# Disturbance -------------------------------------------------------------

Ks <- c(6, 10, 14)

res_delay <- readRDS("rds/simu_supp_delay.rds")
res_undetect <- readRDS("rds/simu_supp_undetect.rds")
res_adjustu <- readRDS("rds/simu_supp_adjustu.rds")

labx <- function(s) {
  switch(s,
         delay = "Delay",
         undetect = "Undetected proportion",
         adjustu = "Risk reduction r")
}

# top panel: distance between schedules
for (s in c("delay", "undetect", "adjustu")) {
  res_dist <- gather(get(paste0("res_", s)),
                     prop, val, -`0`, -round, -K) %>%
    group_by(round, K, prop) %>%
    summarise(dist = lqdist_sched(`0`, val, M = 16, q = 1))
  p <- ggplot(res_dist) +
    geom_boxplot(aes(prop, dist)) +
    facet_wrap( ~ K, labeller = label_both) +
    theme_bw() +
    {
      if (s == "adjustu")
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    } +
    labs(x = labx(s),
         y = 'Wasserstein distance')
  ggsave(paste0("fig/simu-supp-dist-", s, ".pdf"),
         plot = p,
         width = 8, height = 2.8)
}


# bottom panel: schedule of one round
for (s in c("delay", "undetect", "adjustu")) {
  p <- get(paste0("res_", s)) %>%
    filter(round == 100) %>%
    mutate(No = row_number()) %>%
    gather(prop, val, -No, -K, -round) %>%
    mutate(prop = as.numeric(prop)) %>%
    ggplot(aes(prop, val, group = No)) +
    geom_point(size = .7) +
    geom_line(alpha = .25) +
    theme_bw() +
    facet_wrap(~K, labeller = label_both) +
    ylim(0, 16) +
    labs(y = "Schedule time (year)",
         x = labx(s))
  ggsave(paste0("fig/simu-supp-sched-", s, ".pdf"),
         plot = p,
         width = 8, height = 2.8)
}
