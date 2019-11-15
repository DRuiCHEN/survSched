library(tidyverse)


## Main result, evaluation on quantiles ----------------------------------------
res_all <- read_rds("rds/simu_1to3_40020200.rds")

tidy_eval_res <- do.call(bind_rows, lapply(res_all, function(r) r$eval_tib)) %>%
  select(lag_q1, loss_q1, lag_q3, loss_q3, K:setting)

tidy_eval_res %>%
  filter(K >= 4, K <= 20) %>%
  gather(key, val, lag_q1, loss_q1, lag_q3, loss_q3) %>%
  mutate(type = ifelse(startsWith(key, "lag"), "lag", "loss"),
         quantity = ifelse(endsWith(key, "q1"), "1st quantile", "3rd quantile")) %>%
  group_by(K, method, setting, key, type, quantity) %>%
  summarise(mid = mean(val),
            hi = quantile(val, .975),
            lo = quantile(val, .025)) %>%
  ggplot(aes(K, mid, ymin = lo, ymax = hi)) +
  geom_hline(aes(yintercept = val), alpha = 0,
             data = data.frame(type = c("loss", "lag"),
                               val = c(0, 0.6, 1, 0),
                               setting = 1)) +
  geom_ribbon(aes(fill = method, group = interaction(key, method)),
              alpha = .22) +
  geom_point(aes(color = method, shape = quantity),
             size = 1.2, alpha = .8) +
  geom_line(aes(color = method, linetype = method, group = interaction(key, method)),
            alpha = .8) +
  scale_color_manual(values=c("unif95"="#D55E00", "unif90"="#E69F00", "DP"="#0088CC")) +
  scale_fill_manual(values=c("unif95"="#D55E00", "unif90"="#E69F00", "DP"="#0088CC")) +
  facet_grid(type ~ setting,
             scales = "free_y", switch = "y",
             labeller = as_labeller(
               c(loss= "Expected Loss",
                 lag = "Lag Time",
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

ggsave("fig/simu-eval-quantile.pdf", width = 8, height = 6)


## High order ----------------------------------------



## Disturbance ----------------------------------------

