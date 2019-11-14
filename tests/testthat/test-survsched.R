library(magrittr)

test_that("Optimal schedule found", {
  n <- 300
  K <- 4
  t <- rweibull(n, 5, 5 / gamma(1 + 1.2))
  u1 <- rweibull(n, 5, 8 / gamma(1 + 1.2))
  u <- ifelse(u1 <= t, u1, t + rexp(n) * 4)
  cen <- runif(n, 0, 16)
  minTU <- pmin(t, u)
  event1 <- minTU < cen
  time1 <- pmin(minTU, cen)
  event <- u < cen
  Stime <- pmin(u, cen)
  survdat <- TPmsm::survTP(time1, event1, Stime, event)

  tp <- compute_tp(survdat, t_gap = .8)
  sched_dp <- dp_proc(tp, K) %>% find_schedule()

  sched_all <- combn(length(tp$time), K)
  loss_all <- apply(sched_all, 2,
                    function(idx) eval_loss(idx, tp))
  sched_opt <- sched_all[, which.min(loss_all)]

  expect_equal(sched_opt, sched_dp$idx)

})
