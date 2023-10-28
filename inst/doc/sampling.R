## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 7,
  fig.align = "center",
  cache.path = ".sampling_cache/"
)

## ----setup, message=FALSE-----------------------------------------------------
library(mixvlmc)
library(ggplot2)
library(data.table) ## for frollapply

## -----------------------------------------------------------------------------
dts <- sample(0:1, 500, replace = TRUE)

## -----------------------------------------------------------------------------
best_vlmc <- tune_vlmc(dts)
draw(as_vlmc(best_vlmc))

## -----------------------------------------------------------------------------
dts_sample <- simulate(as_vlmc(best_vlmc), 600)[-(1:100)]

## -----------------------------------------------------------------------------
dts_sample_2 <- simulate(as_vlmc(best_vlmc), 10, init = c(0L, 0L))

## -----------------------------------------------------------------------------
dts_sample_2

## -----------------------------------------------------------------------------
dts_sample_3 <- simulate(as_vlmc(best_vlmc), 10, burnin = 20)

## -----------------------------------------------------------------------------
CAC_raw <- as.data.frame(EuStockMarkets)$CAC

## -----------------------------------------------------------------------------
CAC_rel_evol <- diff(CAC_raw) / CAC_raw[-length(CAC_raw)]
CAC_dts <- factor(
  ifelse(CAC_rel_evol >= 0.005, "Up",
    ifelse(CAC_rel_evol <= -0.005, "Down", "Stay")
  ),
  levels = c("Down", "Stay", "Up")
)

## ----cache=TRUE---------------------------------------------------------------
CAC_vlmc <- tune_vlmc(CAC_dts, criterion = "AIC")
CAC_model <- as_vlmc(CAC_vlmc)

## ----fig.height=4-------------------------------------------------------------
CAC_rle <- rle(as.integer(CAC_dts))
CAC_rle_df <- data.frame(value = CAC_rle$values, length = CAC_rle$lengths)
ggplot(CAC_rle_df, aes(x = length)) +
  geom_bar() +
  labs(
    title = "Distribution of the lengths of constant subsequences",
    x = "Length",
    y = "Count"
  )

## -----------------------------------------------------------------------------
long_sequence <- function(dts) {
  dts_int <- as.integer(dts)
  ## RLE cannot be used directly as we need to account for overlapping
  ## patterns
  dts_freq <- frollapply(dts_int, 10, \(x) max(rle(x)$length) >= 5)
  mean(dts_freq, na.rm = TRUE)
}

## ----cache=TRUE---------------------------------------------------------------
bootstrap_samples <- vector(100, mode = "list")
burn_in_time <- 50 * depth(CAC_model)
for (b in seq_along(bootstrap_samples)) {
  bootstrap_samples[[b]] <- simulate(CAC_model,
    length(CAC_dts),
    burin = burn_in_time,
    seed = b
  )
}

## ----cache=TRUE---------------------------------------------------------------
bootstrap_ls <- sapply(bootstrap_samples, long_sequence)

## ----fig.height=4-------------------------------------------------------------
ggplot(mapping = aes(x = bootstrap_ls)) +
  geom_density() +
  geom_rug() +
  labs(
    title = "Bootstrap distribution of the probability of long constant subsequences",
    x = "Probability"
  ) +
  geom_vline(xintercept = long_sequence(CAC_dts), color = "red")

## ----cache=TRUE---------------------------------------------------------------
CAC_covariates <- as.data.frame(EuStockMarkets)[c("DAX", "SMI", "FTSE")][-1, ]
CAC_covlmc <- tune_covlmc(CAC_dts, CAC_covariates, criterion = "AIC")
CAC_comodel <- as_covlmc(CAC_covlmc)

## -----------------------------------------------------------------------------
draw(CAC_comodel, model = "full", with_state = TRUE)

## ----cache=TRUE---------------------------------------------------------------
co_sample <- simulate(CAC_comodel,
  length(CAC_dts),
  seed = 0,
  init = CAC_dts[1:depth(CAC_comodel)],
  covariate = CAC_covariates
)

## ----fig.height=5-------------------------------------------------------------
pc_week <- powerconsumption[powerconsumption$week == 12, ]
elec <- pc_week$active_power
ggplot(pc_week, aes(x = date_time, y = active_power)) +
  geom_line() +
  xlab("Date") +
  ylab("Activer power (kW)")

## -----------------------------------------------------------------------------
elec_dts <- cut(elec, breaks = c(0, 1.2, 8), labels = c("background", "active"))
elec_cov <- data.frame(day = (pc_week$hour >= 7 & pc_week$hour <= 17))

## -----------------------------------------------------------------------------
elec_covlmc_tune <- tune_covlmc(elec_dts, elec_cov, criterion = "AIC")
best_elec_covlmc <- as_covlmc(elec_covlmc_tune)
draw(best_elec_covlmc, model = "full", time_sep = " | ", p_value = FALSE)

## -----------------------------------------------------------------------------
day_only <- simulate(best_elec_covlmc,
  length(elec_dts),
  seed = 0,
  covariate = data.frame(day = rep(TRUE, length(elec_dts)))
)

## -----------------------------------------------------------------------------
night_only <- simulate(best_elec_covlmc,
  length(elec_dts),
  seed = 0,
  covariate = data.frame(day = rep(FALSE, length(elec_dts)))
)

