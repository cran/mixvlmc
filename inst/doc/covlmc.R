## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 7,
  fig.align = "center"
)

## ----setup--------------------------------------------------------------------
library(mixvlmc)
library(ggplot2)

## ----active_power_week_15_16, fig.height=5------------------------------------
pc_week_15_16 <- powerconsumption[powerconsumption$week %in% c(15, 16), ]
elec <- pc_week_15_16$active_power
ggplot(pc_week_15_16, aes(x = date_time, y = active_power)) +
  geom_line() +
  xlab("Date") +
  ylab("Activer power (kW)")

## -----------------------------------------------------------------------------
elec_dts <- cut(elec, breaks = c(0, 0.4, 2, 8), labels = c("low", "typical", "high"))

## -----------------------------------------------------------------------------
elec_vlmc_tune <- tune_vlmc(elec_dts)
best_elec_vlmc <- as_vlmc(elec_vlmc_tune)
draw(best_elec_vlmc)

## -----------------------------------------------------------------------------
set.seed(0)
nb_obs <- 200
covariates <- data.frame(y = runif(nb_obs))
x <- 0
for (k in 2:nb_obs) {
  ## we induce a simple dependency to the covariate
  ## and an order 1 memory
  if (covariates$y[k - 1] < 0.5) {
    if (x[k - 1] == 0) {
      x[k] <- sample(c(0, 1), 1, prob = c(0.7, 0.3))
    } else {
      x[k] <- sample(c(0, 1), 1, prob = c(0.3, 0.7))
    }
  } else {
    if (x[k - 1] == 0) {
      x[k] <- sample(c(0, 1), 1, prob = c(0.1, 0.9))
    } else {
      x[k] <- sample(c(0, 1), 1, prob = c(0.5, 0.5))
    }
  }
}

## -----------------------------------------------------------------------------
model <- covlmc(x, covariates)
model

## -----------------------------------------------------------------------------
draw(model, model = "full", p_value = FALSE)

## -----------------------------------------------------------------------------
model_tune <- tune_covlmc(x, covariates)
model_tune

## -----------------------------------------------------------------------------
draw(as_covlmc(model_tune), model = "full", p_value = FALSE)

## -----------------------------------------------------------------------------
model_tune_2 <- tune_covlmc(x, covariates, min_size = 2)
model_tune_2

## -----------------------------------------------------------------------------
model_tune_10 <- tune_covlmc(x, covariates, min_size = 10)
model_tune_10

## -----------------------------------------------------------------------------
elec_cov <- data.frame(day = (pc_week_15_16$hour >= 7 & pc_week_15_16$hour <= 18))

## -----------------------------------------------------------------------------
elec_tune <- tune_covlmc(elec_dts, elec_cov)
elec_tune

## ----fig.height=4-------------------------------------------------------------
ggplot(elec_tune$results, aes(x = alpha, y = BIC)) +
  geom_line() +
  geom_point()

## -----------------------------------------------------------------------------
print(autoplot(elec_tune))

## -----------------------------------------------------------------------------
draw(as_covlmc(elec_tune), model = "full", p_value = FALSE, with_state = TRUE)

## -----------------------------------------------------------------------------
elec_tune_3 <- tune_covlmc(elec_dts, elec_cov, min_size = 3)
elec_tune_3

## -----------------------------------------------------------------------------
elec_tune_10 <- tune_covlmc(elec_dts, elec_cov, min_size = 10)
elec_tune_10

## -----------------------------------------------------------------------------
elec_model <- as_covlmc(elec_tune)
states(elec_model)
depth(elec_model)
covariate_depth(elec_model)
context_number(elec_model)

## -----------------------------------------------------------------------------
logLik(elec_model)
AIC(elec_model)
BIC(elec_model)

## -----------------------------------------------------------------------------
draw(elec_model)

## -----------------------------------------------------------------------------
contexts(elec_model, model = "coef")

## -----------------------------------------------------------------------------
contexts(elec_model, model = "full")

## -----------------------------------------------------------------------------
ctxs <- contexts(elec_model)
ctxs

## -----------------------------------------------------------------------------
model(ctxs[[2]])
model(ctxs[[3]], type = "full")

