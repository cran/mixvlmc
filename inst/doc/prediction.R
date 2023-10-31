## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  cache.path = ".prediction_cache/"
)

## ----setup--------------------------------------------------------------------
library(mixvlmc)
library(ggplot2)

## -----------------------------------------------------------------------------
sun_activity <- as.factor(ifelse(sunspot.year >= median(sunspot.year), "high", "low"))

## -----------------------------------------------------------------------------
sun_model_tune <- tune_vlmc(sun_activity)
sun_model_tune

## -----------------------------------------------------------------------------
sun_model <- as_vlmc(sun_model_tune)
sun_model_predictions <- predict(sun_model, sun_activity)

## -----------------------------------------------------------------------------
table(sun_model_predictions[-length(sun_model_predictions)], sun_activity)

## -----------------------------------------------------------------------------
sun_model_tune_aic <- tune_vlmc(sun_activity, criterion = "AIC")
sun_model_tune_aic

## -----------------------------------------------------------------------------
table(
  predict(as_vlmc(sun_model_tune_aic), sun_activity, final_pred = FALSE),
  sun_activity
)

## -----------------------------------------------------------------------------
first_half <- 1:(length(sun_activity) %/% 2)
sun_model_tune_aic_half <- tune_vlmc(sun_activity[first_half], criterion = "AIC")
sun_model <- as_vlmc(sun_model_tune_aic_half)
table(
  predict(sun_model, sun_activity[-first_half], final_pred = FALSE),
  sun_activity[-first_half]
)

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
CAC_covariates <- as.data.frame(EuStockMarkets)[c("DAX", "SMI", "FTSE")][-1, ]
CAC_covlmc <- tune_covlmc(CAC_dts, CAC_covariates, criterion = "AIC")
CAC_comodel <- as_covlmc(CAC_covlmc)

## -----------------------------------------------------------------------------
CAC_pred <- predict(CAC_comodel, CAC_dts, CAC_covariates, final_pred = FALSE)

## -----------------------------------------------------------------------------
table(CAC_pred, CAC_dts)

## ----cache=TRUE---------------------------------------------------------------
CAC_probs <- predict(CAC_comodel, CAC_dts, CAC_covariates, final_pred = FALSE, type = "probs")
CAC_probs[1:10, ]

## -----------------------------------------------------------------------------
entropies <- data.frame(entropy = apply(CAC_probs, 1, \(x) -sum(x * log(x))))
ggplot(entropies, aes(x = entropy)) +
  geom_density() +
  geom_rug(alpha = 0.1) +
  geom_vline(xintercept = -log(1 / 3), col = 2)

## -----------------------------------------------------------------------------
sun_probs <- predict(as_vlmc(sun_model_tune_aic), sun_activity,
  final_pred = FALSE,
  type = "probs"
)
sun_entropies <- data.frame(entropy = apply(
  sun_probs, 1,
  \(x) -sum(x * log(x), na.rm = TRUE)
))
ggplot(sun_entropies, aes(x = entropy)) +
  geom_histogram(bins = 50) +
  geom_rug(alpha = 0.5) +
  geom_vline(xintercept = -log(1 / 2), col = 2)

## -----------------------------------------------------------------------------
sun_metrics <- metrics(as_vlmc(sun_model_tune_aic))
sun_metrics

## ----fig.width=6, fig.height=6------------------------------------------------
plot(sun_metrics$roc)

## -----------------------------------------------------------------------------
CAC_metrics <- metrics(CAC_comodel)
CAC_metrics

