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
library(geodist) ## used in the earth quake example
library(ggplot2) ## ditto

## -----------------------------------------------------------------------------
California_centre <- data.frame(longitude = -119.449444, latitude = 37.166111)
distances <- geodist(globalearthquake[, c("longitude", "latitude")],
  California_centre,
  measure = "geodesic"
)
California_earth_quakes <- globalearthquake[distances < 2e6, ] ## distances are in meters
California_weeks <- rep(0, max(globalearthquake$nbweeks))
California_weeks[California_earth_quakes$nbweeks] <- 1
California_weeks_earth_quakes_model <- tune_vlmc(California_weeks,
  initial = "truncated",
  save = "all"
)
model <- as_vlmc(California_weeks_earth_quakes_model)
draw(model, prob = FALSE)

## -----------------------------------------------------------------------------
loglikelihood(model, initial = "truncated")

## -----------------------------------------------------------------------------
loglikelihood(model, initial = "specific")

## -----------------------------------------------------------------------------
loglikelihood(model, initial = "extended")

## ----fig.height=4-------------------------------------------------------------
CA_models <- c(
  list(California_weeks_earth_quakes_model$saved_models$initial),
  California_weeks_earth_quakes_model$saved_models$all
)
CA_extended <- data.frame(
  cutoff = sapply(CA_models, \(x) x$cutoff),
  loglikelihood = sapply(CA_models, loglikelihood,
    initial = "extended"
  )
)
ggplot(CA_extended, aes(cutoff, loglikelihood)) +
  geom_line() +
  xlab("Cut off (native scale)") +
  ylab("Log likelihood") +
  ggtitle("Extended log likelihood")

## ----fig.height=4-------------------------------------------------------------
CA_specific <- data.frame(
  cutoff = sapply(CA_models, \(x) x$cutoff),
  loglikelihood = sapply(CA_models, loglikelihood,
    initial = "specific"
  )
)
ggplot(CA_specific, aes(cutoff, loglikelihood)) +
  geom_line() +
  xlab("Cut off (native scale)") +
  ylab("Log likelihood") +
  ggtitle("Specific log likelihood")

## ----fig.height=4-------------------------------------------------------------
CW_combined <- rbind(
  CA_extended[c("cutoff", "loglikelihood")],
  CA_specific[c("cutoff", "loglikelihood")]
)
CW_combined[["Likelihood function"]] <- rep(c("extended", "specific"), times = rep(nrow(California_weeks_earth_quakes_model$results), 2))
ggplot(
  CW_combined,
  aes(cutoff, loglikelihood, color = `Likelihood function`)
) +
  geom_line() +
  xlab("Cut off (native scale)") +
  ylab("Log likelihood") +
  ggtitle("Log likelihood")

## -----------------------------------------------------------------------------
CA_model_extented <- tune_vlmc(California_weeks, initial = "extended")
model_extended <- as_vlmc(CA_model_extented)
draw(model_extended, prob = FALSE)

## -----------------------------------------------------------------------------
sun_activity <- as.factor(ifelse(sunspot.year >= median(sunspot.year), "high", "low"))
sun_model_tune_truncated <- tune_vlmc(sun_activity, initial = "truncated")
draw(as_vlmc(sun_model_tune_truncated))

## -----------------------------------------------------------------------------
sun_model_tune_extended <- tune_vlmc(sun_activity, initial = "extended")
draw(as_vlmc(sun_model_tune_extended))

## -----------------------------------------------------------------------------
TM0 <- matrix(c(0.7, 0.3, 0.4, 0.6),
  ncol = 2,
  byrow = TRUE
)
TM1 <- matrix(c(0.4, 0.6, 0.8, 0.2),
  ncol = 2,
  byrow = TRUE
)
init <- c(0, 1)
dts <- c(init, rep(NA, 500))
set.seed(0)
for (i in 3:length(dts)) {
  if (dts[i - 1] == 0) {
    probs <- TM0[dts[i - 2] + 1, ]
  } else {
    probs <- TM1[dts[i - 2] + 1, ]
  }
  dts[i] <- sample(0:1, 1, prob = probs)
}

## -----------------------------------------------------------------------------
MC_extended <- tune_vlmc(dts, initial = "extended", save = "all")
draw(as_vlmc(MC_extended))

## -----------------------------------------------------------------------------
MC_truncated <- tune_vlmc(dts, initial = "truncated", save = "all")
draw(as_vlmc(MC_truncated))

## -----------------------------------------------------------------------------
MC_specific <- tune_vlmc(dts, initial = "specific", save = "all")
draw(as_vlmc(MC_specific))

