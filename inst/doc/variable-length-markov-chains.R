## ---- include = FALSE---------------------------------------------------------
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
library(ggplot2) ## used in the earth quake example

## ---- echo=FALSE--------------------------------------------------------------
homc <- data.frame(m = 1:10)
homc$parameters <- 3 * (4^homc$m)
homc

## ----echo=FALSE---------------------------------------------------------------
bin_mark <- cbind(
  data.frame(Probablity = c(0.1, 0.2, 0.1, 0.3, 0.1, 0.4, 0.1, 0.3)),
  expand.grid("n-1" = 0:1, "n-2" = 0:1, "n-3" = 0:1)
)

## ---- echo=FALSE--------------------------------------------------------------
bin_mark

## -----------------------------------------------------------------------------
set.seed(0)
x <- sample(c(0L, 1L, 2L), 200, replace = TRUE)
model <- vlmc(x)
model

## -----------------------------------------------------------------------------
model_theo <- vlmc(x, cutoff = log(length(x)))
model_theo

## -----------------------------------------------------------------------------
model_large <- vlmc(x, cutoff = 0.5 * log(length(x)))
model_large
model_cutoff <- cutoff(model_large, mode = "native")
model_cutoff

## -----------------------------------------------------------------------------
model_medium <- prune(model_large, cutoff = model_cutoff[1])
model_medium

## -----------------------------------------------------------------------------
model_small <- prune(model_large, cutoff = model_cutoff[2])
model_small

## -----------------------------------------------------------------------------
model_tune <- tune_vlmc(x)
model_opt <- as_vlmc(model_tune)
model_opt

## -----------------------------------------------------------------------------
California_centre <- data.frame(longitude = -119.449444, latitude = 37.166111)
distances <- geodist(globalearthquake[, c("longitude", "latitude")],
  California_centre,
  measure = "geodesic"
)
California_earth_quakes <- globalearthquake[distances < 2e6, ] ## distances are in meters

## -----------------------------------------------------------------------------
California_weeks <- rep(0, max(globalearthquake$nbweeks))
California_weeks[California_earth_quakes$nbweeks] <- 1

## -----------------------------------------------------------------------------
California_weeks_earth_quakes_model <- tune_vlmc(California_weeks)

## -----------------------------------------------------------------------------
draw(as_vlmc(California_weeks_earth_quakes_model))

## -----------------------------------------------------------------------------
summary(California_weeks_earth_quakes_model)

## ----fig.height=4-------------------------------------------------------------
ggplot(California_weeks_earth_quakes_model$results, aes(x = alpha, y = BIC)) +
  geom_line() +
  geom_point()

## -----------------------------------------------------------------------------
states(model_large)
depth(model_large)
context_number(model_large)

## -----------------------------------------------------------------------------
logLik(model_large)
AIC(model_large)
BIC(model_large)

## -----------------------------------------------------------------------------
draw(model_large)

## -----------------------------------------------------------------------------
contexts(model_large, cutoff = "native")

## -----------------------------------------------------------------------------
contexts(model_large, cutoff = "quantile", reverse = FALSE, frequency = "detailed")

