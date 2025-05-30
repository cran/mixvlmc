glm_warning_ignore <- function(w) {
  to_ignore <- stringr::coll(
    c(
      gettext(c(
        "glm.fit: fitted probabilities numerically 0 or 1 occurred",
        "glm.fit: algorithm did not converge"
      ), domain = "R-stats"),
      ## keep the English text to circumvent inconsistency in locale setting
      "glm.fit: fitted probabilities numerically 0 or 1 occurred",
      "glm.fit: algorithm did not converge"
    )
  )
  for (msg in to_ignore) {
    if (stringr::str_detect(conditionMessage(w), msg)) {
      rlang::cnd_muffle(w)
    }
  }
}

vgam_warning_ignore <- function(w) {
  to_ignore <- list(
    ".* diagonal elements of the working weights variable 'wz' have been replaced by",
    stringr::coll("fitted values close to 0 or 1"),
    stringr::coll("fitted probabilities numerically 0 or 1 occurred"),
    stringr::coll("some quantities such as z, residuals, SEs may be inaccurate due to convergence at a half-step"),
    stringr::coll("iterations terminated because half-step sizes are very small"),
    stringr::coll("there are NAs here in slot linkinv"),
    stringr::coll("there are NAs in eta in slot inverse")
  )
  for (msg in to_ignore) {
    if (stringr::str_detect(conditionMessage(w), msg)) {
      rlang::cnd_muffle(w)
    }
  }
}

multinom_warning_ignore_generator <- function(target, target_dist) {
  missing <- levels(target)[target_dist == 0L]
  if (length(missing) > 0) {
    the_msg <- stringr::fixed(
      stringr::str_replace(
        ngettext(length(missing),
          "group %s is empty",
          "groups %s are empty",
          domain = "R-nnet"
        ), stringr::fixed("%s"),
        paste(sQuote(missing), collapse = " ")
      )
    )
    function(w) {
      if (stringr::str_detect(conditionMessage(w), the_msg)) {
        rlang::cnd_muffle(w)
      }
    }
  } else {
    function(w) { }
  }
}

fit_glm <- function(target, mm, nb_vals, control) {
  assertthat::assert_that(nrow(mm) > 0)
  engine <- options()[["mixvlmc.predictive"]]
  assertthat::assert_that(engine %in% c("glm", "multinom"))
  target_dist <- table(target)
  non_empty <- target_dist[target_dist > 0]
  if (length(non_empty) == 1) {
    ## degenerate case
    constant_model(target, mm, nb_vals, control$pseudo_obs)
  } else {
    if (engine == "glm") {
      if (nb_vals == 2) {
        if (ncol(mm) > 0) {
          mm$target <- target
          try_glm <- try(
            withCallingHandlers(
              warning = glm_warning_ignore,
              result <-
                stats::glm(target ~ .,
                  data = mm, family = stats::binomial(),
                  x = FALSE, y = FALSE,
                  model = FALSE, control = stats::glm.control(maxit = options()[["mixvlmc.maxit"]])
                )
            ),
            silent = TRUE
          )
          if (inherits(try_glm, "try-error")) {
            err_cond <- as.character(attr(try_glm, "condition"))
            if (stringr::str_detect(
              err_cond,
              stringr::coll(gettext("contrasts can be applied only to factors with 2 or more levels",
                domain = "R-stats"
              ))
            ) || stringr::str_detect(
              err_cond, "contrasts can be applied only to factors with 2 or more levels"
            )) {
              ## fake result, interpreted as a low rank result
              result <- structure(list(coefficients = c(NA), ll = NA, rank = 0, target = NA), class = "constant_model")
            } else {
              stop(attr(try_glm, "condition"))
            }
          }
        } else {
          result <-
            stats::glm(target ~ 1,
              family = stats::binomial(),
              x = FALSE, y = FALSE,
              model = FALSE
            )
        }
        if (inherits(result, "glm")) {
          if (!is_glm_low_rank(result)) {
            ## check only convergence for full rank models
            if (!result$converged) {
              warning("glm.fit did not converge")
            }
          } else {
            ## signal non convergence
            if (!result$converged) {
              message("glm.fit did not converge for a discarded low rank model")
            }
          }
        }
      } else {
        if (ncol(mm) > 0) {
          mm$target <- target
          try_vglm <- try(
            withCallingHandlers(
              warning = vgam_warning_ignore,
              result <- VGAM::vglm(target ~ .,
                data = mm, family = VGAM::multinomial(refLevel = 1),
                x.arg = FALSE, y.arg = FALSE, model = FALSE,
                control = VGAM::vglm.control(maxit = options()[["mixvlmc.maxit"]])
              )
            ),
            silent = TRUE
          )
          if (inherits(try_vglm, "try-error")) {
            err_cond <- as.character(attr(try_vglm, "condition"))
            err_to_ignore <- list(
              stringr::coll("vglm() only handles full-rank models (currently)"),
              stringr::coll("Error in if (min(ans) == 0 || max(ans) == 1)"),
              stringr::coll("there are NAs in eta in @inverse")
            )
            catched <- FALSE
            for (msg in err_to_ignore) {
              if (stringr::str_detect(err_cond, msg)) {
                ## "fake" result, interpreted as a low rank result
                result <- structure(list(coefficients = c(NA), ll = NA, rank = 0, target = NA, class = "constant_model"))
                catched <- TRUE
                break
              }
            }
            if (!catched) {
              stop(attr(try_vglm, "condition"))
            }
          }
        } else {
          try_vglm <- try(
            withCallingHandlers(
              warning = vgam_warning_ignore,
              result <- VGAM::vglm(target ~ 1,
                data = mm, family = VGAM::multinomial(refLevel = 1),
                x.arg = FALSE, y.arg = FALSE, model = FALSE,
                control = VGAM::vglm.control(maxit = options()[["mixvlmc.maxit"]])
              )
            ),
            silent = TRUE
          )
          if (inherits(try_vglm, "try-error")) {
            stop(attr(try_vglm, "condition"))
          }
        }
        if (inherits(result, "vglm")) {
          if (is_glm_low_rank(result)) {
            if (result@iter >= options()[["mixvlmc.maxit"]]) {
              message("vglm.fit did not converge for a discarded low rank model")
            }
          } else {
            if (result@iter >= options()[["mixvlmc.maxit"]]) {
              warning("vglm.fit did not converge")
            }
          }
        }
      }
      result
    } else if (engine == "multinom") {
      if (ncol(mm) > 0) {
        mm$target <- target
        try_nnet <- try(
          withCallingHandlers(
            warning = multinom_warning_ignore_generator(target, target_dist),
            result <- nnet::multinom(target ~ ., data = mm, trace = FALSE, maxit = options()[["mixvlmc.maxit"]])
          ),
          silent = TRUE
        )
      } else {
        try_nnet <- try(
          withCallingHandlers(
            warning = multinom_warning_ignore_generator(target, target_dist),
            result <- nnet::multinom(target ~ 1, trace = FALSE)
          ),
          silent = TRUE
        )
      }
      if (inherits(try_nnet, "try-error")) {
        stop(attr(try_nnet, "condition"))
      }
      result$rank <- length(stats::coef(result))
      result
    } else {
      ## should not happen
      NULL
    }
  }
}
