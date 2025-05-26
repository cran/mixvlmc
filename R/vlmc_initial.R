rec_no_context <- function(ct, ctx) {
  if (is.null(ct) || is.null(ct$children) || length(ct$children) == 0) {
    ## a leaf is a context
    NULL
  } else {
    if (nb_sub_tree(ct) < length(ct$children)) {
      ## this is a context
      local <- NULL
    } else {
      ## full node, not a context
      local <- list(ctx)
    }
    for (k in seq_along(ct$children)) {
      pre <- rec_no_context(ct$children[[k]], c(ctx, k))
      local <- c(local, pre)
    }
    local
  }
}

# vlmc is here a multi_vlmc
initial_distribution <- function(vlmc) {
  one_id <- function(dts, nb_vals, pos, current) {
    val <- dts[pos] + 1L
    if (is.null(current[["counts"]])) {
      current[["counts"]] <- rep(0L, nb_vals)
    }
    current$counts[val] <- current$counts[val] + 1L
    pos <- pos + 1L
    if (pos <= length(dts)) {
      if (is.null(current[["next"]])) {
        current[["next"]] <- vector(length = nb_vals, mode = "list")
      }
      current[["next"]][[val]] <- one_id(dts, nb_vals, pos, current[["next"]][[val]])
    }
    current
  }
  the_probs <- list()
  for (k in seq_along(vlmc$ix)) {
    ## process one series
    dts <- vlmc$ix[[k]]
    the_probs <- one_id(dts, length(vlmc$vals), 1L, the_probs)
  }
  the_probs
}

nb_prob <- function(the_prob) {
  if (is.null(the_prob)) {
    0L
  } else {
    res <- 1L
    sub_probs <- the_prob[["next"]]
    if (!is.null(sub_probs)) {
      for (k in seq_along(sub_probs)) {
        if (!is.null(sub_probs[[k]])) {
          res <- res + nb_prob(sub_probs[[k]])
        }
      }
    }
    res
  }
}

flatten_prob <- function(init_prob) {
  flatten_rec <- function(ctx, the_prob) {
    if (is.null(the_prob)) {
      NULL
    } else {
      pre_res <- data.frame(
        ctx = I(list(ctx)),
        counts = I(list(the_prob$counts))
      )
      sub_probs <- the_prob[["next"]]
      if (!is.null(sub_probs)) {
        for (k in seq_along(sub_probs)) {
          if (!is.null(sub_probs[[k]])) {
            sub_res <- flatten_rec(c(ctx, k - 1L), sub_probs[[k]])
            pre_res <- rbind(pre_res, sub_res)
          }
        }
      }
      pre_res
    }
  }
  flatten_rec(NULL, init_prob)
}

flatten_prob_no_ctx <- function(init_prob, model) {
  flatten_rec <- function(ctx, the_prob, model) {
    if (is.null(the_prob)) {
      NULL
    } else {
      ## should be optimised
      prob_ctx_node <- find_sequence(model, states(model)[ctx + 1L])
      if (is.null(prob_ctx_node)) {
        pre_res <- NULL
      } else {
        ct_children <- children(prob_ctx_node)
        if (length(ct_children) > 0) {
          pre_res <- data.frame(
            ctx = I(list(ctx)),
            counts = I(list(the_prob$counts))
          )
        } else {
          pre_res <- NULL
        }
      }
      sub_probs <- the_prob[["next"]]
      if (!is.null(sub_probs)) {
        for (k in seq_along(sub_probs)) {
          if (!is.null(sub_probs[[k]])) {
            sub_res <- flatten_rec(c(ctx, k - 1L), sub_probs[[k]], model)
            if (!is.null(pre_res)) {
              pre_res <- rbind(pre_res, sub_res)
            } else {
              pre_res <- sub_res
            }
          }
        }
      }
      pre_res
    }
  }
  flatten_rec(NULL, init_prob, model)
}

prune_ctx_prob <- function(init_prob, model) {
  prune_rec <- function(ctx, the_prob, model) {
    if (is.null(the_prob)) {
      NULL
    } else {
      ## should be optimised
      prob_ctx_node <- find_sequence(model, states(model)[ctx + 1L])
      if (is.null(prob_ctx_node)) {
        pre_res <- NULL
      } else {
        ct_children <- children(prob_ctx_node)
        if (length(ct_children) > 0) {
          pre_res <- list(counts = the_prob$counts)
        } else {
          pre_res <- NULL
        }
      }
      sub_probs <- the_prob[["next"]]
      if (!is.null(sub_probs)) {
        pre_next <- vector(mode = "list", length = length(sub_probs))
        no_next <- TRUE
        for (k in seq_along(sub_probs)) {
          if (!is.null(sub_probs[[k]])) {
            pruned_next <- prune_rec(c(ctx, k - 1L), sub_probs[[k]], model)
            if (!is.null(pruned_next)) {
              pre_next[[k]] <- pruned_next
              no_next <- FALSE
            }
          }
        }
        if (!no_next) {
          if (is.null(pre_res)) {
            pre_res <- list("next" = pre_next)
          } else {
            pre_res[["next"]] <- pre_next
          }
        }
      }
      pre_res
    }
  }
  prune_rec(NULL, init_prob, model)
}

loglikelihood_prob <- function(init_prob) {
  log_rec <- function(the_prob) {
    if (is.null(the_prob)) {
      0
    } else {
      ll <- local_loglikelihood_vlmc(the_prob$counts)
      if (!is.null(the_prob[["next"]])) {
        active_subs <- the_prob[["next"]]
        sub_lls <- sapply(active_subs, log_rec)
        ll <- ll + sum(sub_lls)
      }
      ll
    }
  }
  log_rec(init_prob)
}

loglikelihood_mixed <- function(ix, init_prob, model) {

}

internal_nodes <- function(model) {
  internal_rec <- function(ctx, model) {
    if (is.null(model[["children"]])) {
      NULL
    } else {
      res <- list(ctx)
      for (k in seq_along(model[["children"]])) {
        pre_res <- internal_rec(c(k - 1L, ctx), model[["children"]][[k]])
        res <- c(res, pre_res)
      }
      res
    }
  }
  internal_rec(c(), model)
}

fit_initial_prob <- function(dts_ix, ctxs, nb_vals) {
  the_probs <- vector(mode = "list", length = length(ctxs))
  for (k in seq_along(the_probs)) {
    the_ctx <- ctxs[[k]]
    counts <- rep(0, nb_vals)
    for (l in seq_along(dts_ix)) {
      the_dts <- dts_ix[[l]]
      if (length(the_dts) <= length(the_ctx)) {
        break
      }
      if (length(the_ctx) == 0) {
        the_val <- the_dts[1] + 1L
        counts[the_val] <- counts[the_val] + 1L
      } else {
        if (identical(the_dts[1:length(the_ctx)], the_ctx)) {
          the_val <- the_dts[length(the_ctx) + 1] + 1L
          counts[the_val] <- counts[the_val] + 1L
        }
      }
    }
    the_probs[[k]] <- counts
  }
  data.frame(ctx = I(ctxs), counts = I(the_probs))
}
