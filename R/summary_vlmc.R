#' @exportS3Method
summary.vlmc <- function(object, ...) {
  ctx <- contexts(object, frequency = "detailed")
  all_length <- lengths(ctx$context)
  res <- list(
    state_space = states(object),
    nb = nrow(ctx),
    depth = max(all_length),
    avg_depth = mean(all_length),
    cutoff = object$cutoff,
    alpha = object$alpha
  )
  structure(res, class = "summary.vlmc")
}

print_summary_vlmc <- function(x) {
  cat(paste(" Cutoff: ", signif(x$cutoff, 4), " (quantile: ", signif(x$alpha, 4), ")\n", sep = ""))
  cat(paste(" Number of contexts:", x$nb, "\n"))
  cat(paste(" Maximum context length:", x$depth, "\n"))
  cat(paste(" Average context length:", signif(x$avg_depth, 4), "\n"))
}

#' @exportS3Method
print.summary.vlmc <- function(x, ...) {
  cat(paste(
    "VLMC context tree on",
    paste(x$state_space, collapse = ", ")
  ), "\n")
  print_summary_vlmc(x)
  invisible(x)
}
