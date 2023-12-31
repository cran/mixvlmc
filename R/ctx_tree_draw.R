#' Control parameters for `draw`
#'
#' This function returns a list used to fine tune the [draw()] function behaviour.
#'
#' @param root character used for the root node.
#' @param first_node characters used for the first child of a node.
#' @param next_node characters used for other children of a node.
#' @param vbranch characters used to represent a branch in a vertical way.
#' @param hbranch characters used to represent a branch in a horizontal was.
#' @param open_ct characters used to start each node specific text representation.
#' @param close_ct characters used to end each node specific text representation.
#'
#' @returns a list
#' @export
#'
#' @examples
#' draw_control(open_ct = "[", close_ct = "]")
draw_control <- function(root = "*",
                         first_node = "+",
                         next_node = "'",
                         vbranch = "|",
                         hbranch = "--",
                         open_ct = "(",
                         close_ct = ")") {
  list(
    root = root,
    first_node = first_node,
    next_node = next_node,
    vbranch = vbranch,
    hbranch = hbranch,
    open_ct = open_ct,
    close_ct = close_ct
  )
}

cat_with_prefix <- function(label, prefix, node_txt, control) {
  node_txt_lines <- unlist(stringr::str_split(node_txt, "\n"))
  cat(" ", control$open_ct, node_txt_lines[1], sep = "")
  if (length(node_txt_lines) > 1) {
    local_prefix <- stringr::str_pad(prefix, stringr::str_length(label) + 1 + stringr::str_length(control$open_ct), side = "right")
    for (k in seq_along(node_txt_lines)[-1]) {
      cat("\n", local_prefix, node_txt_lines[k], sep = "")
    }
  }
  cat(control$close_ct)
}

rec_draw <- function(label, prefix, ct, vals, control, node2txt, params) {
  cat(label)
  if (!is.null(node2txt)) {
    node_txt <- node2txt(ct, params)
    if (!is.null(node_txt)) {
      cat_with_prefix(label, prefix, node_txt, control)
    }
  }
  cat("\n")
  if (!is.null(ct$children)) {
    nst <- nb_sub_tree(ct)
    if (nst > 1) {
      c_symbol <- control$first_node
    } else {
      c_symbol <- control$next_node
    }
    idx <- 1
    for (v in seq_along(ct$children)) {
      child <- ct$children[[v]]
      if (length(child) > 0) {
        c_prelabel <- stringr::str_c(c_symbol, control$hbranch, " ")
        if (idx < nst) {
          c_prefix <- control$vbranch
        } else {
          c_prefix <- stringr::str_pad("", stringr::str_length(control$vbranch))
        }
        c_prefix <- stringr::str_pad(c_prefix, stringr::str_length(c_prelabel), side = "right")
        ## recursive call
        rec_draw(
          stringr::str_c(prefix, c_prelabel, vals[v]),
          stringr::str_c(prefix, c_prefix), child, vals, control, node2txt, params
        )
        ## prepare for next child
        c_symbol <- control$next_node
        idx <- idx + 1
      }
    }
  }
}

#' Text based representation of a context tree
#'
#' This function 'draws' a context tree as a text.
#'
#' The function uses basic "ascii art" to represent the context tree. Characters
#' used to represent the structure of the tree, e.g. branches, can be modified
#' using [draw_control()].
#'
#' In addition to the structure of the context tree, `draw` can represent
#' information attached to the node (contexts and partial contexts). This is
#' controlled by additional parameters depending on the type of the context
#' tree.
#'
#' @param ct a context tree.
#' @param control a list of low level control parameters of the text
#'   representation. See details and [draw_control()].
#' @param ... additional arguments for draw.
#' @returns the context tree (invisibly).
#' @examples
#' dts <- sample(c(0, 1), 100, replace = TRUE)
#' ctree <- ctx_tree(dts, min_size = 10, max_depth = 2)
#' draw(ctree)
#' dts_c <- sample(c("A", "B", "CD"), 100, replace = TRUE)
#' ctree_c <- ctx_tree(dts_c, min_size = 10, max_depth = 2)
#' draw(ctree_c, draw_control(root = "x"))
#' @export
draw <- function(ct, control = draw_control(), ...) {
  UseMethod("draw")
  invisible(ct)
}

ctx_tree_node2txt <- function(ct, params) {
  if (is.null(ct[["f_by"]])) {
    NULL
  } else {
    if (!is.null(params[["frequency"]])) {
      if (params$frequency == "detailed") {
        stringr::str_c(ct[["f_by"]], collapse = ",")
      } else if (params$frequency == "total") {
        as.character(sum(ct[["f_by"]]))
      } else {
        NULL
      }
    } else {
      NULL
    }
  }
}

#' @inherit draw
#' @param frequency this parameter controls the display of node level
#'   information in the tree. The default `NULL` value does not include
#'   anything. Setting `frequency` to `"total"` includes the frequency of the
#'   (partial) context of the node, while `"detailed"` includes the frequency of
#'   the states that follow the context (as in [contexts.ctx_tree()]).
#' @examples
#' dts_c <- sample(c("A", "B", "CD"), 100, replace = TRUE)
#' ctree_c <- ctx_tree(dts_c, min_size = 10, max_depth = 2)
#' draw(ctree_c, frequency = "total")
#' draw(ctree_c, frequency = "detailed")
#' @export
draw.ctx_tree <- function(ct, control = draw_control(), frequency = NULL, ...) {
  if (is.null(frequency)) {
    rec_draw(control$root, "", ct, ct$vals, control, NULL, list(...))
  } else {
    frequency <- match.arg(frequency, c("total", "detailed"))
    rec_draw(control$root, "", ct, ct$vals, control, ctx_tree_node2txt, c(list(frequency = frequency), list(...)))
  }
  invisible(ct)
}
