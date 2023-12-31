% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ctx_tree_draw.R
\name{draw}
\alias{draw}
\title{Text based representation of a context tree}
\usage{
draw(ct, control = draw_control(), ...)
}
\arguments{
\item{ct}{a context tree.}

\item{control}{a list of low level control parameters of the text
representation. See details and \code{\link[=draw_control]{draw_control()}}.}

\item{...}{additional arguments for draw.}
}
\value{
the context tree (invisibly).
}
\description{
This function 'draws' a context tree as a text.
}
\details{
The function uses basic "ascii art" to represent the context tree. Characters
used to represent the structure of the tree, e.g. branches, can be modified
using \code{\link[=draw_control]{draw_control()}}.

In addition to the structure of the context tree, \code{draw} can represent
information attached to the node (contexts and partial contexts). This is
controlled by additional parameters depending on the type of the context
tree.
}
\examples{
dts <- sample(c(0, 1), 100, replace = TRUE)
ctree <- ctx_tree(dts, min_size = 10, max_depth = 2)
draw(ctree)
dts_c <- sample(c("A", "B", "CD"), 100, replace = TRUE)
ctree_c <- ctx_tree(dts_c, min_size = 10, max_depth = 2)
draw(ctree_c, draw_control(root = "x"))
}
