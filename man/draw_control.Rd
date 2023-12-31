% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ctx_tree_draw.R
\name{draw_control}
\alias{draw_control}
\title{Control parameters for \code{draw}}
\usage{
draw_control(
  root = "*",
  first_node = "+",
  next_node = "'",
  vbranch = "|",
  hbranch = "--",
  open_ct = "(",
  close_ct = ")"
)
}
\arguments{
\item{root}{character used for the root node.}

\item{first_node}{characters used for the first child of a node.}

\item{next_node}{characters used for other children of a node.}

\item{vbranch}{characters used to represent a branch in a vertical way.}

\item{hbranch}{characters used to represent a branch in a horizontal was.}

\item{open_ct}{characters used to start each node specific text representation.}

\item{close_ct}{characters used to end each node specific text representation.}
}
\value{
a list
}
\description{
This function returns a list used to fine tune the \code{\link[=draw]{draw()}} function behaviour.
}
\examples{
draw_control(open_ct = "[", close_ct = "]")
}
