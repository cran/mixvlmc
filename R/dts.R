## Currently the dts class is a trick to avoid printing the seed attribute
## of results of simulate.*

#' @export
print.dts <- function(x, ...) {
  attr(x, "seed") <- NULL
  NextMethod()
}
