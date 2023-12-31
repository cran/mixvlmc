## a simple function that calls draw with different parameters
exercise_draw <- function(ctx) {
  draw(ctx)
  draw(ctx, frequency = "total")
  draw(ctx, frequency = "detailed")
  draw(ctx, control = draw_control(
    root = "x",
    open_ct = "[", close_ct = "]",
    first_node = "*",
    next_node = "°",
    vbranch = "^",
    hbranch = "->"
  ), frequency = "detailed")
}
