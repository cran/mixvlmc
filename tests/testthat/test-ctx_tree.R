test_that("the constructor keeps the state space", {
  vals <- c("1", "A", "Q")
  expect_identical(states(new_ctx_tree(vals)), vals)
})

test_that("the empty constructor returns a zero depth tree", {
  expect_identical(depth(new_ctx_tree(c(0, 1))), 0L)
})

test_that("the depth is calculated correctly", {
  expect_identical(depth(build_demo_tree(1:3, 4)), 4L)
  expect_identical(depth(build_demo_tree(1:4, 2)), 2L)
  expect_identical(depth(build_demo_tree(c("a", "b"), 3)), 3L)
})

test_that("the context number is calculted correctly", {
  expect_identical(context_number(build_demo_tree(1:3, 4)), as.integer(3^4))
  expect_identical(context_number(build_demo_tree(1:4, 2)), as.integer(4^2))
  expect_identical(context_number(build_demo_tree(c("a", "b"), 3)), as.integer(2^3))
})

test_that("the context tree has the correct class", {
  dts <- c(0, 1, 1, 1, 0, 0, 1, 0, 1, 0)
  dts_ctree <- ctx_tree(dts, min_size = 1, max_depth = 2)
  expect_s3_class(dts_ctree, "ctx_tree", exact = TRUE)
})

test_that("a warning is given when the state space is large", {
  dts <- rep(1:20, 2)
  expect_warning({
    dts_ctree <- ctx_tree(dts, min_size = 1, max_depth = 2)
  })
})

test_that("printing", {
  dts <- c(0, 1, 1, 1, 0, 0, 1, 0, 1, 0)
  dts_ctree <- ctx_tree(dts, min_size = 1, max_depth = 2)
  expect_snapshot_output(print(dts_ctree))
})
