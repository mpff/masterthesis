test_that("get arc length parametrisation points on circle", {
  data_curve <- data.frame(x1 = sin(1:6), x2 = cos(1:6))
  expect_equal(get_arc_length_param(data_curve), 0:5/5)
})

test_that("srv_from_points und point_from_srv are inverse",{
  data_curve <- data.frame(x1 = 1:6*sin(1:6), x2 = cos(1:6))
  data_curve_new <- get_points_from_srv(get_srv_from_points(data_curve))
  expect_equal(apply(data_curve_new, 2, diff),
               apply(data_curve, 2, diff))
  data_curve$t <- 0:5/5
  data_curve_new2 <- get_points_from_srv(get_srv_from_points(data_curve))
  expect_equal(apply(data_curve_new, 2, diff),
               apply(data_curve_new2, 2, diff))
})

test_that("input checking", {
  data_curve <- data.frame(x1 = sin(1:6), x2 = cos(1:6), t = 1:6)
  expect_error(get_srv_from_points(data_curve),
               "Parametrisation must be starting at 0 and ending at 1")
  data_curve$t <- c(0,0.2,0.4,0.6,0.5,1)
  expect_error(get_srv_from_points(data_curve),
               "Parametrisation needs to be strictly increasing")
})
