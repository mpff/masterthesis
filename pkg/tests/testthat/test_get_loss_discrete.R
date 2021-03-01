test_that("discrete_loss for example values",{
  data_curve <- data.frame(x1 = 1:6*sin(1:6), x2 = cos(1:6))
  srv_data <- get_srv_from_points(data_curve)
  expect_equal(get_loss_discrete(t = 0:5/5, srv_data, srv_data), 10.24,
               tolerance=1e-2)
  expect_equal(get_loss_discrete(t = c(0, rep(0.1, 3), 0.5, 1),
                                 srv_data, srv_data), 6.932, tolerance = 1e-2)
  expect_error(get_loss_discrete(t = srv_data$t, srv_data, srv_data))
})
