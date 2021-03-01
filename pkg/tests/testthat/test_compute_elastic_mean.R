test_that("initial param does not matter too much",{
  data_curve1 <- data.frame(x1 = sin(1:7/4*pi), x2 = cos(1:7/4*pi))
  data_curve2 <- data_curve <- data.frame(x1 = sin(1:15/8*pi), x2 = cos(1:15/8*pi))
  data_curves <- list(data_curve1, data_curve2)
  mean1 <- compute_elastic_mean(data_curves)
  mean2 <- compute_elastic_mean(data_curves, type = "poly", eps = 10^-5,
                                max_iter = 10)

  data_curves[[1]]$t <- c(0, 0.2, 0.3, 0.4, 0.6, 0.8, 1)
  mean3 <- compute_elastic_mean(data_curves, eps = 10^-5, max_iter = 10)
  mean4 <- compute_elastic_mean(data_curves, type = "poly")

  expect_equal(mean1$coefs, mean3$coefs, tolerance=1e-1)
  expect_equal(mean2$coefs, mean4$coefs, tolerance=1e-1)
})

test_that("0 and 1 iterations of open and closed means",{
  data_curve1 <- data.frame(x1 = sin(1:7/4*pi), x2 = cos(1:7/4*pi))
  data_curve2 <- data_curve <- data.frame(x1 = sin(1:15/8*pi), x2 = cos(1:15/8*pi))
  data_curves <- list(data_curve1, data_curve2)
  expect_warning(compute_elastic_mean(data_curves, max_iter = 1))
  expect_warning(compute_elastic_mean(data_curves, max_iter = 1, eps = 0.0001, closed = TRUE))
  expect_warning(compute_elastic_mean(data_curves, closed = TRUE, max_iter = 0),
                 regexp = NA)
})


test_that("polygon closed mean is closed",{
  data_curve1 <- data.frame(x1 = sin(0:7/4*pi), x2 = cos(0:7/4*pi))
  data_curve2 <- data_curve <- data.frame(x1 = sin(0:15/8*pi), x2 = cos(0:15/8*pi))
  data_curves <- list(data_curve1, data_curve2)
  elastic_mean <- compute_elastic_mean(data_curves, closed = TRUE,
                                       type = "polygon", knots = seq(0,1,0.1))
  points <- get_evals(elastic_mean)
  expect_equal(as.numeric(points[1,]), as.numeric(points[nrow(points),]), tolerance = 1e-2)
})


