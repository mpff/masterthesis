test_that("plot runs without warning",{
  curve <- function(t){
    rbind(t*cos(13*t), t*sin(13*t))
  }
  set.seed(18)
  data_curves <- lapply(1:4, function(i){
    m <- sample(10:15, 1)
    delta <- abs(rnorm(m, mean = 1, sd = 0.05))
    t <- cumsum(delta)/sum(delta)
    data.frame(t(curve(t)) + 0.07*t*matrix(cumsum(rnorm(2*length(delta))),
                                           ncol = 2))
  })

  #compute elastic means
  knots <- seq(0,1, length = 11)
  smooth_elastic_mean <- compute_elastic_mean(data_curves, knots = knots)
  expect_warning(plot(smooth_elastic_mean), regexp = NA)
})


test_that("plot gives error if more than two dim",{
  data_curve1 <- data.frame(x1 = sin(1:7/4*pi), x2 = cos(1:7/4*pi),
                            x3 = -sin(1:7/4*pi))
  data_curve2 <- data_curve <- data.frame(x1 = sin(1:15/8*pi), x2 = cos(1:15/8*pi),
                                          x3 = 1:15/8*pi)
  data_curves <- list(data_curve1, data_curve2)

  #compute elastic means
  knots <- seq(0,1, length = 11)
  elastic_mean <- compute_elastic_mean(data_curves, knots = knots, type = "polygon")
  expect_error(plot(elastic_mean), "Plotting option only for planar curves")
})
