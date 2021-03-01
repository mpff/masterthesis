test_that("get evals default",{
  curve <- function(t){
    c(sin(2*pi*t), cos(2*pi*t))
  }
  expect_equal(as.numeric(get_evals(curve, seq(0,1,0.1))[11,]), c(0,1))
  expect_equal(get_evals(curve, seq(0,1,0.01)), get_evals(curve))
})

test_that("get evals data.frame",{
  data_frame <- data.frame(s = c(0, 0.3, 0.6, 1),
                           x1 = c(1, 0.5, -1, -1),
                           x2 = c(1, -0.5, -1, 1))
  expect_error(get_evals(data_frame), "Parametrisation t must be in the first column!")
  names(data_frame)[1] <- "t"
  expect_equal(nrow(get_evals(data_frame)), 101)
})


test_that("get evals elastic mean",{
  data_curves <- list(data.frame(x1 = c(1, 0.5, -1, -1), x2 = c(1, -0.5, -1, 1)),
                      data.frame(x1 = c(0.1,0.7)*sin(1:6), x2 = cos(1:6)))
  elastic_mean <- compute_elastic_mean(data_curves, max_iter = 10)
  expect_equal(colnames(get_evals(elastic_mean)), c("x1","x2"))

  elastic_mean <- compute_elastic_mean(data_curves, type = "polygon")
  expect_equal(nrow(get_evals(elastic_mean)), 5)
})

test_that("get evals elastic mean with given t_grid",{
  data_curves <- list(data.frame(x1 = c(1, 0.5, -1, -1), x2 = c(1, -0.5, -1, 1)),
                      data.frame(x1 = c(0.1,0.7)*sin(1:6), x2 = cos(1:6)))
  elastic_mean <- compute_elastic_mean(data_curves, max_iter = 10)
  expect_equal(nrow(get_evals(elastic_mean, t_grid = 0:20/20)), 21)

  elastic_mean <- compute_elastic_mean(data_curves, type = "polygon")
  expect_equal(nrow(get_evals(elastic_mean, t_grid = 0:10/10)), 11)
})
