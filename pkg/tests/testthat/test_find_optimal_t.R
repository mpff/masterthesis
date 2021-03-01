test_that("discrete and smooth give same result if problem is nice enough",{
  data_curve1 <- as.data.frame(rbind(c(0,0), c(2,1), c(1,2), c(0,0)))
  data_curve1$t <- c(0, 1/3, 2/3, 1)
  srv_vec1 <- get_srv_from_points(data_curve1)
  srv_curve_1 <- function(t){
    sapply(t, function(t) t(srv_vec1)[-1,findInterval(t, c(srv_vec1$t,1), rightmost.closed = T)])
  }
  data_curve2 <- as.data.frame(rbind(c(0,0), c(2,0), c(0.5,0.5), c(0,2), c(0,0)))
  srv_vec2 <- get_srv_from_points(data_curve2)

  expect_equal(as.numeric(find_optimal_t_discrete(c(srv_vec1$t,1), t(srv_vec1[,-1]),
                                       c(srv_vec2$t,1), t(srv_vec2[,-1]))),
               as.numeric(find_optimal_t(srv_curve_1, c(srv_vec2$t,1), t(srv_vec2[,-1]),
                              eps = 0.1)),
               tolerance=1e-3)
})

test_that("grad is null at maximiser",{
  srv_curve <- function(t){
    sapply(t, function(t) t(c(sin(5*t), t*cos(5*t))))
  }
  data_curve <- as.data.frame(rbind(c(0,0), c(2,0), c(-1,-1)))
  srv_vec <- get_srv_from_points(data_curve)
  t_optim <- find_optimal_t(srv_curve, c(srv_vec$t,1), t(srv_vec[,-1]))[[2]]
  expect_equal(get_grad(t_optim, srv_curve, c(srv_vec$t,1), t(srv_vec[,-1])), 0,
               tolerance=1e-3)
})

test_that("symmetric shape",{
  srv_curve <- function(t){
    sapply(t, function(t) t(c(cos(pi*t), -sin(pi*t))))
  }

  data_curve <- as.data.frame(rbind(c(0,1), c(1,0), c(0,-1)))
  data_curve$t <- c(0, 0.3, 1)
  srv_vec1 <- get_srv_from_points(data_curve)
  t_optim1 <- find_optimal_t(srv_curve, c(srv_vec1$t,1), t(srv_vec1[,-1]))

  data_curve$t <- c(0, 0.9, 1)
  srv_vec2 <- get_srv_from_points(data_curve)
  t_optim2 <- find_optimal_t(srv_curve, c(srv_vec2$t,1), t(srv_vec2[,-1]))
  expect_equal(as.numeric(t_optim1), c(0, 0.5, 1), tolerance=1e-5)
  expect_equal(as.numeric(t_optim1), as.numeric(t_optim2), tolerance=1e-5)
})
