source("/home/mnl/Statistik/masterthesis/code/modify_compute_elastic_mean.R")

#########################################
### TESTS ON DIGIT3 AND SIMULATED SPIRALS

set.seed(18)

###########################################
# Functions for rotating and scaling curves

rand.rotate <- function(x){
  # rotate dataframe of 2D vectors randomly
  theta <- 2*pi*runif(1)
  mat <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2)
  x.rot <- as.matrix(x) %*% t(mat)
  as.data.frame(x.rot)
}
rand.scale <- function(x){
  # scale dataframe of 2D vectors randomly
  beta <- 0.5 + 0.5*runif(1)
  beta * x
}


# Load digits3 dataset
d3_curves <- shapes::digit3.dat
d3_curves <- apply(d3, MARGIN = 3, FUN = function(i){
  data.frame(X1 = i[,1], X2 = i[,2])
})
d3_curves <- lapply(d3_curves, center_curve)

# Simulate open planar curves
curve <- function(t){
  rbind(t*cos(13*t), t*sin(13*t))
}
data_curves <- lapply(1:4, function(i){
  m <- sample(10:15, 1)
  delta <- abs(rnorm(m, mean = 1, sd = 0.05))
  t <- cumsum(delta)/sum(delta)
  data.frame(t(curve(t)) + 0.07*t*matrix(cumsum(rnorm(2*length(delta))), ncol = 2))
})


data_curves <- lapply(data_curves, rand.rotate)
data_curves <- lapply(data_curves, rand.scale)
data_curves <- lapply(data_curves, center_curve)



### SIMULATED SPIRALS: compute elastic means
knots <- seq(0,1, length = 11)
smooth_elastic_mean <- compute_elastic_mean(data_curves, knots = knots)

#plot result
plot.new( )
plot.window( xlim=c(-1,1), ylim=c(-1,1), asp = 1)
lapply(data_curves, lines, col = "gray")
lines(get_evals(smooth_elastic_mean), type = "l", col = "red", lwd = 2)


if (TRUE) {
  ### DIGITS 3: compute elastic means
  knots <- seq(0,1, length = 11)
  smooth_elastic_mean <- compute_elastic_mean(d3_curves, knots = knots)
  
  #plot result
  plot.new( )
  plot.window( xlim=c(-15,15), ylim=c(-15,15), asp = 1)
  lapply(d3_curves, lines, col = "gray")
  lines(get_evals(smooth_elastic_mean), type = "l", col = "red", lwd = 2)
}