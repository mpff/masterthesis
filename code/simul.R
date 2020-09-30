library(elasdics)
 
set.seed(18)
 
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
data_curves <- lapply(data_curves, center_curve)
  
# Plot data curves
plot.new( )
plot.window( xlim=c(-1,1), ylim=c(-1,1), asp = 1)
lapply(data_curves, lines, col = "gray")
 
#compute elastic means
knots <- seq(0,1, length = 11)
smooth_elastic_mean <- compute_elastic_mean(data_curves, knots = knots)
lines(get_evals(smooth_elastic_mean), type = "l", col = "red", lwd = 2)
 
knots <- seq(0,1, length = 15)
polygon_elastic_mean <- compute_elastic_mean(data_curves, knots = knots, type = "poly")
lines(get_evals(polygon_elastic_mean), type = "l", col = "blue", lwd = 2)
 
 
### NEW STUFF HERE
# Rotate and scale curves
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
data_curves_new <- lapply(data_curves, rand.rotate)
data_curves_new <- lapply(data_curves_new, rand.scale)

# Plot data curves
plot.new( )
plot.window( xlim=c(-1,1), ylim=c(-1,1), asp = 1)
lapply(data_curves_new, lines, col = "gray")

#compute elastic means
knots <- seq(0,1, length = 11)
smooth_elastic_mean <- compute_elastic_mean(data_curves_new, knots = knots)
lines(get_evals(smooth_elastic_mean), type = "l", col = "red", lwd = 2)
 
knots <- seq(0,1, length = 15)
polygon_elastic_mean <- compute_elastic_mean(data_curves_new, knots = knots, type = "poly")
lines(get_evals(polygon_elastic_mean), type = "l", col = "blue", lwd = 2)