library(elasdics)
library(mgcv)
library(dplyr)
library(shapeboost)


# MODIFY FIT_MEAN FUNCTION!
fit_procrustes_mean <- function(srv_data_curves, knots, max_iter, type, eps)
{
  t_optims <- lapply(srv_data_curves, function(srv_data_curve) {
    c(srv_data_curve$t, 1)
  })
  coefs <- 0
  
  # NEW: requirements (as of now)
  require(mgcv)
  require(dplyr)
  require(shapeboost)
  
  # NEW: Build id column, super hacky way. THIS NEEDS REWORK!
  ids <- do.call(c, lapply(t_optims, function(x) length(x)-1))  #note: t_optims has 1 more
  ids <- rbind(1:length(ids), ids)
  ids <- apply(ids, 2, function(x) rep(x[1], times = x[2]))
  ids <- do.call(c, as.list(ids))  # coerce as.list for curves of equal length(t)
  
  # NEW: Build arg grid for evaluation of smoothed cov.
  arg.grid <- seq(0,1, len = 200)
  
  for (i in 1:max_iter) {
    model_data <- get_model_data(t_optims, srv_data_curves, knots, type)
    coefs_old <- coefs
    
    if (i == 2){
      # NEW: Fit Procrustes Mean
      # x,y to complex, add id column
      model_data_complex <- complex(re = model_data[,2], im = model_data[,3]) %>% 
          matrix(nrow = dim(model_data)[1])
      model_data_complex <- data.frame(id = ids, m_long = model_data$m_long, q_m_long = model_data_complex)
      
      # build response on s,t-grid per curve
      cov_dat <- lapply(split(model_data_complex, model_data_complex$id), function(x) {
        combs <- combn(1:nrow(x), 2)
        data.frame(
          qq = x$q_m_long[combs[1,]] * Conj(x$q_m_long[combs[2,]]),
          t = x$m_long[combs[1,]],
          s = x$m_long[combs[2,]]
        )
      })
      cov_dat <- do.call(rbind, cov_dat)
      
      # fit covariance surface. ToDo: Check these params: k, knots, cyclic !
      cov_fit_re <- bam( Re(qq) ~ s(t, s, bs = "sps", k = length(knots)),
                         data = cov_dat )
      cov_fit_im <- bam( Im(qq) ~ -1 + s(t, s, bs = "sps", k = length(knots), xt = list(skew = TRUE)),
                         data = cov_dat )
      
      # predict smoothed covariance surface on arg.grid
      cov_dat <- expand.grid(t = arg.grid, s = arg.grid)
      yy <- matrix( complex(
          real = predict(cov_fit_re, newdata = cov_dat),
          imaginary = predict(cov_fit_im, newdata = cov_dat) ),
        ncol = length(arg.grid) )
      
      # calculate procrustes mean as leading eigenvector
      ei <- eigen(yy)
      mu <- data.frame(t = arg.grid, X1 = Re(ei$vectors[,1]), X2 = Im(ei$vectors[,1]))
      
      # calculate procrustes fit of srv data curves
      ## IS THIS NECESSARY?
      
      # Diagnostic plots.
      plot.new( )
      plot.window( xlim=c(-0.002,0.002), ylim=c(-0.002,0.002), asp = 1)
      mu_test <- get_points_from_srv(mu)
      lines(mu_test, type = "l", col = "blue", lwd = 2)
    }
    
    coefs <- apply(model_data[, -1], 2, function(q_m_x_long) {
      q_m_x_long[!is.finite(q_m_x_long)] <- NA
      coef(lm(q_m_x_long ~ -1 + make_design(model_data[,1], knots = knots, closed = FALSE, type = type)))
    })
    
  
    # TAKE MEAN AND CALCULATE THE APPROPRIATE COEFFICIENTS
    
    stop_crit <- sum((coefs - coefs_old)^2)/sum(coefs^2)
    if (stop_crit < eps | max_iter == 0) {
      rownames(coefs) <- NULL
      colnames(coefs) <- colnames(srv_data_curves[[1]][,-1])
      return(list(type = type, coefs = coefs, knots = knots, 
                  t_optims = t_optims))
    }
    if (type == "smooth") {
      pfun <- function(t) {
        t(make_design(t, knots = knots, closed = FALSE, 
                      type = type) %*% coefs)
      }
      t_optims <- lapply(1:length(srv_data_curves), function(j) {
        t_optim <- find_optimal_t(srv_curve = pfun, s = c(srv_data_curves[[j]]$t, 
                                                          1), q = t(srv_data_curves[[j]][, -1]), initial_t = t_optims[[j]], 
                                  eps = eps * 100/i)
        attr(t_optim, "dist_to_mean") <- attr(t_optim, 
                                              "dist")
        attr(t_optim, "dist") <- NULL
        t_optim
      })
    }
    else {
      t_optims <- lapply(1:length(srv_data_curves), function(j) {
        t_optim <- find_optimal_t_discrete(r = knots, 
                                           p = t(coefs), s = c(srv_data_curves[[j]]$t, 
                                                               1), q = t(srv_data_curves[[j]][, -1]), initial_t = t_optims[[j]])
        attr(t_optim, "dist_to_mean") <- attr(t_optim, 
                                              "dist")
        attr(t_optim, "dist") <- NULL
        t_optim
      })
    }
  }
  warning("Stopping criteria eps has not been reached! Consider more iterations max_iter")
  rownames(coefs) <- NULL
  colnames(coefs) <- colnames(srv_data_curves[[1]][, -1])

  return(list(type = type, coefs = coefs, knots = knots, t_optims = t_optims))
}

environment(fit_procrustes_mean) <- asNamespace('elasdics')
assignInNamespace("fit_mean", fit_procrustes_mean, ns = "elasdics")



#########################################
### TESTS ON DIGIT3 AND SIMULATED SPIRALS

set.seed(18)

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
#data_curves <- lapply(data_curves, rand.rotate)
#data_curves <- lapply(data_curves, rand.scale)
data_curves <- lapply(data_curves, center_curve)


### SIMULATED SPIRALS: compute elastic means
knots <- seq(0,1, length = 11)
smooth_elastic_mean <- compute_elastic_mean(data_curves, knots = knots)

#plot result
plot.new( )
plot.window( xlim=c(-1,1), ylim=c(-1,1), asp = 1)
lapply(data_curves, lines, col = "gray")
lines(get_evals(smooth_elastic_mean), type = "l", col = "red", lwd = 2)


### DIGITS 3: compute elastic means
knots <- seq(0,1, length = 11)
smooth_elastic_mean <- compute_elastic_mean(d3_curves, knots = knots)

#plot result
plot.new( )
plot.window( xlim=c(-15,15), ylim=c(-15,15), asp = 1)
lapply(d3_curves, lines, col = "gray")
lines(get_evals(smooth_elastic_mean), type = "l", col = "red", lwd = 2)