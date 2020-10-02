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

# Plot data curves
plot.new( )
plot.window( xlim=c(-1,1), ylim=c(-1,1), asp = 1)
lapply(data_curves, lines, col = "gray")



library(elasdics)

fit_procrustes_mean <- function(srv_data_curves, knots, max_iter, type, eps)
{
  t_optims <- lapply(srv_data_curves, function(srv_data_curve) {
    c(srv_data_curve$t, 1)
  })
  coefs <- 0
  for (i in 1:max_iter) {
    model_data <- get_model_data(t_optims, srv_data_curves, knots, type)
    coefs_old <- coefs
    
    if (i == 2){
      # Build id column, super hacky way. This is bullshit (but it works)
      ids <- do.call(c, lapply(t_optims, function(x) length(x)-1))  #note: t_optimns has 1 more
      ids <- rbind(1:length(ids), ids)
      ids <- apply(ids, 2, function(x) rep(x[1], times = x[2]))
      ids <- do.call(c, ids)

      # BUILD RESPONSE AND DESIGN MATRIX
      model_data_complex <- complex(re = model_data[,2], im = model_data[,3]) %>% 
          matrix(nrow = dim(model_data)[1])
      model_data_complex <- data.frame(id = ids, m_long = model_data$m_long, q_m_long_complex = model_data_complex)
      
      print(model_data_complex)
      
    }
    
    # TAKE MODEL_DATA AND COMPUTE PROCRUSTES MEAN
    

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


#compute elastic means
knots <- seq(0,1, length = 11)
smooth_elastic_mean <- compute_elastic_mean(data_curves, knots = knots)
lines(get_evals(smooth_elastic_mean), type = "l", col = "red", lwd = 2)