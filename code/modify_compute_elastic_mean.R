###################################################
# Add proc2d parameter to 'compute_elastic_mean'. #
#   if true: use 'fit_mean_proc2d'                #
###################################################
compute_elastic_mean <- function(data_curves, 
                                        knots = seq(0, 1, len = 5),
                                        type = c("smooth", "polygon"), 
                                        closed = FALSE,
                                        proc2d = FALSE,
                                        eps = 0.01,
                                        pen_factor = 100,
                                        max_iter = 50)
{

  lapply(data_curves, function(data_curve) {
    if ("t" %in% names(data_curve)) 
      check_param(data_curve, closed)
  })
  
  if (closed) {
    data_curves <- lapply(data_curves, function(data_curve) {
      check_closed(data_curve)
    })
  }
  
  # Create arc length parametrization if curves don't have 't' column.
  data_curves <- lapply(data_curves, function(data_curve) {
    if (!("t" %in% colnames(data_curve))) {
      data.frame(t = get_arc_length_param(data_curve), 
                 data_curve)
    }
    else {
      param <- data_curve$t
      data_curve$t <- NULL
      data.frame(t = param, data_curve)
    }
  })
  
  # Get SRV data curves
  srv_data_curves <- lapply(data_curves, get_srv_from_points)
  
  type <- match.arg(type, c("smooth", "polygon"))
  
  # Calculate Elastic Mean!
  if (!closed) {
    # open curves.
    if (!proc2d){
      # This is new.
      elastic_mean <- fit_mean(srv_data_curves = srv_data_curves, 
                               knots = knots, max_iter = max_iter, type = type, 
                               eps = eps)
    }
    else if (proc2d){
      elastic_mean <- fit_mean_proc2d(srv_data_curves = srv_data_curves, 
                                      knots = knots, max_iter = max_iter, type = type, 
                                    eps = eps)
    }
    # This just prepares output
    data_curves <- lapply(1:length(data_curves), function(j) {
      data_curves[[j]]$t_optim <- elastic_mean$t_optims[[j]]
      attributes(data_curves[[j]]$t_optim) <- NULL
      data_curve <- data_curves[[j]][, c(1, ncol(data_curves[[j]]), 
                                         2:(ncol(data_curves[[j]]) - 1))]
      attr(data_curve, "dist_to_mean") <- attr(elastic_mean$t_optims[[j]], 
                                               "dist_to_mean")
      data_curve
    })
  }
  else if (closed) {
    if (proc2d) {
        warning("Option 'proc2d' does not work with closed curves. Calculating normal mean.")
    }
    # closed curves.
    elastic_mean <- fit_mean_closed(srv_data_curves = srv_data_curves, 
                                    knots = knots, max_iter = max_iter, type = type, 
                                    eps = eps, pen_factor = pen_factor)
    
    # This just prepares output I think. What's the additional stuff????
    data_curves <- lapply(1:length(data_curves), function(j) {
      t_optim <- elastic_mean$t_optims[[j]][-length(elastic_mean$t_optims[[j]])]
      data_curve <- data_curves[[j]][-nrow(data_curves[[j]]), 
                                     ]
      part2_idx <- 1:(length(t_optim) - elastic_mean$shift_idxs[j] + 
                        1)
      data_curve$t_optim <- c(t_optim[-part2_idx], t_optim[part2_idx])
      data_curve <- data_curve[, c(1, ncol(data_curve), 
                                   2:(ncol(data_curve) - 1))]
      data_curve <- rbind(data_curve, data_curve[1, ])
      data_curve$t[nrow(data_curve)] <- 1
      attr(data_curve, "dist_to_mean") <- attr(elastic_mean$t_optims[[j]], 
                                               "dist_to_mean")
      data_curve
    })
  }
  
  # Prepare final output
  elastic_mean$data_curves <- data_curves
  elastic_mean$closed <- closed
  elastic_mean$shift_idxs <- NULL
  elastic_mean$t_optims <- NULL
  if (proc2d){
      # This is new.
      elastic_mean$procrustes_curves <- NULL
      elastic_mean$G_optims <- NULL  # Rotation
      elastic_mean$b_optims <- NULL  # Scale
  }
  class(elastic_mean) <- "elastic_mean"
  elastic_mean
}





############################
# Define 'fit_mean_proc2d' #
############################
fit_mean_proc2d <- function(srv_data_curves, knots, max_iter, type, eps)
{
    
  # NEW: requirements (as of now)
  require(mgcv)
  require(dplyr)
  require(shapeboost)
    
    
  t_optims <- lapply(srv_data_curves, function(srv_data_curve) {
    c(srv_data_curve$t, 1)
  })
  coefs <- 0
  

  
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
    procrustes_fits <- lapply(srv_data_curves, function(x) {
      mu_eval_x1 <- splinefun(x=mu$t, y=mu$X1, method="natural",  ties = mean)(x$t)
      mu_eval_x2 <- splinefun(x=mu$t, y=mu$X2, method="natural",  ties = mean)(x$t)
      mu_eval <- complex(real = mu_eval_x1, imaginary = mu_eval_x2)
      x$complex <- complex(real = x[,2], imaginary = x[,3])
      qq <- Conj(x$complex) %*% x$complex
      qm <- Conj(x$complex) %*% mu_eval
      x$complex <- as.vector(qm) * x$complex / as.vector(qq)
      data.frame(t = x$t, X1 = Re(x$complex), X2 = Im(x$complex))
    })
    
    # Diagnostic plots.
    #plot.new( )
    #plot.window( xlim=c(-0.002,0.002), ylim=c(-0.002,0.002), asp = 1)
    #mu_test <- get_points_from_srv(mu)
    #data_curves_procrustes <- lapply(procrustes_fits, get_points_from_srv)
    #lapply(data_curves_procrustes, lines, col = "gray")
    #lines(mu_test, type = "l", col = "blue", lwd = 2)
    
    # Prepare for warping (coefs + srv_data_curves)
    mu_knots_x1 <- splinefun(x=mu$t, y=mu$X1, method="natural",  ties = mean)(knots)
    mu_knots_x2 <- splinefun(x=mu$t, y=mu$X2, method="natural",  ties = mean)(knots)
    coefs <- as.matrix(data.frame(q_m_long.X1 = mu_knots_x1, q_m_long.X2 = mu_knots_x2))
    srv_data_curves <- procrustes_fits
    ### PROCRUSTES MEAN FINISH!
    
    # Old mean calculation!
    #coefs <- apply(model_data[, -1], 2, function(q_m_x_long) {
    #  q_m_x_long[!is.finite(q_m_x_long)] <- NA
    #  coef(lm(q_m_x_long ~ -1 + make_design(model_data[,1], knots = knots, closed = FALSE, type = type)))
    #})
    stop_crit <- sum((coefs - coefs_old)^2)/sum(coefs^2)
    if (stop_crit < eps | max_iter == 0) {
      
      # NEW: Diagnostic plots.
      plot.new()
      plot.window( xlim=c(-0.002,0.002), ylim=c(-0.002,0.002), asp = 1)
      mu_test <- get_points_from_srv(mu)
      data_curves_procrustes <- lapply(procrustes_fits, get_points_from_srv)
      lapply(data_curves_procrustes, lines, col = "gray")
      lines(mu_test, type = "l", col = "blue", lwd = 2)
      # FINISH
      
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



##################################################################
# Switch out original compute_elastic_mean with modified version #
##################################################################

# Add both to namespace.
environment(compute_elastic_mean) <- asNamespace('elasdics')
environment(fit_mean_proc2d) <- asNamespace('elasdics')
