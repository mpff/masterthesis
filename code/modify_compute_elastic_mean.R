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
      elastic_mean <- fit_mean(srv_data_curves = srv_data_curves, 
                               knots = knots, max_iter = max_iter, type = type, 
                               eps = eps)
    }
    else if (proc2d){
      # This is new.
      elastic_mean <- fit_mean_proc2d(srv_data_curves = srv_data_curves, 
                                      knots = knots, max_iter = max_iter, type = type, 
                                    eps = eps)
      elastic_mean$procrustes_fits <- lapply(1:length(data_curves), function(j){
          get_proc2d_fit(data_curves[[1]][,-1], 
                         G = elastic_mean$G_optims[[j]], 
                         b = elastic_mean$b_optims[[j]])
      })
    }
    # This just prepares output
    data_curves <- lapply(1:length(data_curves), function(j) {
      data_curves[[j]]$t_optim <- elastic_mean$t_optims[[j]]
      attributes(data_curves[[j]]$t_optim) <- NULL
      data_curve <- data_curves[[j]][, c(1, ncol(data_curves[[j]]), 
                                         2:(ncol(data_curves[[j]]) - 1))]
      attr(data_curve, "dist_to_mean") <- attr(elastic_mean$t_optims[[j]], "dist_to_mean")
      data_curve
    })
  }
  else if (closed) {
    if (proc2d) {
        warning("Option 'proc2d' does not (yet) work with closed curves. Calculating non-procrustes mean.")
    }
    # closed curves.
    elastic_mean <- fit_mean_closed(srv_data_curves = srv_data_curves, 
                                    knots = knots, max_iter = max_iter, type = type, 
                                    eps = eps, pen_factor = pen_factor)
    
    # This just prepares output I think. What's the additional stuff??
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
    
  # OLD: initialize t_optimns and coefs.
  t_optims <- lapply(srv_data_curves, function(srv_data_curve) {
    c(srv_data_curve$t, 1)
  })
  coefs_knots <- 0
    
  # NEW: initialize G_optimns as 0 and b_optimns as 1 (Rotation and Scale).
  # ToDo: calculate G,b instead of directly calculating the procrustes fits.
  G_optims <- as.list(rep(0, length(srv_data_curves)))
  b_optims <- as.list(rep(1, length(srv_data_curves)))
    
  # NEW: Build id column, super hacky way. THIS NEEDS REWORK!
  ids <- do.call(c, lapply(t_optims, function(x) length(x)-1))  #note: t_optims has 1 more
  ids <- rbind(1:length(ids), ids)
  ids <- apply(ids, 2, function(x) rep(x[1], times = x[2]))
  ids <- do.call(c, as.list(ids))  # coerce as.list for curves of equal length(t)
  
  # NEW: Build arg grid for evaluation of smoothed cov.
  arg.grid <- seq(0,1, len = 200)
  
  for (i in 1:max_iter) {
    
    # OLD: Evaluate the srv_data_curves at points t_optims. Get q(m_optim).
    model_data <- get_model_data(t_optims, srv_data_curves, knots, type)
    coefs_old <- coefs_knots
    
    
    #######################
    # Fit Procrustes Mean #
    #######################
      
    # x,y to complex, add id column
    model_data_complex <- complex(re = model_data[,2], im = model_data[,3]) %>% matrix(nrow = dim(model_data)[1])
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
    
    # fit covariance surface. 
    # Note: The parameters here have a huge influence on the results!!!
    # ToDo: Find good standard parameters!
    cov.m = ifelse(type == "smooth", 2, 1)
    cov.knots = c(rep(0,cov.m+1), knots, rep(1,cov.m+1))
    cov.k = length(cov.knots) - cov.m - 2  # Why use this?
    cov_fit_re <- bam(Re(qq) ~ s(t, s, bs = "sps", k = cov.k, m = c(cov.m,0),
                                 fx = FALSE, xt = list(skew = FALSE)),
                      data = cov_dat, knots=list(t = cov.knots, s = cov.knots), method = "REML")
    cov_fit_im <- bam(Im(qq) ~ -1 + s(t, s, bs = "sps", k = cov.k, m = c(cov.m,0),
                                      fx = FALSE, xt = list(skew = TRUE)),
                      data = cov_dat, knots=list(t = cov.knots, s = cov.knots), method = "REML")
    
    # predict smoothed covariance surface on arg.grid
    cov_dat <- expand.grid(t = arg.grid, s = arg.grid)
    yy <- matrix( complex(
      real = predict(cov_fit_re, newdata = cov_dat, outer.ok=TRUE),
      imaginary = predict(cov_fit_im, newdata = cov_dat, outer.ok=TRUE) ),
      ncol = length(arg.grid) )
    
    # calculate procrustes mean (leading eigenvector)
    ei <- eigen(yy)
    coefs_smooth <- as.matrix(data.frame(q_m_long.X1 = Re(ei$vectors[,1]), q_m_long.X2 = Im(ei$vectors[,1])))
    coefs_knots <- make_design(knots, knots = arg.grid, closed = FALSE, type = type) %*% coefs_smooth
    coefs_knots <- as.matrix(data.frame(q_m_long.X1 = coefs_knots[,1], q_m_long.X2 = coefs_knots[,2]))
    
      
    ###################################
    # Calculate Procrustes SRV Curves #
    ###################################
      
    procrustes_fits <- lapply(1:length(srv_data_curves), function(j) {
        # Grab warped srv_data_curve from model_data.
        x <- model_data[model_data_complex$id == j,]
        # !!!! Calculate overlap of arg.grid and t_optims (is this ok???)
        idx <- findInterval(arg.grid, x$m_long)
        idx.bool <- which(idx > 0 & idx < length(x$m_long))
        arg.grid.x <- arg.grid[idx.bool]
        # Treat srv_data_curve as function and evaluate on overlap.
        # Note: how to smooth srv_data_curve? -> depends on "type"
        q_coefs <- as.matrix(x[,-1])
        q_eval <- make_design(arg.grid.x, knots = x$m_long, closed = FALSE, type = type) %*% q_coefs
        q_eval <- complex(real = q_eval[,1], imaginary = q_eval[,2])
        # Calculate qm and qq over arg.grid.x
        qm <- Conj(q_eval) %*% ei$vectors[idx.bool,1]
        qq <- Conj(q_eval) %*% q_eval
        # Calculate G and b (Note: using "<<-" is not good practice...)
        G_optims[j] <<- Arg(qm)
        b_optims[j] <<- (Mod(qm) / Re(qq))^2  # Note: squared to adjust for SRV framework. b is on data_curve level!!!
        # Calculate procrustes fit of original srv_data_curve
        # Note: warping is performed on the original srv_data_curves not on model_data!
        srv_complex = complex(real = srv_data_curves[[j]][,2], imaginary = srv_data_curves[[j]][,3])
        pfit <- as.vector(qm) * srv_complex / as.vector(qq)
        data.frame(t = srv_data_curves[[j]][,1], X1 = Re(pfit), X2 = Im(pfit))
    })
         
      
    # OLD: Stopping criteria step.
    stop_crit <- sum((coefs_knots - coefs_old)^2)/sum(coefs_knots^2)
    if (stop_crit < eps | max_iter == 0) {
      print(paste("Number of basis functions: ",length(cov_fit_re$coefficients)))
      print("Knots (from):")
      print(cov_fit_re$smooth[[1]]$knots)
      print("-------------------------------")
      print("Printing summary and plotting fit")
      print(summary(cov_fit_re))
      plot(cov_fit_re)
      print("------------END----------------")

      # NEW: calculate procrustes data curves on basis of srv.
      data_curves_procrustes <- lapply(procrustes_fits, get_points_from_srv)
      data_curves_procrustes <- lapply(data_curves_procrustes, center_curve)
      # Prepare output.
      rownames(coefs_smooth) <- NULL
      colnames(coefs_smooth) <- colnames(srv_data_curves[[1]][,-1])
      rownames(coefs_knots) <- NULL
      colnames(coefs_knots) <- colnames(srv_data_curves[[1]][,-1])
      return(list(type = type, coefs = coefs_knots, knots = knots,  
                  t_optims = t_optims, G_optims = G_optims, b_optims = b_optims,
                  coefs_smooth = coefs_smooth, knots_smooth = arg.grid,
                  procrustes_fits_on_srv_basis = data_curves_procrustes,
                  model_re = cov_fit_re, model_im = cov_fit_im))
    }
    # NEW: Calculate warping on the procrustes fits!
    if (type == "smooth") {
      pfun <- function(t) {
        t(make_design(t, knots = knots, closed = FALSE, type = type) %*% coefs_knots)  # knots = arg.grid here!
      }
      t_optims <- lapply(1:length(srv_data_curves), function(j) {
        t_optim <- find_optimal_t(srv_curve = pfun,
                                  s = c(procrustes_fits[[j]]$t,1), 
                                  q = t(procrustes_fits[[j]][, -1]), 
                                  initial_t = t_optims[[j]], 
                                  eps = eps * 100/i
                                 )
        attr(t_optim, "dist_to_mean") <- attr(t_optim, "dist")
        attr(t_optim, "dist") <- NULL
        t_optim
      })
    }
    # ToDo: Check here! knots and coefs might have to fit!
    else {
      t_optims <- lapply(1:length(srv_data_curves), function(j) {
        t_optim <- find_optimal_t_discrete(r = knots, 
                                           p = t(coefs_knots), 
                                           s = c(procrustes_fits[[j]]$t, 1),
                                           q = t(procrustes_fits[[j]][, -1]),
                                           initial_t = t_optims[[j]]
                                          )
        attr(t_optim, "dist_to_mean") <- attr(t_optim, "dist")
        attr(t_optim, "dist") <- NULL
        t_optim
      })
    }
  }
  warning("Stopping criteria eps has not been reached! Consider more iterations max_iter")
               
  # NEW: calculate procrustes data curves on the basis of srv curves.
  data_curves_procrustes <- lapply(procrustes_fits, get_points_from_srv)
  data_curves_procrustes <- lapply(data_curves_procrustes, center_curve)
  # NEW: evaluate procrustes mean at knots (necessarry for elasdics::get_evals() !)
  coefs_knots <- make_design(knots, knots = arg.grid, closed = FALSE, type = type) %*% coefs
  coefs_knots <- as.matrix(data.frame(q_m_long.X1 = coefs_knots[,1], q_m_long.X2 = coefs_knots[,2]))
  rownames(coefs_smooth) <- NULL
  colnames(coefs_smooth) <- colnames(srv_data_curves[[1]][,-1])
  rownames(coefs_knots) <- NULL
  colnames(coefs_knots) <- colnames(srv_data_curves[[1]][,-1])
  # NEW: also return procrustes fits.
  data_curves_procrustes <- lapply(procrustes_fits, get_points_from_srv)
  return(list(type = type, coefs = coefs_knots, knots = knots,  # Note: also returns procrustes_fits. TODO: G,b!
              t_optims = t_optims, G_optims = G_optims, b_optims = b_optims,
              coefs_smooth = coefs_smooth, knots_smooth = arg.grid, 
              procrustes_fits_on_srv_basis = data_curves_procrustes))
}

               
# Helper function to calculate procrustes fits from G_optim and b_optim.          
get_proc2d_fit <- function(data_curve, G, b)
{
    # Procrustes fit: beta * e^(i*G) * curve
    names <- colnames(data_curve)
    mat <- matrix(c(cos(G), sin(G), - sin(G), cos(G)), nrow = 2, ncol = 2)
    data_curve.rot <- as.matrix(data_curve) %*% t(mat)
    data_curve.rot.scale <- b * data_curve.rot
    data_curve <- as.data.frame(data_curve.rot.scale)
    colnames(data_curve) <- names
    return(data_curve)
}

               

##################################################################
# Switch out original compute_elastic_mean with modified version #
##################################################################

# Add both to namespace.
environment(compute_elastic_mean) <- asNamespace('elasdics')
environment(fit_mean_proc2d) <- asNamespace('elasdics')
environment(get_proc2d_fit) <- asNamespace('elasdics')
