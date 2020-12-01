compute_elastic_proc2d_mean <- function(data_curves, 
                                        knots = seq(0, 1, len = 5),
                                        type = c("smooth", "polygon"), 
                                        closed = FALSE,
                                        eps = 0.01,
                                        pen_factor = 100,
                                        max_iter = 50) {

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
    # open curves. Calculate mean!
    elastic_mean <- fit_mean(srv_data_curves = srv_data_curves, 
                             knots = knots, max_iter = max_iter, type = type, 
                             eps = eps)
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
    # closed curves. Calculate CLOSED mean!
    elastic_mean <- fit_mean_closed(srv_data_curves = srv_data_curves, 
                                    knots = knots, max_iter = max_iter, type = type, 
                                    eps = eps, pen_factor = pen_factor)
    
    # This juut prepares output I think. What's the additional stuff????
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
  
  # Prepare Final Output!
  elastic_mean$data_curves <- data_curves
  elastic_mean$closed <- closed
  elastic_mean$shift_idxs <- NULL
  elastic_mean$t_optims <- NULL
  class(elastic_mean) <- "elastic_mean"
  elastic_mean
}