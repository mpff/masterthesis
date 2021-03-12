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
  require(sparseFLMM)
  require(orthogonalsplinebasis)

  # OLD: initialize t_optimns and coefs.
  t_optims <- lapply(srv_data_curves, function(srv_data_curve) {
    c(srv_data_curve$t, 1)
  })
  coefs <- 0

  # NEW: initialize G_optimns as 0 and b_optimns as 1 (Rotation and Scale).
  # ToDo: calculate G,b instead of directly calculating the procrustes fits.
  G_optims <- as.list(rep(0, length(srv_data_curves)))
  b_optims <- as.list(rep(1, length(srv_data_curves)))

  # NEW: Build id column, super hacky way. THIS NEEDS REWORK!
  ids <- do.call(c, lapply(t_optims, function(x) length(x)-1))  #note: t_optims has 1 more
  ids <- rbind(1:length(ids), ids)
  ids <- apply(ids, 2, function(x) rep(x[1], times = x[2]))
  ids <- do.call(c, as.list(ids))  # coerce as.list for curves of equal length(t)

  # Parameters for covariance smoothing
  # Note: The parameters here have a huge influence on the results!!!
  cov.m = ifelse(type == "smooth", 0, -1)
  cov.d = 1 # penalty
  knotl = 1 / ( length(knots) - 1 )  # mean length of a knot
  cov.knots = c(rep(-knotl,cov.m+1), knots, rep(1+knotl,cov.m+1))  # asusmes cov.m is 0 or -1.
  cov.k = length(cov.knots) - cov.m - 2  # Why use this?

  # NEW: Build arg grid for evaluation of scalar products..
  arg.grid <- seq(0,1, len = 101)

  # Number of curves
  n <- length(srv_data_curves)

  # Define mean basis.
  smooth1D <- smooth.construct(s(t, bs="ps", k = cov.k, m = c(cov.m, cov.d),
                                 fx = FALSE),
                               data = list(t=arg.grid), knots=list(t=cov.knots))

  # Get Gram and trafo matrix for mean basis.
  G <- build_gram_matrix(smooth1D)
  Ginv <- solve(G)
  L <- chol(Ginv)


  for (i in 1:max_iter) {

    # OLD: Evaluate the srv_data_curves at points t_optims. Get q(m_optim).
    model_data <- get_model_data(t_optims, srv_data_curves, knots, type)
    coefs_old <- coefs


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

    # Smooth covariance surface
    cov_fit_re <- bam(Re(qq) ~ s(t, s, bs="symm", k = cov.k, m = c(cov.m, cov.d),
                                 fx = FALSE, xt = list(skew = FALSE)),
                      data = cov_dat, method = "REML", knots=list(t = cov.knots, s = cov.knots))
    cov_fit_im <- bam(Im(qq) ~ -1 + s(t, s, bs="symm", k = cov.k, m = c(cov.m, cov.d),
                                      fx = FALSE, xt = list(skew = TRUE)),
                      data = cov_dat, method = "REML", knots=list(t = cov.knots, s = cov.knots))

    # Get coefficient matrix
    beta.mat.re <- get_coef_matrix(cov_fit_re)
    beta.mat.im <- get_coef_matrix(cov_fit_im)
    beta.mat <- matrix(
        complex(real = as.vector(beta.mat.re), imaginary = as.vector(beta.mat.im)),
        ncol = cov.k
    )

    # Calculate inverse coefficint matrix for calculation of procrustes fits.
    beta.mat.inv <- solve(beta.mat)

    # Calculate largest eigenvector
    pca <- eigen(t(L) %*% beta.mat %*% L)
    coefs.compl <- pca$vectors[,1]
    coefs <- as.matrix(data.frame(q_m_long.X1 = Re(coefs.compl), q_m_long.X2 = Im(coefs.compl)))



    ###################################
    # Calculate Procrustes SRV Curves #
    ###################################

    srv_procrustes_curves <- lapply(1:length(srv_data_curves), function(j) {
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
        mean_eval <- make_design(arg.grid.x, knots=knots, closed = FALSE, type = type) %*% coefs
        mean_eval <- complex(real = mean_eval[,1], imaginary = mean_eval[,2])
        # Calculate qm and qq over arg.grid.x
        qm <- Conj(q_eval) %*% mean_eval
        qq <- Conj(q_eval) %*% q_eval
        # Calculate G and b (Note: using "<<-" is not good practice...)
        G_optims[j] <<- Arg(qm)
        # Note: squared to adjust for SRV framework. b is on data_curve level!!!
        b_optims[j] <<- (Mod(qm) / Re(qq))^2
        # Calculate procrustes fit of original srv_data_curve
        # Note: warping is performed on the original srv_data_curves not on model_data!
        srv_complex = complex(real = srv_data_curves[[j]][,2], imaginary = srv_data_curves[[j]][,3])
        pfit <- as.vector(qm) * srv_complex / as.vector(qq)
        data.frame(t = srv_data_curves[[j]][,1], X1 = Re(pfit), X2 = Im(pfit))
    })

    srv_procrustes_curves <- lapply(1:length(srv_data_curves), function(j) {
        # Grab warped srv_data_curve from model_data.
        q <- model_data_complex[model_data_complex$id == j,]
        # Create design matrix and coefs in mean basis for q.
        q_B <- make_design(q$m_long, knots)
        q_coefs <- solve(t(Conj(q_B)) %*% q_B + t(beta.mat.inv)) %*% t(Conj(q_B)) %*% q$q_m_long
        # Calculate qm and qq
        qm <- t(Conj(q_coefs)) %*% G %*% Conj(coefs.compl)
        qq <- t(Conj(q_coefs)) %*% G %*% q_coefs
        # Calculate G and b (Note: using "<<-" is not good practice...)
        G_optims[j] <<- Arg(qm)
        # Note: squared to adjust for SRV framework. b is on data_curve level!!!
        b_optims[j] <<- (Mod(qm) / Re(qq))^2
        # Calculate procrustes fit of original srv_data_curve
        #srv_complex = complex(real = srv_data_curves[[j]][,2], imaginary = srv_data_curves[[j]][,3])
        #pfit <- c(qm/qq) * srv_complex
        pfit <- c(qm/qq) * q[,3]
        #data.frame(t = srv_data_curves[[j]][,1], X1 = Re(pfit), X2 = Im(pfit))
        data.frame(t = q[,2], X1 = Re(pfit), X2 = Im(pfit))
    })

    # OLD: Stopping criteria step.
    stop_crit <- sum((coefs - coefs_old)^2)/sum(coefs^2)
    if (stop_crit < eps | max_iter == 0) {

      # NEW: calculate procrustes data curves on basis of srv.
      procrustes_curves <- lapply(srv_procrustes_curves, get_points_from_srv)
      procrustes_curves <- lapply(procrustes_curves, center_curve)

      # Prepare output.
      rownames(coefs) <- NULL
      colnames(coefs) <- colnames(srv_data_curves[[1]][,-1])
      return(list(type = type, coefs = coefs, knots = knots,
                  t_optims = t_optims, G_optims = G_optims, b_optims = b_optims,
                  procrustes_curves = procrustes_curves))
    }


    # NEW: Calculate warping on the procrustes fits!
    if (type == "smooth") {
      pfun <- function(t) {
        t(make_design(t, knots = knots, closed = FALSE, type = type) %*% coefs)
      }
      t_optims <- lapply(1:length(srv_procrustes_curves), function(j) {
        t_optim <- find_optimal_t(srv_curve = pfun,
                                  s = c(srv_procrustes_curves[[j]]$t,1),
                                  q = t(srv_procrustes_curves[[j]][, -1]),
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
      t_optims <- lapply(1:length(srv_procrustes_curves), function(j) {
        t_optim <- find_optimal_t_discrete(r = knots,
                                           p = t(coefs),
                                           s = c(srv_procrustes_curves[[j]]$t, 1),
                                           q = t(srv_procrustes_curves[[j]][, -1]),
                                           initial_t = t_optims[[j]]
                                          )
        attr(t_optim, "dist_to_mean") <- attr(t_optim, "dist")
        attr(t_optim, "dist") <- NULL
        t_optim
      })
    }
  }
  warning("Stopping criteria eps has not been reached! Consider more iterations max_iter")

  # NEW: calculate procrustes data curves on basis of srv.
  procrustes_curves <- lapply(srv_procrustes_curves, get_points_from_srv)
  procrustes_curves <- lapply(procrustes_curves, center_curve)

  # Prepare output.
  rownames(coefs) <- NULL
  colnames(coefs) <- colnames(srv_data_curves[[1]][,-1])
  return(list(type = type, coefs = coefs, knots = knots,
              t_optims = t_optims, G_optims = G_optims, b_optims = b_optims,
              procrustes_curves = procrustes_curves))
}

build_gram_matrix <- function(smooth){
    # Takes mgcv smooth and calculates Gram matrix.
    # Only works for b-splines.
    order = smooth$m[1] + 2  # degree + 1
    knots = smooth$knots  # inner+outer knots
    if( order == 1 ){
        diag(smooth$bs.dim)
    } else {
        osb_smooth = orthogonalsplinebasis::SplineBasis(knots,order=order)
        orthogonalsplinebasis::GramMatrix(osb_smooth)
    }
}

get_coef_matrix <- function(model){
    F <- model$smooth[[1]]$bs.dim
    beta <- model$smooth[[1]]$Z %*% model$coefficients
    matrix(beta, nrow=F, ncol=F)
}

make_design <- function (t, knots, closed = FALSE, type = "smooth")
{
    deg <- ifelse(type == "smooth", 1, 0)
    knotl = 1 / ( length(knots) - 1 )  # mean length of a knot
    # Switch add outer knots outside of [0,1] etc. !
    design_mat <- splineDesign(
        knots = c(rep(-knotl,deg), knots, rep(1+knotl,deg)),  # assumes deg is 0 or 1.
        x = t,
        ord = deg + 1
    )
    if (closed == TRUE & type == "smooth") {
        design_mat[, 1] <- design_mat[, 1] + design_mat[, ncol(design_mat)]
        design_mat <- design_mat[, -ncol(design_mat)]
    }
    design_mat
}

##################################################################
# Switch out original compute_elastic_mean with modified version #
##################################################################

# Add both to namespace.
environment(compute_elastic_mean) <- asNamespace('elasdics')
environment(fit_mean_proc2d) <- asNamespace('elasdics')
environment(get_proc2d_fit) <- asNamespace('elasdics')
environment(get_coef_matrix) <- asNamespace('elasdics')
environment(build_gram_matrix) <- asNamespace('elasdics')
environment(make_design) <- asNamespace('elasdics')
