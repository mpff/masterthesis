#' Compute a elastic mean for a collection of curves
#' @name compute_elastic_mean
#' @description Computes a Fréchet mean for the curves stored in \code{data_curves}) with respect
#' to the elastic distance. Constructor function for class \code{elastic_mean}.
#' @param data_curves list of \code{data.frame}s with observed points in each row. Each
#' variable is one coordinate direction. If there is a variable \code{t},
#' it is treated as the time parametrisation, not as an additional coordinate.
#' @param knots set of knots for the mean spline curve
#' @param type if "smooth" linear srv-splines are used which results in a differentiable mean curve
#' if "polygon" the mean will be piecewise linear.
#' @param closed \code{TRUE} if the curves should be treated as closed.
#' @param eps the algorithm stops if L2 norm of coefficients changes less
#' @param pen_factor penalty factor forcing the mean to be closed
#' @param max_iter maximal number of iterations
#' @return an object of class \code{elastic_mean}, which is a \code{list}
#' with entries
#'   \item{type}{"smooth" if mean was modeled using linear srv-splines or
#'   "polygon" if constant srv-splines are used}
#'   \item{coefs}{spline coeffiecients}
#'   \item{knots}{spline knots}
#'   \item{data_curves}{list of \code{data.frame}s with observed points in each row.
#'   First variable \code{t} gives the initial parametrisation, second variable \code{t_optim}
#'   the optimal parametrisation when the curve is aligned to the mean.}
#'   \item{closed}{\code{TRUE} if the mean is supposed to be a closed curve.}
#' @export
#' @exportClass elastic_mean
#' @importFrom splines splineDesign
#' @examples
#' curve <- function(t){
#'   rbind(t*cos(13*t), t*sin(13*t))
#' }
#' set.seed(18)
#' data_curves <- lapply(1:4, function(i){
#'   m <- sample(10:15, 1)
#'   delta <- abs(rnorm(m, mean = 1, sd = 0.05))
#'   t <- cumsum(delta)/sum(delta)
#'   data.frame(t(curve(t)) + 0.07*t*matrix(cumsum(rnorm(2*length(delta))),
#'              ncol = 2))
#' })
#'
#' #compute elastic means
#' knots <- seq(0,1, length = 11)
#' smooth_elastic_mean <- compute_elastic_mean(data_curves, knots = knots)
#' plot(smooth_elastic_mean)
#'
#' knots <- seq(0,1, length = 15)
#' polygon_elastic_mean <- compute_elastic_mean(data_curves, knots = knots, type = "poly")
#' lines(get_evals(polygon_elastic_mean), col = "blue", lwd = 2)
#'
#' #compute closed smooth mean, takes a little longer
#' knots <- seq(0,1, length = 11)
#' closed_elastic_mean <- compute_elastic_mean(data_curves, knots = knots, closed = TRUE)
#' plot(closed_elastic_mean)


compute_elastic_mean <- function(data_curves, knots = seq(0, 1, len = 5),
                                 type = c("smooth", "polygon"), closed = FALSE,
                                 eps = 0.01, pen_factor = 100, max_iter = 50) {
  #input checking
  stopifnot(all(sapply(data_curves, is.data.frame)))
  # remove duplicated points
  data_curves <- lapply(data_curves, remove_duplicate)
  if(sum(sapply(data_curves, function(curve){attributes(curve)$points_rm}) > 0)){
    warning("Duplicated points in data curves have been removed!")
  }
  data_curves <- lapply(data_curves, function(curve){
    attr(curve, "points_rm") <- NULL
    curve
  })
  # input checking given parametrisation t
  lapply(data_curves, function(data_curve){
    if("t" %in% names(data_curve)) check_param(data_curve, closed)
  })
  # input checking for closed curves
  if(closed){
    data_curves <- lapply(data_curves, function(data_curve){
      check_closed(data_curve)
    })
  }
  # parametrisation with respect to arc length if not given,
  # after this, parametrisation is always in the first column
  data_curves <- lapply(data_curves, function(data_curve){
    if(!("t" %in% colnames(data_curve))){
      data.frame("t" = get_arc_length_param(data_curve), data_curve)
    } else {
      param <- data_curve$t
      data_curve$t <- NULL
      data.frame("t" = param, data_curve)
    }
  })

  srv_data_curves <- lapply(data_curves, get_srv_from_points)

  type <- match.arg(type, c("smooth", "polygon"))
  if(!closed){
    elastic_mean <- fit_mean(srv_data_curves = srv_data_curves, knots = knots,
                           max_iter = max_iter, type = type, eps = eps)
    data_curves <- lapply(1:length(data_curves), function(j){
      data_curves[[j]]$t_optim <- elastic_mean$t_optims[[j]]
      attributes(data_curves[[j]]$t_optim) <- NULL
      data_curve <- data_curves[[j]][, c(1, ncol(data_curves[[j]]), 2:(ncol(data_curves[[j]]) - 1))]
      attr(data_curve, "dist_to_mean") <- attr(elastic_mean$t_optims[[j]], "dist_to_mean")
      data_curve
    })

  } else if(closed){
    elastic_mean <- fit_mean_closed(srv_data_curves = srv_data_curves, knots = knots,
                             max_iter = max_iter, type = type, eps = eps, pen_factor = pen_factor)
    data_curves <- lapply(1:length(data_curves), function(j){
      t_optim <- elastic_mean$t_optims[[j]][-length(elastic_mean$t_optims[[j]])]
      data_curve <- data_curves[[j]][-nrow(data_curves[[j]]), ]
      part2_idx <- 1:(length(t_optim) - elastic_mean$shift_idxs[j] + 1)
      data_curve$t_optim <- c(t_optim[-part2_idx], t_optim[part2_idx])
      data_curve <- data_curve[, c(1, ncol(data_curve), 2:(ncol(data_curve) - 1))]
      data_curve <- rbind(data_curve, data_curve[1,])
      data_curve$t[nrow(data_curve)] <- 1

      attr(data_curve, "dist_to_mean") <- attr(elastic_mean$t_optims[[j]], "dist_to_mean")
      data_curve
    })
  }

  elastic_mean$data_curves <- data_curves
  elastic_mean$closed <- closed
  elastic_mean$shift_idxs <- NULL
  elastic_mean$t_optims <- NULL
  class(elastic_mean) <- "elastic_mean"
  elastic_mean
}
