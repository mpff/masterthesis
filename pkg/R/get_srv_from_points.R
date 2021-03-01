#' Helper functions for curve data measured at discrete points
#' @name get_srv_from_points
#' @aliases get_arc_length_param
#' @aliases get_points_from_srv
#' @description Compute the square-root-velocity transformation or the
#' parametrisation with respect to arc length for a curve
#' observed at discrete points.
#' @param data_curve A \code{data.frame} with observed points on a curve.
#' Each row is one point, each variable one coordinate direction. If there is a variable \code{t},
#' it is treated as the time parametrisation, not as an additional coordinate.
#' @param srv_data A \code{data.frame} with
#' first column \code{t} corresponding to the parametrisation and square-root-velocity
#' vectors in the remaining columns.
#' @examples
#' data_curve1 <- data.frame(x1 = 1:6*sin(1:6), x2 = cos(1:6))
#' get_arc_length_param(data_curve1) #same parametrisation as in
#' get_srv_from_points(data_curve1)
#'
#' data_curve2 <- data.frame(t = seq(0,1, length = 6), data_curve1)
#' plot(data_curve2[,2:3], type = "l", xlim = c(-6, 2), ylim = c(-2, 1))
#' srv_data <- get_srv_from_points(data_curve2)
#' #back transformed curve starts at (0,0)
#' lines(get_points_from_srv(srv_data), col = "red")


NULL

#' @describeIn get_srv_from_points Compute square-root-velocity transformation
#' for curve data measured at discrete points. The inverse transformation can
#' be computed with \code{get_points_from_s}
#' @export
#' @return \code{get_srv_from_points} returns a \code{data.frame} with
#' first column \code{t} corresponding to the parametrisation and square-root-velocity
#' vectors in the remaining columns. If no parametrisation is given, the curve will
#' be parametrised with respect to arc length. This parametrisation will be
#' computed by a call to \code{get_arc_length_param} as well.


get_srv_from_points <- function(data_curve){
  # parametrisation with respect to arg length if not given
  if(!("t" %in% colnames(data_curve))){
    data_curve <- data.frame("t" = get_arc_length_param(data_curve), data_curve)
  }
  #input checking of parametrisation
  if(!(data_curve$t[1] == 0 & data_curve$t[nrow(data_curve)] == 1))
    stop("Parametrisation must be starting at 0 and ending at 1")
  if(!all(diff(data_curve$t) > 0))
    stop("Parametrisation needs to be strictly increasing")

  data_points <- as.matrix(data_curve[, names(data_curve) != "t"])
  diff_points <- diff(data_points)
  factor <- 1/(sqrt(diff(data_curve$t))*rowSums(diff_points^2)^0.25)
  srv_vectors <- factor*diff_points
  srv_data <- data.frame("t" = data_curve$t[-nrow(data_curve)], srv_vectors)
  srv_data
}

#' @export
#' @describeIn get_srv_from_points The inverse transformation to
#' \code{get_srv_from_points}. Transforms square-root-velocity data to
#' points representing a curve (with no parametrisation).

get_points_from_srv <- function(srv_data){
  # extract srv vectors
  srv_vectors <- as.matrix(srv_data[, names(srv_data) != "t"])
  # extract parametrisation
  t <- c(srv_data$t, 1)
  norm <- sqrt(apply(srv_vectors^2, 1, sum))
  v <- apply(srv_vectors, 2, function(x) x*norm)
  path <- apply(v, 2, function(x) x*diff(t))
  points <- rbind(c(0,0), apply(path, 2, cumsum))
  data.frame(points)
}

#' @export
#' @describeIn get_srv_from_points Compute arc length parametrisation.

get_arc_length_param <- function(data_curve){
  # remove given parametrisation if it exists
  try(data_curve <- data_curve[,-t], silent = TRUE)

  t_arc_length <- c(0, cumsum(sqrt(rowSums(apply(t(data_curve), 1, diff) ^ 2))))
  t_arc_length <- t_arc_length/max(t_arc_length)
  t_arc_length
}


