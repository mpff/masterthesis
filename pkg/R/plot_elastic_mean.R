#' Plot method for planar elastic mean curves
#' @description Plots objects of class \code{elastic_mean}.
#' @param x object of class \code{elastic_mean},
#' usually a result of a call to \code{\link{compute_elastic_mean}}
#' @param asp numeric, giving the aspect ratio of the two coordinates,
#' see \code{\link{plot.window}} for details.
#' @param col color of the mean curve.
#' @param ... further plotting parameters.
#' @importFrom graphics plot lines
#' @export
#'
#' @seealso For examples see documentation of \code{\link{compute_elastic_mean}}.

plot.elastic_mean <- function(x, asp = 1, col = "red", ...){
  if(ncol(x$coefs) != 2){
    stop("Plotting option only for planar curves")
  }
  data_curves <- lapply(x$data_curves, center_curve)
  data_curves <- lapply(data_curves, function(data) data[,colnames(x$coefs)])
  data_all <- do.call("rbind", data_curves)

  #empty plot
  plot(NULL, xlim = range(data_all[,1]), ylim = range(data_all[,2]), xlab = colnames(x$coefs)[1],
       ylab = colnames(x$coefs)[2], asp = asp)
  #plot data

  lapply(data_curves, lines, col = "gray")

  #plot mean
  lines(get_evals(x), col = col, lwd = 2)
}
