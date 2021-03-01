#' Plot method for aligned curves
#' @description Plots objects of class \code{aligned_curves}.
#' Points of same colour correspond after the second curve is optimally aligned to the first curve.
#' @param x object of class \code{aligned_curves}, usually a result of a call to \code{\link{align_curves}}
#' @param points_col which colour palette is used for points on the curves,
#' default is rainbow, see \code{\link{rainbow}} for further options.
#' @param ... further plotting parameters.
#' @importFrom grDevices rainbow
#' @importFrom graphics par  plot  points
#' @export
#'
#' @seealso For examples see documentation of \code{\link{align_curves}}.

plot.aligned_curves <- function(x, points_col = rainbow, ...){
  if(x$closed){
    x <- get_start_end(x)
  }

  warping <- get_warping(x)
  points1 <- get_evals(x$data_curve1, warping$t)
  x$data_curve2_aligned$t_optim <- NULL
  points2 <- get_evals(x$data_curve2_aligned, warping$gamma_t)

  par(mfrow = c(1,2))
  plot(points1, col = "gray", type = "l", ...)
  points(points1, pch = 16, col = points_col(nrow(points1)), ...)
  plot(points2, col = "gray", type = "l", ...)
  points(points2, pch = 16, col = points_col(nrow(points2)), ...)
  par(mfrow = c(1,1))
}

get_warping <- function(aligned_curves){
  srv_curve1 <- get_srv_from_points(aligned_curves$data_curve1)
  srv_curve2 <- get_srv_from_points(aligned_curves$data_curve2_aligned[,-2])
  t_optim <- aligned_curves$data_curve2_aligned$t_optim
  srv_curve2$t <- t_optim[-length(t_optim)]
  srv_curve2 <- srv_curve2[order(srv_curve2$t),]

  t_all <- unique(sort(c(aligned_curves$data_curve2_aligned$t_optim, aligned_curves$data_curve1$t)))

  idx1 <- sapply(t_all, findInterval, vec = srv_curve1$t)
  srv_curve1_all <- srv_curve1[idx1,]
  srv_curve1_all$t <- t_all

  idx2 <- sapply(t_all, findInterval, vec = srv_curve2$t)
  idx2 <- ifelse(idx2 == 0, max(idx2), idx2)
  srv_curve2_all <- srv_curve2[idx2,]
  srv_curve2_all$t <- t_all

  inner_prods <- sapply(rowSums(srv_curve1_all[,-1]*srv_curve2_all[,-1]), max, 0)^2
  inner_prod_parts <- cbind("idx" = idx2[-length(idx2)], diff(t_all)*inner_prods[-length(inner_prods)])
  idx_missing <- setdiff(1:(length(t_optim)-1), inner_prod_parts[,1])
  try(inner_prod_parts <- rbind(inner_prod_parts, cbind("idx" = idx_missing, 0)), silent = TRUE)
  inner_prod_parts <- inner_prod_parts[order(inner_prod_parts[,1]),]
  t_all <- sort(c(t_all, t_optim[idx_missing]))

  delta_gamma <- unlist(sapply(unique(inner_prod_parts[,1]), function(idx){
    x <- inner_prod_parts[inner_prod_parts[,1] == idx, 2]
    if(sum(x) == 0){
      diff(aligned_curves$data_curve2_aligned$t)[idx]
    } else {
      diff(aligned_curves$data_curve2_aligned$t)[idx]*x/sum(x)
    }
  }))

  data.frame("t" = t_all, "gamma_t" = c(0, cumsum(delta_gamma)))
}

get_start_end <- function(aligned_curves){
  srv_curve2 <- get_srv_from_points(aligned_curves$data_curve2_aligned[,-2])
  t_optim <- aligned_curves$data_curve2_aligned$t_optim
  t_optim <- ifelse(t_optim > 1, t_optim - 1, t_optim)
  srv_curve2$t_optim <- t_optim[-length(t_optim)]
  srv_curve2 <- srv_curve2[order(srv_curve2$t_optim),]
  srv_curve2$t[-length(srv_curve2$t)] <- ifelse(srv_curve2$t[-length(srv_curve2$t)] == 0,
                                                1, srv_curve2$t[-length(srv_curve2$t)])
  if(srv_curve2$t_optim[1] == 0){
    start <- aligned_curves$data_curve2_aligned[round(aligned_curves$data_curve2_aligned$t_optim, 10) == 0,]
  } else {
    delta_s <- diff(srv_curve2$t[c(nrow(srv_curve2), 1)])
    q <- srv_curve2[nrow(srv_curve2), -1]
    q$t_optim <- NULL

    srv_curve1 <- get_srv_from_points(aligned_curves$data_curve1)
    t_all <- unique(sort(c(srv_curve2$t_optim[c(nrow(srv_curve2), 1)], aligned_curves$data_curve1$t)))
    idx1 <- sapply(t_all, findInterval, vec = srv_curve1$t)
    srv_curve1_all <- srv_curve1[idx1,]
    srv_curve1_all$t <- t_all
    p_right <- srv_curve1_all[srv_curve1_all$t <= srv_curve2$t_optim[1], ]
    p_left <- srv_curve1_all[srv_curve1_all$t >= srv_curve2$t_optim[length(srv_curve2$t_optim)], ]

    inner_prods_right <- sapply(as.matrix(p_right[-nrow(p_right),-1])%*%t(q), max, 0)^2*diff(p_right[,1])
    inner_prods_left <- sapply(as.matrix(p_left[-nrow(p_left),-1])%*%t(q), max, 0)*diff(p_left[,1])

    delta_gamma_left <- delta_s*inner_prods_left/(inner_prods_right + inner_prods_left)
    start_s <- srv_curve2$t[1] - delta_gamma_left
    start <- get_evals(aligned_curves$data_curve2_aligned[,-2], start_s)
    start <- data.frame("t" = start_s, "t_optim" = 0, start)
  }

  #remove first point
  data_curve2 <- aligned_curves$data_curve2_aligned[-1,]
  #change order and add start
  data_curve2 <- unique(rbind(start, data_curve2[order(data_curve2$t_optim),]))
  data_curve2$t <- data_curve2$t - data_curve2$t[1]
  data_curve2$t <- ifelse(data_curve2$t < 0, data_curve2$t + 1, data_curve2$t)
  #add last point
  data_curve2 <- rbind(data_curve2, c("t"= 1, "t_optim" = 1, data_curve2[1, -(1:2)]))

  aligned_curves$data_curve2_aligned <- data_curve2
  aligned_curves
}
