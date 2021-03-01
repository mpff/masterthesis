
# Subroutines for randomly scaling and rotating curves, and for centering curves

rand_rotate_curve2d <- function(x){
    # rotate dataframe of 2D vectors randomly
    names <- colnames(x)
    theta <- 2*pi*runif(1)
    mat <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2)
    x.rot <- as.matrix(x) %*% t(mat)
    x <- as.data.frame(x.rot)
    colnames(x) <- names
    return(x)
}

rand_scale_curve2d <- function(x){
    # scale dataframe of 2D vectors randomly
    beta <- 0.5 + runif(1)
    beta * x
}

center_curve <- function (data_curve) 
{
    coord_idx <- !(colnames(data_curve) %in% c("t", "id"))
    data_curve[, coord_idx] <- data_curve[, coord_idx] - matrix(colMeans(data_curve[, coord_idx]), 
                                                                nrow = nrow(data_curve), 
                                                                ncol = ncol(data_curve[, coord_idx]), 
                                                                byrow = TRUE
                                                               )
    data_curve
}


# Datasets

curves.spiral <- function(n_curves = 4, rotate = FALSE, scale = FALSE, center = TRUE){
    # Dataset: Simulated spirals with optional random rotation and scaling.
    
    # Define spiral shape.
    curve <- function(t){
      rbind(t*cos(13*t), t*sin(13*t))
    }
    
    # Transform to data curves format.
    data_curves <- lapply(1:n_curves, function(i){
      m <- sample(10:15, 1)
      delta <- abs(rnorm(m, mean = 1, sd = 0.05))
      t <- cumsum(delta)/sum(delta)
      data.frame(t(curve(t)) + 0.07*t*matrix(cumsum(rnorm(2*length(delta))), ncol = 2))
    })

    # Apply random rotation, random scaling, centering.
    if(rotate){
        data_curves <- lapply(data_curves, rand_rotate_curve2d)
    }
    if(scale){
        data_curves <- lapply(data_curves, rand_scale_curve2d)
    }
    if(center){
        data_curves <- lapply(data_curves, center_curve)
    }
    
    
    return(data_curves)
}


curves.digit3 <- function(rotate = FALSE, scale = FALSE, center = TRUE){
    # Dataset : Handwritten digits 3.
    
    data_curves <- shapes::digit3.dat
    data_curves <- apply(data_curves, MARGIN = 3, FUN = function(i){
      data.frame(X1 = i[,1], X2 = i[,2])
    })
    
    # Apply random rotation, random scaling, centering.
    if(rotate){
        data_curves <- lapply(data_curves, rand_rotate_curve2d)
    }
    if(scale){
        data_curves <- lapply(data_curves, rand_scale_curve2d)
    }
    if(center){
        data_curves <- lapply(data_curves, center_curve)
    }
    
    return(data_curves)
}

