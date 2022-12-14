getPLOT <- function(fn, resolution = 500L) {
  design <- moPLOT::generateDesign(fn, points.per.dimension = resolution, evaluate.points = TRUE)
  
  gradients <- computeGradientFieldGrid(design, impute.boundary = TRUE)
  
  divergence <- computeDivergenceGrid(gradients$multi.objective, design$dims, design$step.sizes)
  
  less <- localEfficientSetSkeleton(
    design, gradients, divergence, integration = "fast", verbose = FALSE
  )
  
  plot_dec <- ggplotPLOT(design$dec.space, design$obj.space, less$sinks, less$height,
                         colorscale.efficient = viridis::turbo(500, begin = 0.05, end = 0.95)) +
    theme_bw() +
    coord_fixed() +
    theme(legend.position = "none")
  
  plot_obj <- ggplotPLOTObjSpace(design$obj.space, less$sinks, less$height) +
    scale_fill_gradientn(colors = viridis::turbo(500, begin = 0.05, end = 0.95), trans = "log") +
    lims(x = c(0,1), y = c(0,1)) +
    theme_bw() +
    coord_fixed() +
    theme(legend.position = "none")
  
  list(plot_dec = plot_dec, plot_obj = plot_obj)
}

get_peak_metadata <- function(npeaks, dimension, topology,
                              randomSeed, rotated, peakShape) {
  
  xopt <- getAllPeaks(npeaks, dimension, topology,
                        randomSeed, rotated, peakShape)
  colnames(xopt) <- paste0("x", 1:dimension)
  
  cov <- getCovarianceMatrices(npeaks, dimension, topology,
                               randomSeed, rotated, peakShape) %>%
    lapply(function(x) matrix(x[[1]], nrow = dimension, ncol = dimension))
  
  height <- getAllHeights(npeaks, dimension, topology,
                          randomSeed, rotated, peakShape)
  
  shape <- getAllShapes(npeaks, dimension, topology,
                        randomSeed, rotated, peakShape)
  
  radius <- getAllRadii(npeaks, dimension, topology,
                        randomSeed, rotated, peakShape)
  
  peak_fns <- lapply(1:nrow(xopt), function(i) {
    create_peak_function(cov[[i]], xopt[i,], height[i], shape[i], radius[i])
  })
  
  fn <- function(x) {
    min(sapply(peak_fns, function(f) f(x)))
  }
  
  list(
    xopt = xopt,
    cov = cov,
    height = height,
    shape = shape,
    radius = radius,
    peak_fns = peak_fns,
    fn = fn
  )
}

create_peak_function <- function(cov, xopt, height, shape, radius) {
  function(x) {
    if (is.matrix(x)) {
      dx <- t(apply(x, 1, function(row) (row - xopt))) %*% chol(cov)
      md <- apply(dx, 1, function(row) sqrt(sum(row**2)))
    } else {
      md <- sqrt(t(x - xopt) %*% cov %*% (x - xopt))
      # dx <- (x - xopt) %*% chol(cov)
      # md <- sqrt(sum(dx**2))
    }
    g <- height / (1 + md**shape / radius)
    return(1 - g)
  }
}

compute_set_2o <- function(Cov1, Cov2, xopt1, xopt2, resolution = 1000L) {
  t(sapply(seq(0, 1, length.out = resolution), function(t) {
    solve((1 - t) * Cov1 + t * Cov2) %*% ((1 - t) * Cov1 %*% xopt1 + t * Cov2 %*% xopt2)
  }))
}

compute_set_3o <- function(Cov1, Cov2, Cov3, xopt1, xopt2, xopt3, resolution = 100L) {
  points <- seq(0, 1, length.out = resolution)
  points <- expand_grid(s = points, t = points)
  
  t(apply(points, 1, function(row) {
    s <- row[1]
    t <- row[2]
    
    if (s + t > 1) {
      s <- 1 - s
      t <- 1 - t
    }
    
    solve((1 - s - t) * Cov1 + t * Cov2 + s * Cov3) %*%
      ((1 - s - t) * Cov1 %*% xopt1 + t * Cov2 %*% xopt2 + s * Cov3 %*% xopt3)
  }))
}

get_optimistic_y <- function(y_nondominated) {
  y_nondominated <- y_nondominated[order(y_nondominated[,1]),]
  y_optimistic <- matrix(nrow = nrow(y_nondominated) - 1, ncol = ncol(y_nondominated))
  
  y_optimistic[,1] <- y_nondominated[1:(nrow(y_nondominated) - 1),1]
  y_optimistic[,2] <- y_nondominated[2:nrow(y_nondominated),2]
  
  return(y_optimistic)
}

get_pessimistic_y <- function(y_nondominated) {
  y_nondominated <- y_nondominated[order(y_nondominated[,1]),]
  y_pessimistic <- matrix(nrow = nrow(y_nondominated) - 1, ncol = ncol(y_nondominated))
  
  y_pessimistic[,1] <- y_nondominated[2:nrow(y_nondominated),1]
  y_pessimistic[,2] <- y_nondominated[1:(nrow(y_nondominated) - 1),2]
  
  return(y_pessimistic)
}

hv_gap <- function(y_nondominated, ref_point = c(1, 1)) {
  hv <- ecr::computeHV(t(y_nondominated), ref.point = ref_point)
  
  y_optimistic <- get_optimistic_y(y_nondominated)
  
  hv_optimistic <- ecr::computeHV(t(y_optimistic), ref.point = ref_point)
  
  list(
    hv = hv,
    gap = hv_optimistic - hv
  )
}

approximate_hv <- function(peak_data_1, peak_data_2, initial_resolution = 4L,
                           gap_target = 1e-5) {
  gap <- Inf
  hv <- 0
  
  set_resolution <- initial_resolution
  
  peak_combinations <- expand_grid(p1 = 1:nrow(peak_data_1$xopt),
                                   p2 = 1:nrow(peak_data_2$xopt))
  
  while (gap / hv > gap_target) {
    sets_individual <- lapply(1:nrow(peak_combinations), function(i) {
      p1 <- peak_combinations[i,]$p1
      p2 <- peak_combinations[i,]$p2
      
      set_x <- compute_set_2o(
        peak_data_1$cov[[p1]], peak_data_2$cov[[p2]],
        peak_data_1$xopt[p1,], peak_data_2$xopt[p2,],
        resolution = set_resolution
      )
      set_y <- t(apply(set_x, 1, function(x) c(peak_data_1$peak_fns[[p1]](x),
                                               peak_data_2$peak_fns[[p2]](x))))
      
      colnames(set_x) <- paste0("x", 1:ncol(set_x))
      colnames(set_y) <- paste0("y", 1L:2L)
      
      cbind(set_x, set_y, p1 = p1, p2 = p2)
    })
    sets <- do.call(rbind, sets_individual)

    sets_nondominated <- sets[moleopt:::nondominated(sets[,c("y1", "y2")]),]
    y_nondominated <- as.matrix(sets_nondominated[,c("y1", "y2")])
    
    y_optimistic <- matrix(nrow = 0, ncol = 2)
    peak_combinations <- tibble(p1 = integer(), p2 = integer())
    
    for (set in sets_individual) {
      optimistic <- get_optimistic_y(set[,c("y1","y2")])
      optimistic_nondom <- moleopt:::nondominated(rbind(y_nondominated, optimistic))
      if (any(optimistic_nondom[-c(1:nrow(y_nondominated))])) {
        peak_combinations <- rbind(peak_combinations, set[1,c("p1", "p2")])
        y_optimistic <- rbind(y_optimistic, optimistic[optimistic_nondom[-c(1:nrow(y_nondominated))],])
      }
    }
    
    colnames(peak_combinations) <- c("p1", "p2")
    # print(nrow(peak_combinations))
    
    # hv_data <- hv_gap(y_nondominated, ref_point = c(1,1))
    # if (set_resolution >= min_final_resolution && gap - hv_data$gap < gap_target) break()
    # gap <- hv_data$gap
    # hv <- hv_data$hv
    
    hv <- ecr::computeHV(t(y_nondominated), ref.point = c(1,1))
    hv_optimistic <- ecr::computeHV(t(y_optimistic), ref.point = c(1,1))
    gap <- hv_optimistic - hv
    
    cat(sprintf("HV: %f (gap: %f) w/ resolution %d, sets remaining: %d\n", hv, gap, set_resolution, nrow(peak_combinations)))
    set_resolution <- 2 * set_resolution
  }
  
  list(
    hv = hv,
    gap = gap,
    set_resolution = set_resolution,
    n_relevant_sets = nrow(peak_combinations),
    sets_nondominated = sets_nondominated
  )
}

compute_hv_for_cutoffs <- function(logged_obj_space, cutoffs, ref_point = c(1, 1)) {
  sapply(cutoffs, function(c) {
    ecr::computeHV(logged_obj_space[, 1:c, drop = F], ref.point = ref_point)
  })
}
