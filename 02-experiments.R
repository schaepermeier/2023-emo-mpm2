library(smoof)
library(moPLOT)
library(tidyverse)

theme_set(theme_bw())

source("01-utils.R")

# ==== Problem Definition ====

dimensions <- 40L
peaks_1 <- 256L
peaks_2 <- 256L
topology_1 <- "random"
topology_2 <- "random"
seed_1 <- 1L # 667172L
seed_2 <- 1001L # 540835L
rotated_1 <- TRUE
rotated_2 <- TRUE
peak_shape_1 <- "ellipse"
peak_shape_2 <- "ellipse"

fn <- moPLOT::makeBiObjMPM2Function(
  dimensions = dimensions,
  n.peaks.1 = peaks_1, topology.1 = topology_1, seed.1 = seed_1,
  n.peaks.2 = peaks_2, topology.2 = topology_2, seed.2 = seed_2
)

if (dimensions == 2L) {
  # ps <- getPLOT(fn, resolution = 500L)
  # gridExtra::grid.arrange(ps$plot_dec, ps$plot_obj, ncol = 2)
}

# ==== Extraction of optima, cov for first objective ====

peak_data_1 <- get_peak_metadata(
  dimension = dimensions,
  npeaks = peaks_1,
  topology = topology_1,
  randomSeed = seed_1,
  rotated = rotated_1,
  peakShape = peak_shape_1
)
all_peaks_1 <- nrow(peak_data_1$xopt)

# ==== Extraction of optima, cov for second objective ====

peak_data_2 <- get_peak_metadata(
  dimension = dimensions,
  npeaks = peaks_2,
  topology = topology_2,
  randomSeed = seed_2,
  rotated = rotated_2,
  peakShape = peak_shape_2
)
all_peaks_2 <- nrow(peak_data_2$xopt)

# ==== Joint visualization with optima ====

if (dimensions == 2L) {
  gridExtra::grid.arrange(
    ps$plot_dec +
      geom_point(mapping = aes(x1, x2), data = as.data.frame(peak_data_1$xopt)) +
      geom_point(mapping = aes(x1, x2), data = as.data.frame(peak_data_2$xopt)),
    ps$plot_obj,
    ncol = 2
  )
}

# ==== Extract theoretical sets ====

set_resolution <- 100L
peak_combinations <- expand_grid(p1 = 1:all_peaks_1, p2 = 1:all_peaks_2)

sets <- lapply(1:nrow(peak_combinations), function(i) {
  p1 <- peak_combinations[i,]$p1
  p2 <- peak_combinations[i,]$p2
  
  set_x <- compute_set_2o(
    peak_data_1$cov[[p1]], peak_data_2$cov[[p2]],
    peak_data_1$xopt[p1,], peak_data_2$xopt[p2,],
    resolution = set_resolution
  )
  
  set_y <- t(apply(set_x, 1, function(x) c(peak_data_1$peak_fns[[p1]](x),
                                           peak_data_2$peak_fns[[p2]](x))))
  
  colnames(set_x) <- paste0("x", 1:dimensions)
  colnames(set_y) <- paste0("y", 1L:2L)
  
  cbind(set_x, set_y, p1 = p1, p2 = p2)
}) %>% do.call(rbind, .)

set_alpha <- ifelse(moleopt:::nondominated(sets[,c("y1", "y2")]), 1, 0.5)

if (dimensions == 2L) {
  gridExtra::grid.arrange(
    ps$plot_dec +
      geom_line(mapping = aes(x1, x2, group = paste0(p1, ", ", p2)), alpha = set_alpha, data = as.data.frame(sets)),
    ps$plot_obj +
      geom_line(mapping = aes(y1, y2, group = paste0(p1, ", ", p2)), alpha = set_alpha, data = as.data.frame(sets)),
    ncol = 2
  )
} else {
  ggplot(as.data.frame(sets), aes(y1, y2, group = paste(p1, ", ", p2))) +
    geom_path() +
    theme_bw() +
    coord_fixed() +
    theme(legend.position = "none") +
    lims(x = c(0,1), y = c(0,1))
}

# ==== Approximate Hypervolume ====

hv_data <- approximate_hv(peak_data_1, peak_data_2, initial_resolution = 4L, gap_target = 1e-03)
sets_nondominated <- hv_data$sets_nondominated

if (dimensions == 2L) {
  gridExtra::grid.arrange(
    ps$plot_dec +
      geom_point(mapping = aes(x1, x2), data = as.data.frame(sets_nondominated)),
    ps$plot_obj +
      geom_point(mapping = aes(y1, y2), data = as.data.frame(sets_nondominated)),
    ncol = 2
  )
} else {
  # ggplot(as.data.frame(sets_nondominated), aes(x = y1, y = y2, color = factor(paste(p1, p2)))) +
  #   geom_point() +
  #   theme_bw() +
  #   coord_fixed() +
  #   theme(legend.position = "none") +
  #   lims(x = c(0,1), y = c(0,1))
  ggplot(as.data.frame(sets), aes(y1, y2)) +
    geom_path(aes(group = paste(p1, ", ", p2)), color = "lightgray") +
    geom_point(data = as.data.frame(sets_nondominated), size = 0.05) +
    theme_bw() +
    coord_fixed() +
    theme(legend.position = "none") +
    scale_x_continuous(limits = c(0, 1), labels = c(0, "", "", "", 1)) +
    scale_y_continuous(limits = c(0, 1), labels = c(0, "", "", "", 1)) +
    labs(x = expression(y[1]),
         y = expression(y[2]))
  # ggsave(paste0("figures/sets-", topology_1, "-", dimensions, "d.pdf"),
  #        width = 2, height = 2)
}

