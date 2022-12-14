library(smoof)
library(moPLOT)
library(tidyverse)

theme_set(theme_bw())

source("01-utils.R")

# ==== Problem Definition ====

dimensions <- 2L
peaks_1 <- 2L
peaks_2 <- 2L
topology_1 <- "random"
topology_2 <- "random"
seed_1 <- 667172L # as.integer(runif(1, min = 0, max = 1e6)) # 4L # 667172L
seed_2 <- 540835L # as.integer(runif(1, min = 0, max = 1e6)) # 8L # 540835L
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
  ps <- getPLOT(fn, resolution = 500L)
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
  ggplot(as.data.frame(sets), aes(y1, y2, alpha = set_alpha)) +
    geom_point() +
    theme_bw() +
    coord_fixed() +
    theme(legend.position = "none") +
    lims(x = c(0,1), y = c(0,1))
}

# ==== Visualization per Peak ====

design <- moPLOT::generateDesign(fn, points.per.dimension = 500L)

for (p_1 in 1:2) {
  for (p_2 in 1:2) {
    # p_1 <- 1L
    # p_2 <- 1L
    
    y1 <- apply(design$dec.space, 1, peak_data_1$peak_fns[[p_1]])
    y2 <- apply(design$dec.space, 1, peak_data_2$peak_fns[[p_2]])
    
    single_set_data <- as.data.frame(sets) %>%
      filter(p1 == p_1 & p2 == p_2)
    
    g_dec <- ggplot(single_set_data, aes(x1, x2)) +
      geom_contour(aes(x1, x2, z = y1), color = "darkgray", breaks = quantile(y1, seq(0, 1, by = 0.1)**2),
                   data = cbind.data.frame(design$dec.space, y1 = y1)) +
      geom_contour(aes(x1, x2, z = y2), color = "lightgray", breaks = quantile(y2, seq(0, 1, by = 0.1)**2),
                   data = cbind.data.frame(design$dec.space, y2 = y2)) +
      geom_path(aes(color = paste0(p_1, ", ", p_2))) +
      scale_color_discrete(limits = c("1, 1", "1, 2", "2, 1", "2, 2")) +
      scale_x_continuous(limits = c(0, 1), labels = c(0, "", "", "", 1)) +
      scale_y_continuous(limits = c(0, 1), labels = c(0, "", "", "", 1)) +
      labs(x = expression(x[1]),
           y = expression(x[2])) +
      theme(legend.position = "none") +
      coord_fixed()
    
    g_obj <- ggplot(single_set_data, aes(y1, y2)) +
      geom_path(aes(color = paste0(p_1, ", ", p_2))) +
      scale_color_discrete(limits = c("1, 1", "1, 2", "2, 1", "2, 2")) +
      scale_x_continuous(limits = c(0, 1), labels = c(0, "", "", "", 1)) +
      scale_y_continuous(limits = c(0, 1), labels = c(0, "", "", "", 1)) +
      labs(x = expression(y[1]),
           y = expression(y[2])) +
      theme(legend.position = "none") +
      coord_fixed()
    
    ggsave(
      paste0("figures/bipeak-example-", p_1, "-", p_2, "-dec.pdf"), g_dec,
      width = 2, height = 2
    )
    
    ggsave(
      paste0("figures/bipeak-example-", p_1, "-", p_2, "-obj.pdf"), g_obj,
      width = 2, height = 2
    )
  }
}

set_y_allpeaks <- lapply(1:nrow(sets), function(i) {
  x1 <- sets[i,"x1"]
  x2 <- sets[i,"x2"]
  
  fn(c(x1, x2))
}) %>% do.call(rbind, .)

set_active <- lapply(1:nrow(sets), function(i) {
  y1 <- sets[i,"y1"]
  y2 <- sets[i,"y2"]
  x1 <- sets[i,"x1"]
  x2 <- sets[i,"x2"]
  
  y_allpeaks <- fn(c(x1, x2))
  
  ecr::dominates(c(y1, y2), y_allpeaks + 1e-8)
}) %>% do.call(rbind, .)

# g_bipeak_dec <- ggplot(cbind.data.frame(design$dec.space, design$obj.space), aes(x1, x2)) +
#   geom_contour(aes(z = y1), color = "red", breaks = quantile(design$obj.space[,"y1"], seq(0, 1, by = 0.1)**2)) +
#   geom_contour(aes(z = y2), color = "blue", breaks = quantile(design$obj.space[,"y2"], seq(0, 1, by = 0.1)**2)) +
#   lims(x = c(0, 1),
#        y = c(0, 1)) +
#   labs(x = expression(x[1]),
#        y = expression(x[2])) +
#   geom_path(mapping = aes(x1, x2, group = paste0(p1, ", ", p2, set_active)),
#             linetype = ifelse(set_active, "solid", "dashed"), data = as.data.frame(sets)) +
#   coord_fixed()
# g_bipeak_dec

g_bipeak_dec <- ggplot(cbind.data.frame(design$dec.space, design$obj.space), aes(x1, x2)) +
  geom_contour(aes(z = y1), color = "darkgray", breaks = quantile(design$obj.space[,"y1"], seq(0, 1, by = 0.1)**2)) +
  geom_contour(aes(z = y2), color = "lightgray", breaks = quantile(design$obj.space[,"y2"], seq(0, 1, by = 0.1)**2)) +
  scale_x_continuous(limits = c(0, 1), labels = c(0, "", "", "", 1)) +
  scale_y_continuous(limits = c(0, 1), labels = c(0, "", "", "", 1)) +
  labs(x = expression(x[1]),
       y = expression(x[2]),
       color = "Peaks") +
  geom_path(mapping = aes(x1, x2, group = paste0(p1, ", ", p2, set_active), color = paste0(p1, ", ", p2)),
            linetype = ifelse(set_active, "solid", "dashed"), data = as.data.frame(sets)) +
  theme(legend.position = "none") +
  scale_color_discrete(limits = c("1, 1", "1, 2", "2, 1", "2, 2")) +
  coord_fixed()
g_bipeak_dec

g_bipeak_obj <- ggplot(cbind.data.frame(design$dec.space, design$obj.space), aes(y1, y2)) +
  scale_x_continuous(limits = c(0, 1), labels = c(0, "", "", "", 1)) +
  scale_y_continuous(limits = c(0, 1), labels = c(0, "", "", "", 1)) +
  labs(x = expression(y[1]),
       y = expression(y[2]),
       color = "Peaks") +
  geom_path(mapping = aes(y1, y2, group = paste0(p1, ", ", p2), color = paste0(p1, ", ", p2)),
            alpha = ifelse(set_active, 1, 0.5), data = as.data.frame(sets)) +
  geom_path(mapping = aes(set_y_allpeaks[,1], set_y_allpeaks[,2], group = paste0(p1, ", ", p2, set_active), color = paste0(p1, ", ", p2)),
            linetype = ifelse(set_active, "solid", "dashed"), data = as.data.frame(sets)) +
  scale_color_discrete(limits = c("1, 1", "1, 2", "2, 1", "2, 2")) +
  theme(legend.position = "none") +
  coord_fixed()
g_bipeak_obj

ggsave(
  paste0("figures/bipeak-example-all-dec.pdf"), g_bipeak_dec,
  width = 2, height = 2
)

ggsave(
  paste0("figures/bipeak-example-all-obj.pdf"), g_bipeak_obj,
  width = 2, height = 2
)

ggsave(
  paste0("figures/bipeak-example-plot-dec.pdf"), ps$plot_dec  +
    scale_x_continuous(limits = c(0, 1), labels = c(0, "", "", "", 1)) +
    scale_y_continuous(limits = c(0, 1), labels = c(0, "", "", "", 1)),
  width = 2, height = 2
)

ggsave(
  paste0("figures/bipeak-example-plot-obj.png"), ps$plot_obj  +
    scale_x_continuous(limits = c(0, 1), labels = c(0, "", "", "", 1)) +
    scale_y_continuous(limits = c(0, 1), labels = c(0, "", "", "", 1)),
  width = 2, height = 2
)

# ==== Approximate Hypervolume ====

gap <- Inf
hv <- 0
initial_resolution <- 4L
gap_target <- 1e-4

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
  sets_combined <- do.call(rbind, sets_individual)
  y_hypervolume <- lapply(sets_individual, function(set) {
    new_set <- rbind(set[,c("y1", "y2")], get_pessimistic_y(set[,c("y1", "y2")]) - 1e-6)
    new_set <- rbind(new_set, c(max(set[,"y1"]), min(set[,"y2"])), c(min(set[,"y1"]), max(set[,"y2"])), c(1, 0))
    new_set <- cbind(new_set, p1 = set[1, "p1"], p2 = set[2, "p2"])
    new_set
  }) %>% do.call(rbind, .)
  y_optimistic_hypervolume <- lapply(sets_individual, function(set) {
    new_set <- rbind(set[,c("y1", "y2")], get_optimistic_y(set[,c("y1", "y2")]) + 1e-6)
    new_set <- rbind(new_set, c(max(set[,"y1"]), min(set[,"y2"])), c(min(set[,"y1"]), max(set[,"y2"])), c(1, 0))
    new_set <- cbind(new_set, p1 = set[1, "p1"], p2 = set[2, "p2"])
    new_set
  }) %>% do.call(rbind, .)
  
  ggplot(as.data.frame(sets_combined) %>% filter(p1 == 1 & p2 == 2),
         aes(y1, y2,
             color = paste0(p1, ", ", p2),
             fill = paste0(p1, ", ", p2))) +
    # geom_point(data = as.data.frame(y_optimistic_hypervolume) %>% filter(p1 == 1 & p2 == 2),
    #            color = "black") +
    geom_path(data = as.data.frame(sets) %>% filter(p1 == 1 & p2 == 2)) +
    geom_ribbon(data = as.data.frame(y_hypervolume) %>% filter(p1 == 1 & p2 == 2),
                mapping = aes(ymin = y2, ymax = 1),
                alpha = 0.5, color = "transparent") +
    geom_ribbon(data = as.data.frame(y_optimistic_hypervolume) %>% filter(p1 == 1 & p2 == 2),
                mapping = aes(ymin = y2, ymax = 1),
                alpha = 0.5, color = "transparent") +
    geom_step(data = as.data.frame(sets_combined) %>% filter(p1 == 1 & p2 == 2),
              direction = "vh", color = "black", linetype = "dashed", alpha = 1) +
    geom_point() +
    scale_x_continuous(limits = c(0, 1), labels = c(0, "", "", "", 1)) +
    scale_y_continuous(limits = c(0, 1), labels = c(0, "", "", "", 1)) +
    labs(x = expression(y[1]),
         y = expression(y[2])) +
    scale_color_discrete(limits = c("1, 1", "1, 2", "2, 1", "2, 2")) +
    scale_fill_discrete(limits = c("1, 1", "1, 2", "2, 1", "2, 2")) +
    coord_fixed() +
    theme(legend.position = "none")
  
  ggsave("figures/method-gap-illustration.pdf", width = 3, height = 3)
  
  sets_nondominated <- sets_combined[moleopt:::nondominated(sets_combined[,c("y1", "y2")]),]
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
  
  y_optimistic <- y_optimistic[moleopt:::nondominated(y_optimistic),]
  y_optimistic <- y_optimistic[order(y_optimistic[,1]),]
  colnames(y_optimistic) <- paste0("y", 1:2)
  
  y_with_creators <- cbind(y_optimistic + 1e-6, p1 = c(2, 2, 1, 2), p2 = c(2, 2, 2, 2))
  y_with_creators <- rbind(
    y_with_creators,
    sets_individual[[2]][2:3, c("y1", "y2", "p1", "p2")],
    sets_individual[[4]][,c("y1", "y2", "p1", "p2")]
  )
  
  ggplot(as.data.frame(sets_combined),
         aes(y1, y2,
             color = paste0(p1, ", ", p2),
             fill = paste0(p1, ", ", p2))) +
    geom_ribbon(data = as.data.frame(y_optimistic_hypervolume),
                mapping = aes(ymin = y2, ymax = 1, group = paste0(p1, ", ", p2)),
                alpha = 1, color = "transparent", fill = "gray") +
    geom_ribbon(data = as.data.frame(y_hypervolume),
                mapping = aes(ymin = y2, ymax = 1, group = paste0(p1, ", ", p2)),
                alpha = 1, color = "transparent", fill = "black") +
    geom_step(data = as.data.frame(y_with_creators),
              direction = "vh", linetype = "dashed") +
    geom_point(data = as.data.frame(y_optimistic), aes(y1, y2), color = "black", fill = "transparent") +
    geom_point() +
    scale_x_continuous(limits = c(0, 1), labels = c(0, "", "", "", 1)) +
    scale_y_continuous(limits = c(0, 1), labels = c(0, "", "", "", 1)) +
    labs(x = expression(y[1]),
         y = expression(y[2])) +
    scale_color_discrete(limits = c("1, 1", "1, 2", "2, 1", "2, 2")) +
    scale_fill_discrete(limits = c("1, 1", "1, 2", "2, 1", "2, 2")) +
    coord_fixed() +
    theme(legend.position = "none")
  
  ggsave("figures/method-global-gap-4.pdf", width = 3, height = 3)
  
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


