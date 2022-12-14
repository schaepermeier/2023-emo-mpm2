library(smoof)
library(moPLOT)
library(tidyverse)

source("01-utils.R")

fn_for_setup_only <- moPLOT::makeBiObjMPM2Function()

theme_set(theme_bw())

mpm2_design <- expand_grid(
  topology = c("funnel", "random"), #, "mixed"
  seed = 1L:20L,
  dimensions = c(2L, 3L, 5L, 10L),
  peaks = c(1L, 2L, 4L, 8L, 16L, 32L, 64L, 128L)
)

nrow(mpm2_design)

results <- lapply(1:nrow(mpm2_design), function(i) {
  cat(paste0(i, "/", nrow(mpm2_design), "\r"))
  
  dimensions <- mpm2_design[i,]$dimensions
  peaks_1 <- mpm2_design[i,]$peaks
  peaks_2 <- mpm2_design[i,]$peaks
  topology_1 <- ifelse(mpm2_design[i,]$topology == "mixed", "funnel", mpm2_design[i,]$topology)
  topology_2 <- ifelse(mpm2_design[i,]$topology == "mixed", "random", mpm2_design[i,]$topology)
  seed_1 <- mpm2_design[i,]$seed
  seed_2 <- mpm2_design[i,]$seed + 1000L
  rotated_1 <- TRUE
  rotated_2 <- TRUE
  peak_shape_1 <- "ellipse"
  peak_shape_2 <- "ellipse"
  
  peak_data_1 <- get_peak_metadata(
    dimension = dimensions,
    npeaks = peaks_1,
    topology = topology_1,
    randomSeed = seed_1,
    rotated = rotated_1,
    peakShape = peak_shape_1
  )

  peak_data_2 <- get_peak_metadata(
    dimension = dimensions,
    npeaks = peaks_2,
    topology = topology_2,
    randomSeed = seed_2,
    rotated = rotated_2,
    peakShape = peak_shape_2
  )

  hv_data <- approximate_hv(peak_data_1, peak_data_2, initial_resolution = 4L, gap_target = 1e-4)
  
  cbind(
    mpm2_design[i,],
    hv = hv_data$hv,
    gap = hv_data$gap,
    set_resolution = hv_data$set_resolution,
    n_relevant_sets = hv_data$n_relevant_sets,
    n_points_nondominated = nrow(hv_data$sets_nondominated),
    xmax = max(select(as_tibble(hv_data$sets_nondominated), starts_with("x"))),
    xmin = min(select(as_tibble(hv_data$sets_nondominated), starts_with("x")))
  )
})

results_df <- results %>% Reduce(rbind, .)

write_csv(results_df, "03-generate-data.csv")

results_df <- read_csv("03-generate-data.csv")

ggplot(results_df, aes(peaks, n_relevant_sets)) +
  geom_boxplot(aes(group = factor(peaks))) +
  facet_grid(rows = vars(dimensions), cols = vars(topology)) +
  labs(x = "Peaks",
       y = "Number of Relevant Local Sets") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(trans = "log2") +
  scale_color_viridis_d()

ggsave("figures/problems-relevant-sets.pdf", width = 4, height = 4)

ggplot(results_df, aes(factor(peaks), n_relevant_sets)) +
  facet_grid(rows = vars(dimensions)) +
  labs(x = "Peaks",
       y = "Number of local sets in Pareto set") +
  geom_violin(aes(color = topology))

ggplot(results_df, aes(peaks, hv)) +
  geom_boxplot(aes(group = factor(peaks))) +
  facet_grid(rows = vars(dimensions), cols = vars(topology)) +
  labs(x = "Peaks",
       y = "Hypervolume of Pareto Set") +
  theme(legend.position = "bottom") +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(limits = c(0, 1), labels = c(0, "", "", "", 1)) +
  scale_color_viridis_d()

ggsave("figures/problems-hv.pdf", width = 4, height = 4)

