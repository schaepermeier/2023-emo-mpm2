library(smoof)
library(moPLOT)
library(tidyverse)

source("01-utils.R")

design_viz <- expand_grid(
  topology = c("funnel", "random"), #, "mixed"
  seed = 1L,
  dimensions = c(2L),
  peaks = c(1L, 2L, 4L, 8L, 16L, 32L, 64L, 128L)
)

lapply(1:nrow(design_viz), function(i) {
  cat(paste0(i, "/", nrow(design_viz), "\r"))
  
  dimensions <- design_viz[i,]$dimensions
  peaks_1 <- design_viz[i,]$peaks
  peaks_2 <- design_viz[i,]$peaks
  topology_1 <- ifelse(design_viz[i,]$topology == "mixed", "funnel", design_viz[i,]$topology)
  topology_2 <- ifelse(design_viz[i,]$topology == "mixed", "random", design_viz[i,]$topology)
  seed_1 <- design_viz[i,]$seed
  seed_2 <- design_viz[i,]$seed + 1000L

  fn <- moPLOT::makeBiObjMPM2Function(
    dimensions = dimensions,
    n.peaks.1 = peaks_1, topology.1 = topology_1, seed.1 = seed_1,
    n.peaks.2 = peaks_2, topology.2 = topology_2, seed.2 = seed_2
  )
  
  ps <- getPLOT(fn)
  ggsave(paste0("figures/mpm2-2d-", peaks_1, "-", topology_1, "-plot-dec.png"), ps$plot_dec,
         width = unit(2, "in"), height = unit(2, "in"))
  ggsave(paste0("figures/mpm2-2d-", peaks_1, "-", topology_1, "-plot-obj.png"), ps$plot_obj,
         width = unit(2, "in"), height = unit(2, "in"))
})
