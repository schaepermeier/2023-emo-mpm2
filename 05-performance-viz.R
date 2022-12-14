library(smoof)
library(moPLOT)
library(tidyverse)
library(scmamp) # devtools::install_github("b0rxa/scmamp")

theme_set(theme_bw())

source("01-utils.R")

hv_df <- read_csv("03-generate-data.csv") %>% 
  dplyr::mutate(seed = as.integer(seed),
                dimensions = as.integer(dimensions),
                peaks = as.integer(peaks),
                set_resolution = as.integer(set_resolution),
                n_relevant_sets = as.integer(n_relevant_sets),
                n_points_nondominated = as.integer(n_points_nondominated))

exp_results <- read_csv(paste0("results/", list.files("results/"))) %>% 
  select(-c(hv, gap, set_resolution, n_relevant_sets, n_points_nondominated, xmin, xmax)) %>%
  left_join(hv_df)


# ==== Critical Differences ====

for (d in c(2, 3, 5, 10)) {
  for (t in c("random", "funnel")) {
    # for (p in c(1, 2, 4, 8, 16, 32, 64, 128)) {
      final_performance <- exp_results %>%
        filter(dimensions == d & topology == t) %>%
        select(topology, seed, dimensions, peaks, algorithm, algo_seed, hv, final_hv = t100) %>% 
        group_by(topology, seed, dimensions, peaks, algorithm) %>%
        summarize(mean_final_hv = -mean(hv - final_hv)) %>%
        ungroup() %>%
        pivot_wider(names_from = algorithm, values_from = mean_final_hv) %>% 
        select(nsga2, omnioptimizer, smsemoa, mole, random)
      
      # pdf(paste0("figures/cd-mean-final-hv-", t, "-", d, "d.pdf"), width = 6, height = 1.2)
      scmamp::plotCD(final_performance, cex = 0.87)
      # dev.off()
    # }
  }
}





 exp_results %>%
  filter(peaks %in% c(2, 8, 32, 128) & dimensions == 5) %>% 
  pivot_longer(paste0("t", 1:100), names_to = "t", values_to = "achieved_hv") %>% 
  mutate(t = 100 * as.integer(substring(t, 2)),
                delta_to_target = (hv + gap - achieved_hv) / (hv + gap),
                log_to_target = ifelse(delta_to_target > 0, log10(delta_to_target), log10(min(gap)))) %>% 
  group_by(topology, dimensions, algorithm, peaks, t) %>%
  summarize(mean_to_target = mean(delta_to_target), se_to_target = sd(delta_to_target) / sqrt(n()),
            min_to_target = min(delta_to_target), max_to_target = max(delta_to_target),
            median_to_target = median(delta_to_target)) %>% 
  ggplot(aes(t, mean_to_target, color = algorithm, fill = algorithm)) +
  facet_grid(vars(topology), vars(peaks), scales = "free") +
  scale_x_log10(labels = function(x) format(x, scientific = TRUE)) +
  scale_y_log10() +
  # lims(y = c(0, 1)) +
  labs(
    x = "Function Evaluations",
    y = "Mean Relative Hypervolume Gap",
    color = "Algorithm"
  ) +
  scale_colour_discrete(
    labels = c(
      "mole" = "MOLE",
      "nsga2" = "NSGA-II",
      "omnioptimizer" = "Omni-Optimizer",
      "random" = "Random Search",
      "smsemoa" = "SMS-EMOA"
    )
  ) +
  geom_step() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

ggsave("figures/algorithms-convergence-plot.pdf", width = 8, height = 3.5)

exp_results %>%
  pivot_longer(paste0("t", 1:100), names_to = "t", values_to = "achieved_hv") %>% 
  dplyr::mutate(t = as.integer(substring(t, 2))) %>% 
  group_by(topology, seed, dimensions, algorithm, peaks, t) %>% 
  summarize(mean_hv = mean(achieved_hv), hv = mean(hv)) %>% 
  filter(peaks == 16L & dimensions == 5L) %>% 
  ggplot(aes(t, mean_hv, color = algorithm)) +
  facet_grid(vars(dimensions, topology), vars(peaks, seed)) +
  geom_abline(aes(intercept = hv, slope = 0)) +
  geom_step()

exp_results %>%
  pivot_longer(paste0("t", 1:100), names_to = "t", values_to = "achieved_hv") %>% 
  dplyr::mutate(t = as.integer(substring(t, 2)),
                delta_to_target = (hv + gap - achieved_hv) / (hv + gap),
                log_to_target = ifelse(delta_to_target > 0, log10(delta_to_target), -4)) %>%
  filter(t == max(t)) %>% 
  group_by(dimensions, n_relevant_sets, algorithm) %>%
  summarize(mean_log_to_target = mean(log_to_target), sd_log_to_target = sd(log_to_target)) %>%
  ggplot(aes(factor(n_relevant_sets), mean_log_to_target, color = algorithm, fill = algorithm)) +
  facet_grid(vars(dimensions)) +
  geom_point()



# ==== ERT Analysis ====

ABS_TARGET <- 0.01

ert_df <- exp_results %>%
  pivot_longer(paste0("t", 1:100), names_to = "t", values_to = "achieved_hv") %>% 
  dplyr::mutate(t = 100 * as.integer(substring(t, 2)),
                delta_to_target = (hv + gap - achieved_hv)) %>%
  group_by(topology, seed, algo_seed, dimensions, algorithm, peaks) %>% 
  summarize(
    solved = (sum(delta_to_target < ABS_TARGET) > 0),
    target_hit_time = ifelse(solved, min(t[delta_to_target < ABS_TARGET]), max(t))
  ) %>% 
  group_by(topology, seed, dimensions, algorithm, peaks) %>% 
  summarize(mean_rt = mean(target_hit_time),
            p_solved = mean(solved),
            ert = mean_rt / p_solved)

PENALTY <- 1 * 15 * 10000 # penalty factor * repetitions * evals per run

ert_df_wide <- ert_df %>% 
  mutate(ert = pmin(ert, PENALTY)) %>% 
  group_by(topology, dimensions, peaks, algorithm) %>% 
  summarize(mean_ert = mean(ert)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = algorithm, values_from = mean_ert)

for (TOPOLOGY in c("random", "funnel")) {
  for (DIMENSIONS in c(2L, 3L, 5L, 10L)) {
    # TOPOLOGY <- "random"
    # DIMENSIONS <- 2L
    
    g_ert <- ert_df %>% 
      mutate(ert = pmin(ert, PENALTY)) %>% 
      group_by(topology, dimensions, peaks, algorithm) %>% 
      summarize(mean_ert = mean(ert), se_ert = sd(ert) / sqrt(n())) %>% 
      ungroup() %>% 
      filter(dimensions == DIMENSIONS & topology == TOPOLOGY) %>%
      ggplot(aes(peaks, mean_ert, color = algorithm)) +
      scale_x_continuous(trans = "log2") +
      scale_y_log10(limits = c(10, PENALTY)) +
      scale_colour_discrete(
        labels = c(
          "mole" = "MOLE",
          "nsga2" = "NSGA-II",
          "omnioptimizer" = "Omni-Optimizer",
          "random" = "Random Search",
          "smsemoa" = "SMS-EMOA"
        )
      ) +
      scale_fill_discrete(guide = "none") +
      geom_line() +
      geom_ribbon(aes(ymin = mean_ert - se_ert,
                      ymax = mean_ert + se_ert,
                      fill = algorithm),
                  alpha = 0.5,
                  color = NA) +
      labs(
        x = "Number of Peaks",
        y = "Mean ERT",
        color = "Algorithm"
      )
    
    ggsave(
      paste0("figures/ert-peaks-", TOPOLOGY, "-", DIMENSIONS, "d.pdf"), g_ert,
        width = unit(5, "cm"), height = unit(3, "cm")
      )
  }
}

