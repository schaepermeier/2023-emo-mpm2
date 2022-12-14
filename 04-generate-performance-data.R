library(smoof)
library(moPLOT)
library(tidyverse)

theme_set(theme_bw())

source("01-utils.R")

# Define experimental setup

algo_design <- expand_grid(
  algo_seed = 1L:15L,
  algorithm = c("nsga2", "omnioptimizer", "smsemoa", "mole", "random")
)

hv_df <- read_csv("03-generate-data.csv") %>% 
  dplyr::mutate(seed = as.integer(seed),
                dimensions = as.integer(dimensions),
                peaks = as.integer(peaks),
                set_resolution = as.integer(set_resolution),
                n_relevant_sets = as.integer(n_relevant_sets),
                n_points_nondominated = as.integer(n_points_nondominated))

exp_design <- hv_df %>%
  as_tibble()

nrow(exp_design)
nrow(exp_design) * nrow(algo_design)


# Cluster setup

cl <- parallel::makeCluster(4)
parallel::clusterEvalQ(cl, {
  library(smoof)
  library(moPLOT)
  library(tidyverse)
  library(ecr)

  source("01-utils.R")
  
  # required to load all auxiliary functions:
  moPLOT::makeBiObjMPM2Function()
})
parallel::clusterExport(cl, list("exp_design", "algo_design"))


# Evaluate problems

par_results <- parallel::parLapply(cl, 1:nrow(exp_design), function(i) {
  cat(paste0(i, "/", nrow(exp_design), "\n"))
  dimensions <- exp_design[i,]$dimensions
  peaks_1 <- exp_design[i,]$peaks
  peaks_2 <- exp_design[i,]$peaks
  topology_1 <- exp_design[i,]$topology # ifelse(exp_design[i,]$topology == "mixed", "funnel", exp_design[i,]$topology)
  topology_2 <- exp_design[i,]$topology # ifelse(exp_design[i,]$topology == "mixed", "random", exp_design[i,]$topology)
  seed_1 <- exp_design[i,]$seed
  seed_2 <- exp_design[i,]$seed + 1000L
  rotated_1 <- TRUE
  rotated_2 <- TRUE
  peak_shape_1 <- "ellipse"
  peak_shape_2 <- "ellipse"
  
  time_start <- Sys.time()
  
  filename <- paste0("results-2022-10-23/", dimensions, "d-", topology_1, "-", peaks_1, "p-", seed_1, ".csv")
  if (file.exists(filename)) return(0)
  
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
  
  fn <- smoof::makeMultiObjectiveFunction(
    name = paste0(topology_1, "_", peaks_1, "_", "instance_", seed_1, "_", dimensions, "D"),
    id = paste0(topology_1, "_", peaks_1, "_", "instance_", seed_1, "_", dimensions, "D"),
    fn = function(x) {c(
      peak_data_1$fn(x),
      peak_data_2$fn(x)
    )},
    par.set = ParamHelpers::makeNumericParamSet("x", lower = rep(-0.2, dimensions), upper = rep(1.2, dimensions))
  )
  
  algo_results <- lapply(1:nrow(algo_design), function(j) {
    algorithm <- algo_design[j,]$algorithm
    algo_seed <- algo_design[j,]$algo_seed
    
    set.seed(algo_seed)
    
    print(algorithm)
    
    SAMPLES <- 1e4L
    
    if (algorithm == "nsga2") {
      fn_wrapped <- smoof::addLoggingWrapper(fn, logg.y = TRUE)
      
      # ecr::nsga2(fn, lower = smoof::getLowerBoxConstraints(fn),
      #            upper = smoof::getUpperBoxConstraints(fn),
      #            terminators = list(stopOnEvals(SAMPLES)))
      
      mco::nsga2(fn_wrapped, idim = dimensions, odim = 2,
                 lower.bounds = smoof::getLowerBoxConstraints(fn),
                 upper.bounds = smoof::getUpperBoxConstraints(fn),
                 generations = (SAMPLES / 100L - 1L),
                 popsize = 100L)
      
      log <- smoof::getLoggedValues(fn_wrapped)
      obj_vals <- log$obj.vals
    } else if (algorithm == "smsemoa") {
      fn_wrapped <- smoof::addLoggingWrapper(fn, logg.y = TRUE)
      
      ecr::smsemoa(fn_wrapped, lower = smoof::getLowerBoxConstraints(fn),
                   upper = smoof::getUpperBoxConstraints(fn),
                   ref.point = c(1,1),
                   terminators = list(stopOnEvals(SAMPLES)))
      
      log <- smoof::getLoggedValues(fn_wrapped)
      obj_vals <- log$obj.vals
    } else if (algorithm == "random") {
      X <- flacco::createInitialSample(SAMPLES, dimensions, list("init_sample.lower" = smoof::getLowerBoxConstraints(fn),
                                                                 "init_sample.upper" = smoof::getUpperBoxConstraints(fn)))
      X <- matrix(X, nrow = nrow(X), ncol = ncol(X))
      obj_vals <- apply(X, 1, fn)
    } else if (algorithm == "mole") {
      fn_wrapped <- smoof::addLoggingWrapper(fn, logg.y = TRUE)
      
      starting_points <- lapply(1:(SAMPLES / (2 * dimensions + 1)), function(i) moleopt::runif_box(smoof::getLowerBoxConstraints(fn), smoof::getUpperBoxConstraints(fn))) %>% do.call(rbind, .)
      # starting_points <- rbind((smoof::getLowerBoxConstraints(fn) + smoof::getUpperBoxConstraints(fn)) / 2, starting_points)
      moleopt::run_mole(fn_wrapped, starting_points,
                        max_budget = SAMPLES, logging = "none",
                        refine_hv_target = 1e-3)
      
      log <- smoof::getLoggedValues(fn_wrapped)
      obj_vals <- log$obj.vals
    } else if (algorithm == "omnioptimizer") {
      fn_wrapped <- smoof::addLoggingWrapper(fn, logg.y = TRUE)
      pop_size <- 100L
      n_gens <- SAMPLES / pop_size
      
      omnioptr::omniopt(fn_wrapped, pop.size = pop_size, n.gens = n_gens, verbose = FALSE, seed = runif(1), frequency = 100L)
      
      log <- smoof::getLoggedValues(fn_wrapped)
      obj_vals <- log$obj.vals
    }
    
    cutoffs <- seq(100L, SAMPLES, by = 100L)
    
    hv_values <- compute_hv_for_cutoffs(
      obj_vals,
      cutoffs
    )
    
    hv_values_df <- rbind(hv_values)
    colnames(hv_values_df) <- paste0("t", 1:ncol(hv_values_df))
    
    cbind.data.frame(exp_design[i,], algorithm, algo_seed, hv_values_df)
  })
  
  results <- do.call(rbind, algo_results)
  
  write_csv(results, filename)
  
  as.numeric(Sys.time() - time_start)
})

parallel::stopCluster(cl)
