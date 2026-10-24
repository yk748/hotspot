##################################################################
# Supplementary hotspot simulation: Spurious Hotspot Detection
# Companion to Simulations_hotspot.R (Section 5.2 / Appendix B)
##################################################################
source("Hotspot_simulation_shared.R")

##################################################################
# 1. Global settings
##################################################################
n <- 100
iter_max <- 200
alpha <- 0.05
epsilon <- 0.2
eta <- 0.3
B_ci <- 1000
B_threshold <- 1000
G_values <- c(20, 40)
case_values <- c(1, 2, 4)
num_cores <- get_num_cores()

bi_thresholds <- calibrate_bi_thresholds(n, G_values, B_threshold, alpha)
refinements <- refinements_fixed

##################################################################
# 2. DGP: no true cross-series hotspot in any setting
#    none: both series iid N(0,1)
#    X_only / Y_only: just one series changes at t=50.
#    The changing series uses Cases 1, 2, and 4 in the main DGP.
##################################################################
run_spurious_replication <- function(scenario, case_id, n, G_values,
                                     alpha, epsilon, eta, B_ci) {
  if (scenario == "none") {
    series_X <- rnorm(n)
    series_Y <- rnorm(n)
    
  } else {
    pars <- case_par[case_par$case == case_id, ]
    changed <- generate_piecewise_series(n = n, true_cp = 50L, a = pars$a,
                                         sigma_min = pars$sigma_min, sigma_max = pars$sigma_max)
    unchanged <- rnorm(n)
    
    if (scenario == "X_only") {
      series_X <- changed
      series_Y <- unchanged
      
    } else if (scenario == "Y_only") {
      series_X <- unchanged
      series_Y <- changed
      
    } else stop("Unknown scenario: ", scenario)
  }
  fit_both_windows(series_X, series_Y, n, G_values,
                   alpha, epsilon, eta, B_ci)
}

##################################################################
# 3. Parallel Monte Carlo simulation
##################################################################
cl <- start_hotspot_cluster(num_cores, "run_spurious_replication")
summary_list <- list()
raw_file_list <- character()
counter <- 0L

for (scenario in c("none", "X_only", "Y_only")) {
  this_cases <- if (scenario == "none") 0L else case_values
  for (case_id in this_cases) {
    cat("\nRunning spurious setting: ", scenario, ", Case ", case_id,
        ", iter_max=", iter_max, "\n", sep = "")
    
    scenario_number <- match(scenario, c("none", "X_only", "Y_only"))
    scenario_seed <- 210000L + 1000L * scenario_number + 10L * case_id

    list_iter <- foreach::foreach(iter = seq_len(iter_max), .errorhandling = "pass",
                                  .options.RNG = scenario_seed) %dorng% {
                                    run_spurious_replication(scenario, case_id, n, G_values,
                                                             alpha, epsilon, eta, B_ci)
    }
    check_replications(list_iter, iter_max,paste(scenario, "case", case_id))

    raw_file <- paste0("spurious_raw_", scenario, "_case", case_id, ".RData")
    save(list_iter, scenario, case_id, n, iter_max, file = raw_file)
    raw_file_list <- c(raw_file_list, raw_file)

    for (method_name in names(method_map(G_values))) {
      result_name <- unname(method_map(G_values)[[method_name]])
      detected <- vapply(list_iter, function(res) {
        length(normalize_intervals(res[[result_name]], n)) > 0L
      }, logical(1))
      
      lengths <- vapply(list_iter, function(res) {
        hotspot_total_length(res[[result_name]], n)
      }, numeric(1))
      
      rule <- if (startsWith(method_name, "Thrs")) "Thrs" else "CI"
      G <- as.integer(sub(".*G", "", method_name))
      counter <- counter + 1L
      
      summary_list[[counter]] <- data.frame(
        Scenario = scenario, Case = if (scenario == "none") NA_integer_ else case_id,
        Rule = rule, G = G,
        FHR = mean(detected), FHL = mean(lengths))
    }
    cat("Done.\n")
  }
}
parallel::stopCluster(cl)

##################################################################
# 4. Final results (ONLY iter_max replications)
##################################################################
spurious_summary <- do.call(rbind, summary_list)
spurious_summary$Scenario <- factor(spurious_summary$Scenario,
                                    levels = c("none", "X_only", "Y_only"))
spurious_summary$Rule <- factor(spurious_summary$Rule,
                                levels = c("Thrs", "CI"))
spurious_summary <- spurious_summary[order(spurious_summary$Scenario, spurious_summary$Rule,
                                           spurious_summary$G, spurious_summary$Case), ]
row.names(spurious_summary) <- NULL

write.csv(spurious_summary, "spurious_hotspot_summary.csv", row.names = FALSE)
save(spurious_summary, case_par, bi_thresholds, raw_file_list,
     n, iter_max, alpha, epsilon, eta, B_ci, B_threshold,
     file = "spurious_hotspot_results.RData")
print(spurious_summary)
cat("Saved spurious_hotspot_summary.csv / spurious_hotspot_results.RData\n")
