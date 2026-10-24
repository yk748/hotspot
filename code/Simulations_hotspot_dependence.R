##################################################################
# Supplementary hotspot simulation: Sensitivity to Temporal Dependence
# Companion to Simulations_hotspot.R (Section 5.3 / Appendix B)
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
delta_values <- c(0, 5)
rho_values <- c(0.3, 0.6)
num_cores <- get_num_cores()

bi_thresholds <- calibrate_bi_thresholds(n, G_values, B_threshold, alpha)
refinements <- refinements_fixed

##################################################################
# 2. DGP: stationary AR(1) errors, separate X and Y processes
#    CP_X=50; CP_Y=50+delta; rho is common to the two series.
##################################################################
run_dependence_replication <- function(case_id, delta, rho, n, G_values,
                                       alpha, epsilon, eta, B_ci) {
  pars <- case_par[case_par$case == case_id, ]
  series_X <- generate_piecewise_ar1_series(n = n, true_cp = 50L, a = pars$a,
                                            sigma_min = pars$sigma_min, 
                                            sigma_max = pars$sigma_max, rho = rho)
  series_Y <- generate_piecewise_ar1_series(n = n, true_cp = 50L + delta, a = pars$a,
                                            sigma_min = pars$sigma_min, 
                                            sigma_max = pars$sigma_max, rho = rho)
  fit_both_windows(series_X, series_Y, n, G_values,alpha, epsilon, eta, B_ci)
}

##################################################################
# 3. Parallel Monte Carlo simulation
##################################################################
cl <- start_hotspot_cluster(num_cores, "run_dependence_replication")
summary_list <- list()
raw_file_list <- character()
counter <- 0L

for (rho in rho_values) {
  for (delta in delta_values) {
    for (case_id in case_values) {
      
      cat("\nRunning temporal dependence: rho=", rho,
          ", delta=", delta, ", Case ", case_id,
          ", iter_max=", iter_max, "\n", sep = "")
      scenario_seed <- 310000L + 10000L * as.integer(10 * rho) + 100L * as.integer(delta) + case_id
      
      list_iter <- foreach::foreach(iter = seq_len(iter_max), .errorhandling = "pass",
                                    .options.RNG = scenario_seed) %dorng% {
        run_dependence_replication(case_id, delta, rho, n, G_values,alpha, epsilon, eta, B_ci)
      }
      check_replications(list_iter, iter_max,
                         paste("rho", rho, "delta", delta, "case", case_id))

      truth <- list(X = 50L, Y = 50L + delta)
      raw_file <- paste0("dependence_raw_rho", gsub("\\.", "p", as.character(rho)),
                         "_delta", delta, "_case", case_id, ".RData")
      save(list_iter, truth, rho, delta, case_id, n, iter_max, file = raw_file)
      raw_file_list <- c(raw_file_list, raw_file)

      for (method_name in names(method_map(G_values))) {
        
        result_name <- unname(method_map(G_values)[[method_name]])
        metrics <- t(vapply(list_iter, function(res) {
          evaluate_hotspot(res[[result_name]], truth$Y, n)
        }, numeric(2)))
        
        rule <- if (startsWith(method_name, "Thrs")) "Thrs" else "CI"
        
        G <- as.integer(sub(".*G", "", method_name))
        counter <- counter + 1L
        
        summary_list[[counter]] <- data.frame(
          Rho = rho, Delta = delta, Case = case_id, Rule = rule, G = G,
          HitRate = mean(metrics[, "HR"]),
          HotspotLength = mean(metrics[, "HotspotLength"]))
      }
      cat("Done.\n")
    }
  }
}
parallel::stopCluster(cl)

##################################################################
# 4. Final results (ONLY iter_max replications)
##################################################################
dependence_summary <- do.call(rbind, summary_list)
dependence_summary$Rule <- factor(dependence_summary$Rule,
                                  levels = c("Thrs", "CI"))
dependence_summary <- dependence_summary[order(dependence_summary$Rho, 
                                               dependence_summary$Delta,
                                               dependence_summary$Rule, 
                                               dependence_summary$G, dependence_summary$Case), ]
row.names(dependence_summary) <- NULL
write.csv(dependence_summary, "dependence_hotspot_summary.csv", row.names = FALSE)
save(dependence_summary, case_par, bi_thresholds, raw_file_list,
     n, iter_max, alpha, epsilon, eta, B_ci, B_threshold,
     file = "dependence_hotspot_results.RData")
print(dependence_summary)
cat("Saved dependence_hotspot_summary.csv / dependence_hotspot_results.RData\n")

