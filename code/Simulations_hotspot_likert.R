##################################################################
# Supplementary hotspot simulation: Randomized Likert Transformation
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
M <- 10L                     # randomized transformations per dataset
num_cores <- get_num_cores()

bi_thresholds <- calibrate_bi_thresholds(n, G_values, B_threshold, alpha)
refinements <- refinements_fixed

##################################################################
# 2. DGP: X is continuous; Y* is continuous and converted to 1--5.
#    Within each replication, X and ordinal Y remain fixed for all
#    M transformations. Only W_t^(m) changes across transformations.
#    For CI bootstrap, each G uses an identical random seed for all M,
#    reducing bootstrap randomness in comparisons of transformations.
##################################################################
run_likert_replication <- function(case_id, delta, n, G_values, M,
                                   alpha, epsilon, eta, B_ci) {
  
  pars <- case_par[case_par$case == case_id, ]
  series_X <- generate_piecewise_series(n = n, true_cp = 50L, a = pars$a,
                                        sigma_min = pars$sigma_min, sigma_max = pars$sigma_max)
  latent_Y <- generate_piecewise_series(n = n, true_cp = 50L + delta, a = pars$a,
                                        sigma_min = pars$sigma_min, sigma_max = pars$sigma_max)
  cutpoints <- c(-Inf, qnorm((1:4) / 5), Inf)
  ordinal_Y <- as.integer(cut(latent_Y, breaks = cutpoints, labels = 1:5, include.lowest = TRUE))

  # Generate independent Uniform(0,1) draws up front.  A fixed CI
  # bootstrap seed per G helps separate transformation sensitivity
  # from arbitrary changes in bootstrap RNG across repetitions.
  W <- matrix(runif(n * M), nrow = n, ncol = M)
  bootstrap_seeds <- sample.int(100000000L, length(G_values))
  names(bootstrap_seeds) <- as.character(G_values)

  all_hotspots <- setNames(vector("list", length(method_map(G_values))),
                           unname(method_map(G_values)))
  for (name in names(all_hotspots)) {
    all_hotspots[[name]] <- vector("list", M)
  }
  
  for (m in seq_len(M)) {
    series_Y <- transform_ordinal_with_uniforms(ordinal_Y, W[, m])
    
    for (G in G_values) {
      
      set.seed(bootstrap_seeds[[as.character(G)]])
      fit <- construct_hotspots(series_X, series_Y, G = G, n = n,
                                alpha = alpha, epsilon = epsilon,
                                eta = eta, B_ci = B_ci)
      all_hotspots[[paste0("G", G, "_Thrs")]][[m]] <- fit$threshold
      all_hotspots[[paste0("G", G, "_CI")]][[m]] <- fit$CI
    }
  }

  vapply(all_hotspots, average_pairwise_jaccard, numeric(1), n = n)
}

##################################################################
# 3. Parallel Monte Carlo simulation
##################################################################
cl <- start_hotspot_cluster(num_cores, "run_likert_replication")
summary_list <- list()
raw_file_list <- character()
counter <- 0L

for (delta in delta_values) {
  for (case_id in case_values) {
    
    cat("\nRunning Likert sensitivity: delta=", delta,
        ", Case ", case_id, ", M=", M,
        ", iter_max=", iter_max, "\n", sep = "")
    
    scenario_seed <- 410000L + 100L * as.integer(delta) + case_id
    list_iter <- foreach::foreach(
      iter = seq_len(iter_max), .errorhandling = "pass",
      .options.RNG = scenario_seed) %dorng% {
      run_likert_replication(case_id, delta, n, G_values, M,alpha, epsilon, eta, B_ci)
    }
    check_replications(list_iter, iter_max,paste("delta", delta, "case", case_id))

    raw_file <- paste0("likert_raw_delta", delta, "_case", case_id, ".RData")
    save(list_iter, delta, case_id, n, M, iter_max, file = raw_file)
    raw_file_list <- c(raw_file_list, raw_file)

    for (method_name in names(method_map(G_values))) {
      
      result_name <- unname(method_map(G_values)[[method_name]])
      vals <- vapply(list_iter, function(res) unname(res[[result_name]]), numeric(1))
      rule <- if (startsWith(method_name, "Thrs")) "Thrs" else "CI"
      
      G <- as.integer(sub(".*G", "", method_name))
      counter <- counter + 1L
      summary_list[[counter]] <- data.frame(
        Delta = delta, Case = case_id, Rule = rule, G = G,
        AverageJaccard = mean(vals))
    }
    cat("Done.\n")
  }
}
parallel::stopCluster(cl)

##################################################################
# 4. Final results (ONLY iter_max replications)
##################################################################
likert_summary <- do.call(rbind, summary_list)
likert_summary$Rule <- factor(likert_summary$Rule, levels = c("Thrs", "CI"))
likert_summary <- likert_summary[order(likert_summary$Delta, likert_summary$Rule,
                                       likert_summary$G, likert_summary$Case), ]
row.names(likert_summary) <- NULL
write.csv(likert_summary, "likert_hotspot_summary.csv", row.names = FALSE)
save(likert_summary, case_par, bi_thresholds, raw_file_list,
     n, iter_max, M, alpha, epsilon, eta, B_ci, B_threshold,
     file = "likert_hotspot_results.RData")
print(likert_summary)
cat("Saved likert_hotspot_summary.csv / likert_hotspot_results.RData\n")
