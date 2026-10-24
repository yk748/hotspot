##################################################################
# Numerical Experiment 1: Univariate Change Point Detection
##################################################################

source("Library.R")

library(doParallel)
library(doRNG)
library(foreach)
library(breakfast)
library(changepoint)
library(stepR)
library(mosum)

##################################################################
# 1. Global settings
##################################################################
n <- 100
iter_max <- 500

# Tolerance used for DR(5) and FDR(5)
tol <- 5

# Window sizes for Joint-MOSUM and conventional MOSUM
G_values <- c(20, 40)

# Number of parallel workers
# num_cores <- min(50, max(1, parallel::detectCores() - 1))
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))

if (is.na(num_cores) || num_cores < 1) {
  num_cores <- 1
}

##################################################################
# 2. Change-point configurations
##################################################################
# True change-point locations
true_cp_list <- list(`1` = c(50), `2` = c(40, 60), `3` = c(25, 50, 75))

##################################################################
# 3. Six data-generating settings
##################################################################
case_par <- data.frame(case = 1:6,
                      a = c(2, 2, 2,1, 1, 1),
                      sigma_min = c(0.1, 0.1, 0.4,0.1, 0.1, 0.4),
                      sigma_max = c(0.4, 0.8, 0.8,0.4, 0.8, 0.8))

##################################################################
# 4. Data-generating functions
##################################################################
generate_piecewise_series <- function(n,
                                      true_cp,
                                      a,
                                      sigma_min,
                                      sigma_max){
  
  # Segment boundaries
  boundaries <- c(0, true_cp, n)
  
  # Number of segments
  num_segment <- length(true_cp) + 1
  
  # Random segment-specific parameters
  mu_vec <- runif(num_segment,min = -a,max = a)
  sigma_vec <- runif(num_segment,min = sigma_min,max = sigma_max)
  
  # Generate series
  x <- rep(NA_real_, n)
  for (k in seq_len(num_segment)){
    idx <- (boundaries[k] + 1):boundaries[k + 1]
    x[idx] <- rnorm(length(idx),mean = mu_vec[k],sd = sigma_vec[k])
  }
  return(x)
}

# -------------------------------------------------------------- #
# No-change scenario
#   X_t ~ N(0,1)
# -------------------------------------------------------------- #

generate_null_series <- function(n){
  rnorm(n, mean = 0, sd = 1)
}


##################################################################
# 5. Utility functions
##################################################################
# -------------------------------------------------------------- #
# Clean estimated change-point vector
# -------------------------------------------------------------- #
clean_cpts <- function(x, n){
  
  if (is.null(x) || length(x) == 0){
    return(integer(0))
  }
  x <- as.numeric(x)
  
  # Remove NA / Inf
  x <- x[is.finite(x)]
  
  # Change points must be inside the series
  x <- x[x > 0 & x < n]
  
  # Remove duplicates
  x <- sort(unique(x))
  
  return(x)
}

# -------------------------------------------------------------- #
# DR(5) and FDR(5) for one replication
#
# DR(5):
#   1 if every true change point has at least one
#   estimated change point within +/- 5; 0 otherwise.
#
# FDR(5):
#   Conditional on at least one estimated change point,
#   1 if none of the estimated change points lies within
#   +/- 5 of any true change point; 0 otherwise.
#   Replications with no detected change points are excluded
#   from the FDR denominator.
# -------------------------------------------------------------- #
evaluate_change_detection <- function(est,
                                      truth,
                                      n,
                                      tol = 5){
  
  est <- clean_cpts(est,n)
  truth <- sort(unique(truth))
  
  ##############################################################
  # DR(5): 
  ##############################################################
  detected_truth <- vapply(truth,function(t0){any(abs(est - t0) <= tol)},
                           logical(1))
  DR <- as.numeric(all(detected_truth))
  
  ##############################################################
  # FDR(5): 
  ##############################################################
  if (length(est) == 0){
    # Exclude from denominator
    FDR5 <- NA_real_
    
  } else {
    # Is at least one estimated CP within +/- tol
    # of at least one true CP?
    near_true_cp <- any(vapply(est, function(e0){any(abs(e0 - truth) <= tol)},
                               logical(1))
    )
    # 1 = every estimated CP misses every true CP neighborhood
    # 0 = at least one estimated CP is within +/- tol
    FDR5 <- as.numeric(!near_true_cp)
  }
  
  return(c(DR = DR, FDR5 = FDR5))
}


##################################################################
# 6. Run all change-point detection methods
##################################################################
run_univariate_methods <- function(x, n){
  
  ##############################################################
  # Joint-MOSUM: G = 20
  ##############################################################
  track20 <- single_tracker(x, n = n, Gl = 20, Gr = 20)
  joint20_fit <- refinements(test_stat = track20$Dist,
                              type = "both",
                              n = n,
                              Gl = 20,
                              Gr = 20,
                              alpha = 0.05,
                              epsilon = 0.3,
                              eta = 0.3)
  joint20 <- joint20_fit$cpt_eta
  
  ##############################################################
  # Joint-MOSUM: G = 40
  ##############################################################
  track40 <- single_tracker(x,n = n, Gl = 40, Gr = 40)
  joint40_fit <- refinements(test_stat = track40$Dist,
                              type = "both",
                              n = n,
                              Gl = 40,
                              Gr = 40,
                              alpha = 0.05,
                              epsilon = 0.3,
                              eta = 0.3)
  joint40 <- joint40_fit$cpt_eta
  
  ##############################################################
  # PELT
  ##############################################################
  pelt_fit <- cpt.meanvar(x,
                          penalty = "MBIC",
                          method = "PELT",
                          test.stat = "Normal",
                          class = TRUE,
                          param.estimates = TRUE)
  pelt_cp <- cpts(pelt_fit)
  
  ##############################################################
  # stepR
  ##############################################################
  stepR_fit <- stepFit(x,
                      alpha = 0.05,
                      jumpint = TRUE,
                      confband = TRUE)
  stepR_cp <- setdiff(stepR_fit$leftEnd,c(1, n))
  
  ##############################################################
  # Wild Binary Segmentation
  ##############################################################
  WBS_fit <- breakfast(x,solution.path = "wbs")
  WBS_cp <- WBS_fit$cptmodel.list[[1]]$cpts
  
  ##############################################################
  # Conventional MOSUM: G = 20
  ##############################################################
  MOSUM20_cp <- mosum(x,G = 20)$cpts.info$cpts
  
  ##############################################################
  # Conventional MOSUM: G = 40
  ##############################################################
  MOSUM40_cp <- mosum(x,G = 40)$cpts.info$cpts
  
  ##############################################################
  # Return cleaned change-point estimates
  ##############################################################
  output <- list(Joint_MOSUM_G20 = clean_cpts(joint20, n),
                Joint_MOSUM_G40 = clean_cpts(joint40, n),
                PELT = clean_cpts(pelt_cp, n),
                stepR = clean_cpts(stepR_cp, n),
                WBS = clean_cpts(WBS_cp, n),
                MOSUM_G20 = clean_cpts(MOSUM20_cp, n),
                MOSUM_G40 = clean_cpts(MOSUM40_cp, n))
  return(output)
}


##################################################################
# 7. Set up parallel backend
##################################################################
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# Load Library.R once on every worker
invisible(clusterEvalQ(cl,source("Library.R")))


##################################################################
# 8. No-change simulation
##################################################################
cat("\n")
cat("############################################################\n")
cat("Running no-change simulations\n")
cat("############################################################\n")

null_results <- foreach(
  iter = 1:iter_max,
  
  .packages = c("breakfast","changepoint","stepR","mosum"),
  .export = c("generate_null_series","clean_cpts","run_univariate_methods","n"),
  .options.RNG = 123
  
) %dorng% {
  
  # Generate N(0,1) series
  x <- generate_null_series(n)
  
  # Apply all methods
  run_univariate_methods(x = x,n = n)
}

##################################################################
# 9. Calculate FPR under no change
##################################################################
method_names <- names(null_results[[1]])

FPR_results <- sapply(method_names, function(method){
                                      mean( sapply( null_results,
                                          function(res){
                                            length(res[[method]]) > 0
                                          }
                                        )
                                      )
                                    })

FPR_table <- data.frame(Method = method_names,FPR = as.numeric(FPR_results))

cat("\n")
cat("No-change results:\n")
print(FPR_table)


# Save no-change results
save(null_results,FPR_table,file = "univariate_null_results.RData")
write.csv(FPR_table,file = "univariate_null_results.csv",row.names = FALSE)


##################################################################
# 10. Change-point simulation grid
##################################################################
scenario_grid <- expand.grid(num_change = 1:3,case = 1:6)
simulation_results <- vector(mode = "list",length = nrow(scenario_grid))

##################################################################
# 11. Run change-point simulations
##################################################################
for (s in seq_len(nrow(scenario_grid))){
  
  ##############################################################
  # Current scenario
  ##############################################################
  num_change <- scenario_grid$num_change[s]
  case_id <- scenario_grid$case[s]
  
  true_cp <- true_cp_list[[as.character(num_change)]]
  
  a <- case_par$a[case_id]
  sigma_min <- case_par$sigma_min[case_id]
  sigma_max <- case_par$sigma_max[case_id]
  
  cat("\n")
  cat("############################################################\n")
  cat("Running:",num_change,"change point(s), Case",case_id,"\n")
  cat("True change points:", true_cp, "\n")
  cat("############################################################\n")
  
  ##############################################################
  # Monte Carlo replications
  ##############################################################
  tmp_results <- foreach(iter = 1:iter_max,
    .packages = c("breakfast","changepoint","stepR","mosum"),
    .export = c("generate_piecewise_series","clean_cpts","run_univariate_methods","n",
      "true_cp","a","sigma_min","sigma_max"),
    # Different reproducible random stream for each scenario
    .options.RNG = 1000 + s
  ) %dorng% {
    ############################################################
    # Generate data
    ############################################################
    x <- generate_piecewise_series(n = n, true_cp = true_cp,
                                   a = a, sigma_min = sigma_min, sigma_max = sigma_max)
    
    ############################################################
    # Apply all methods
    ############################################################
    run_univariate_methods(x = x, n = n)
  }
  
  ##############################################################
  # Store scenario
  ##############################################################
  simulation_results[[s]] <- list(num_change = num_change,
                                  case = case_id,
                                  true_cp = true_cp,
                                  a = a,
                                  sigma_min = sigma_min,
                                  sigma_max = sigma_max,
                                  results = tmp_results)
  cat("Done.\n")
}

##################################################################
# 12. Stop parallel backend
##################################################################
stopCluster(cl)

##################################################################
# 13. Calculate DR(5) and FDR(5)
##################################################################
performance_list <- list()
cnt <- 0
for (s in seq_along(simulation_results)){
  ##############################################################
  # Current scenario
  ##############################################################
  sim <- simulation_results[[s]]
  
  true_cp <- sim$true_cp
  results <- sim$results
  
  method_names <- names(results[[1]])
  
  ##############################################################
  # Evaluate each method
  ##############################################################
  for (method in method_names){
    metrics <- t(sapply(results, function(res){ evaluate_change_detection(
                                                est = res[[method]],
                                                truth = true_cp,
                                                n = n,
                                                tol = tol
                                              )
                                            }))
    
    cnt <- cnt + 1
    performance_list[[cnt]] <- data.frame(num_change = sim$num_change,
                                          case = sim$case,
                                          method = method,
                                          
                                          # Proportion of replications recovering all true CPs within +/- 5
                                          DR = mean(metrics[, "DR"]),
                                          
                                          # Conditional on at least one detection
                                          FDR = if (all(is.na(metrics[, "FDR5"]))) {
                                            NA_real_
                                          } else {
                                            mean(metrics[, "FDR5"], na.rm = TRUE)
                                          })
  }
}

##################################################################
# 14. Combine performance results
##################################################################
performance_table <- do.call(rbind,performance_list)

# Order methods
method_order <- c("Joint_MOSUM_G20", "Joint_MOSUM_G40", "PELT", "stepR", "WBS", "MOSUM_G20", "MOSUM_G40")
performance_table$method <- factor(performance_table$method,levels = method_order)
performance_table <- performance_table[order(performance_table$num_change, performance_table$case, performance_table$method),]


# ##################################################################
# # 15. Display results
# ##################################################################
# cat("\n")
# cat("############################################################\n")
# cat("Univariate change-point simulation results\n")
# cat("############################################################\n")
# 
# print(performance_table)


##################################################################
# 16. Save results
##################################################################
save(simulation_results, performance_table,file = "univariate_change_results.RData")
write.csv(performance_table,file = "univariate_change_results.csv",row.names = FALSE)


# ##################################################################
# # 17. Display results by scenario
# ##################################################################
# for (num_change in 1:3){
#   for (case_id in 1:6){
#     cat("\n")
#     cat("# Change points:",num_change,"| Case:",case_id,"\n")
#     tmp <- performance_table[performance_table$num_change == num_change & performance_table$case == case_id,]  
#     print(tmp)
#   }
# }

