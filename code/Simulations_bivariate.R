##################################################################
# Numerical Experiment 1: Bivariate Change Point Detection
##################################################################
source("Library.R")

library(doParallel)
library(doRNG)
library(foreach)

library(mvtnorm)
library(stringr)
library(ecp)
library(hdbinseg)
library(changepoints)

##################################################################
# 1. Global settings
##################################################################
n <- 100
iter_max <- 500

# Tolerance for DR
tol <- 5

# Number of parallel workers
# num_cores <- min(50, max(1, parallel::detectCores() - 1))
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))

if (is.na(num_cores) || num_cores < 1) {
  num_cores <- 1
}

##################################################################
# 2. Bivariate change-point configurations
##################################################################
# Scenario 1:
#   Series 1: change at 60
#   Series 2: change at 40
#
# Scenario 2:
#   Series 1: change at 50
#   Series 2: changes at 40, 60
#
# Scenario 3:
#   Series 1: changes at 40, 60
#   Series 2: changes at 30, 70
scenario_cp <- list(`1` = list(series1 = c(60),series2 = c(40)),
                    `2` = list(series1 = c(50),series2 = c(40, 60)),
                    `3` = list(series1 = c(40, 60),series2 = c(30, 70)))

##################################################################
# 3. Six data-generating settings
##################################################################
case_par <- data.frame(case = 1:6,
                       a = c(2, 2, 2, 1, 1, 1),
                       sigma_min = c(0.1, 0.1, 0.4, 0.1, 0.1, 0.4),  
                       sigma_max = c(0.4, 0.8, 0.8, 0.4, 0.8, 0.8))

##################################################################
# 4. Data-generating function
##################################################################
generate_piecewise_series <- function(n,
                                      true_cp,
                                      a,
                                      sigma_min,
                                      sigma_max){
  
  boundaries <- c(0,true_cp,n)
  num_segment <- length(true_cp) + 1
  
  # Segment-specific means
  mu_vec <- runif(num_segment,min = -a,max = a)
  
  # Segment-specific standard deviations
  sigma_vec <- runif(num_segment,min = sigma_min,max = sigma_max)

  x <- rep(NA_real_,n)
  for (k in seq_len(num_segment)) {
    idx <- (boundaries[k] + 1):boundaries[k + 1]
    x[idx] <- rnorm(length(idx),mean = mu_vec[k],sd = sigma_vec[k])
  }
  
  return(x)
}


##################################################################
# 5. Utility functions
##################################################################
# -------------------------------------------------------------- #
# Clean estimated change-point vector
# -------------------------------------------------------------- #
clean_cpts <- function(x, n){
  
  if (is.null(x) || length(x) == 0) {
    return(integer(0))
  }
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  x <- x[x > 0 & x < n]
  x <- sort(unique(x))
  
  return(x)
}

# -------------------------------------------------------------- #
# DR(5) and FDR(0) for one replication
# DR(5): 1 if every pooled true change point has at least one
#   estimated change point within +/- 5; 0 otherwise.
# FDR(0): Conditional on at least one estimated change point,
#   1 if none of the estimated change points exactly matches
#   any pooled true change point; 0 otherwise.
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
  detected_truth <- vapply(truth, function(t0) {any(abs(est - t0) <= tol)},
                           logical(1))
  DR <- as.numeric(all(detected_truth))
  
  ##############################################################
  # FDR(0): 
  ##############################################################
  if (length(est) == 0) {
    # Exclude from denominator
    FDR0 <- NA_real_
    
  } else {
    # Is there at least one exact match to any true CP?
    exact_match <- any(est %in% truth)
    
    # 1 = all detections miss every true CP exactly
    # 0 = at least one exact match
    FDR0 <- as.numeric(!exact_match)
  }
  return(c(DR = DR,FDR0 = FDR0))
}

##################################################################
# 6. Run all bivariate change-point detection methods
##################################################################
run_bivariate_methods <- function(series1,
                                  series2,
                                  n){
  ##############################################################
  # Bi-MOSUM: G = 20
  ##############################################################
  bi20_fit <- bivariate_mosum(series1,series2,
                              n = n, Gl = 20, Gr = 20,
                              alpha = 0.05, eta = 0.3)
  
  bi20_cp <- sort(unique(c(bi20_fit$YX$cpt_eta,bi20_fit$YX2$cpt_eta,
                           bi20_fit$Y2X$cpt_eta,bi20_fit$Y2X2$cpt_eta)))
  
  ##############################################################
  # Bi-MOSUM: G = 40
  ##############################################################
  bi40_fit <- bivariate_mosum(series1,series2,
                              n = n, Gl = 40, Gr = 40,
                              alpha = 0.05, eta = 0.3)
  
  bi40_cp <- sort(unique(c(bi40_fit$YX$cpt_eta,
                           bi40_fit$YX2$cpt_eta,
                           bi40_fit$Y2X$cpt_eta,
                           bi40_fit$Y2X2$cpt_eta)))
  
  ##############################################################
  # ECP
  ##############################################################
  ecp_alpha1_raw <- e.divisive(cbind(series1,series2),sig.lvl = 0.05,
                               R = 100,alpha = 1,min.size = 10)$estimates
  ecp_alpha2_raw <- e.divisive(cbind(series1,series2),sig.lvl = 0.05,
                               R = 100,alpha = 2,min.size = 10)$estimates
  
  # e.divisive returns starts of segments:
  # e.g., a boundary after observation 40 is returned as 41.
  ecp_alpha1 <- setdiff(ecp_alpha1_raw,c(1,n + 1)) - 1
  ecp_alpha2 <- setdiff(ecp_alpha2_raw,c(1,n + 1)) - 1
  
  ECP_cp <- sort(unique(c(ecp_alpha1,ecp_alpha2)))
  
  
  ##############################################################
  # SBS
  ##############################################################
  SBS_type1 <- sbs.alg(rbind(series1,series2),cp.type = 1,do.parallel = 0)$ecp
  SBS_type2 <- sbs.alg(rbind(series1,series2),cp.type = 2,do.parallel = 0)$ecp
  SBS_cp <- sort(unique(c(SBS_type1,SBS_type2)))
  
  ##############################################################
  # DCBS
  ##############################################################
  DCBS_type1 <- dcbs.alg(rbind(series1,series2),cp.type = 1,do.parallel = 0)$ecp
  DCBS_type2 <- dcbs.alg(rbind(series1,series2),cp.type = 2,do.parallel = 0)$ecp
  DCBS_cp <- sort(unique(c(DCBS_type1,DCBS_type2)))

  ##############################################################
  # Nonparametric multivariate method
  ##############################################################
  intervals <- WBS.intervals(M = 10,lower = 1,upper = n)
  temp <- WBS.multi.nonpar(Y = rbind(series1,series2),
                           W = rbind(series1,series2),1,n,
                           intervals$Alpha,
                           intervals$Beta,
                           h = 5 * (30 * log(n) / n)^(1 / 2),
                           delta = 10)
  tau_NMP <- median(temp$Dval)
  
  if (length(unique(temp$Level)) <= 1) {
    NMP_cp <- temp$S[is.finite(temp$Dval) & temp$Dval > tau_NMP]
    NMP_cp <- sort(unique(NMP_cp))
    
  } else {
    NMP_fit <- thresholdBS(temp,tau_NMP)
    if (is.null(NMP_fit$cpt_hat)) {
      NMP_cp <- integer(0)
    } else {
      NMP_cp <- sort(unique(NMP_fit$cpt_hat$location))
    }
  }
  
  ##############################################################
  # Return cleaned change-point estimates
  ##############################################################
  output <- list(Bi_MOSUM_G20 = clean_cpts(bi20_cp,n),
                 Bi_MOSUM_G40 = clean_cpts(bi40_cp,n),
                 ECP = clean_cpts(ECP_cp,n),
                 SBS = clean_cpts(SBS_cp,n),
                 DCBS = clean_cpts(DCBS_cp,n),
                 NMP = clean_cpts(NMP_cp,n))
  return(output)
}

##################################################################
# 7. Simulation grid
##################################################################
scenario_grid <- expand.grid(scenario = 1:3,
                             case = 1:6)

simulation_results <- vector(mode = "list",length = nrow(scenario_grid))

##################################################################
# 8. Set up parallel backend
##################################################################
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# Load Library.R once on every worker
invisible(clusterEvalQ(cl,source("Library.R")))

##################################################################
# 9. Run bivariate simulations
##################################################################
for (s in seq_len(nrow(scenario_grid))) {
  
  ##############################################################
  # Current scenario
  ##############################################################
  scenario_id <- scenario_grid$scenario[s]
  case_id <- scenario_grid$case[s]
  
  true_cp1 <- scenario_cp[[as.character(scenario_id)]]$series1
  true_cp2 <- scenario_cp[[as.character(scenario_id)]]$series2
  
  # Pooled distinct true change points
  pooled_true_cp <- sort(unique(c(true_cp1,true_cp2)))
  
  # Current data-generating parameters
  a <- case_par$a[case_id]
  sigma_min <- case_par$sigma_min[case_id]
  sigma_max <- case_par$sigma_max[case_id]
  
  ##############################################################
  # Minimal progress information for Slurm log
  ##############################################################
  cat("\n")
  cat("Scenario:",scenario_id,"| Case:",case_id,"| Series 1 CP:",
      paste(true_cp1, collapse = ","),"| Series 2 CP:",
      paste(true_cp2, collapse = ","),"| Pooled CP:",
      paste(pooled_true_cp, collapse = ","),"\n")
  
  ##############################################################
  # Monte Carlo replications
  ##############################################################
  tmp_results <- foreach(iter = 1:iter_max,
                          .packages = c("mvtnorm",
                                        "stringr",
                                        "ecp",
                                        "hdbinseg",
                                        "changepoints"),
                         .export = c("generate_piecewise_series",
                                    "clean_cpts",
                                    "run_bivariate_methods",
                                    "n",
                                    "true_cp1",
                                    "true_cp2",
                                    "a",
                                    "sigma_min",
                                    "sigma_max"),
                         .errorhandling = "stop",
                         .options.RNG = 2000 + s) %dorng% {
    
    ############################################################
    # Generate Series 1
    ############################################################
    series1 <- generate_piecewise_series(n = n,true_cp = true_cp1,
                                         a = a,sigma_min = sigma_min,
                                         sigma_max = sigma_max)
    
    ############################################################
    # Generate Series 2 independently
    ############################################################
    series2 <- generate_piecewise_series(n = n,true_cp = true_cp2,
                                         a = a,sigma_min = sigma_min,
                                         sigma_max = sigma_max)
    
    ############################################################
    # Apply all methods
    ############################################################
    run_bivariate_methods(series1 = series1,series2 = series2,n = n)
  }
  
  ##############################################################
  # Store scenario
  ##############################################################
  simulation_results[[s]] <- list(scenario = scenario_id,
                                  case = case_id,
                                  true_cp1 = true_cp1,
                                  true_cp2 = true_cp2,
                                  pooled_true_cp = pooled_true_cp,
                                  a = a,
                                  sigma_min = sigma_min,
                                  sigma_max = sigma_max,
                                  results = tmp_results)
  cat("Done.\n")
}

##################################################################
# 10. Stop parallel backend
##################################################################
stopCluster(cl)

##################################################################
# 11. Calculate DR(5) and FDR(0)
##################################################################
performance_list <- list()
cnt <- 0

for (s in seq_along(simulation_results)) {
  
  sim <- simulation_results[[s]]
  pooled_true_cp <- sim$pooled_true_cp
  results <- sim$results
  
  method_names <- names(results[[1]])
  for (method in method_names) {
    metrics <- t(sapply(results, function(res) {
          evaluate_change_detection(est = res[[method]],
                                    truth = pooled_true_cp,
                                    n = n,tol = tol)
        }
      )
    )
    
    cnt <- cnt + 1
    
    performance_list[[cnt]] <- data.frame(scenario = sim$scenario,
                                          case = sim$case,
                                          n_cp_series1 = length(sim$true_cp1),
                                          n_cp_series2 = length(sim$true_cp2),
                                          pooled_K = length(pooled_true_cp),
                                          method = method,
                                          
                                          # Proportion of replications recovering all true CPs within +/- 5
                                          DR = mean(metrics[, "DR"]),
                                          
                                          # Conditional on at least one detection
                                          FDR = if (all(is.na(metrics[, "FDR0"]))) {
                                            NA_real_
                                          } else {
                                            mean(metrics[, "FDR0"], na.rm = TRUE)
                                          })
  }
}

##################################################################
# 12. Combine performance results
##################################################################
performance_table <- do.call(rbind,performance_list)

##################################################################
# 13. Order methods
##################################################################
method_order <- c("Bi_MOSUM_G20","Bi_MOSUM_G40","ECP","SBS","DCBS","NMP")
performance_table$method <- factor(performance_table$method,levels = method_order)
performance_table <- performance_table[order(performance_table$scenario,
                                             performance_table$case,
                                             performance_table$method),]

##################################################################
# 14. Save results
##################################################################
save(simulation_results,performance_table,file = "bivariate_change_results.RData")
write.csv(performance_table,file = "bivariate_change_results.csv",row.names = FALSE)

