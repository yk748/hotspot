##################################################################
# Numerical Experiment 2: Hotspot Simulation
##################################################################
source("Library.R")

library(doParallel)
library(doRNG)
library(foreach)

##################################################################
# 1. Global settings
##################################################################
n <- 100
iter_max <- 200

alpha <- 0.05
epsilon <- 0.2
eta <- 0.3

# Number of bootstrap replications for pointwise CIs
B_ci <- 1000

# Number of Monte Carlo replications used to calibrate the
# Bi-MOSUM threshold for each window size
B_threshold <- 1000

G_values <- c(20, 40)
delta_values <- c(0, 5)

# Number of parallel workers supplied by Slurm
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))

if (is.na(num_cores) || num_cores < 1) {
  num_cores <- 1
}

##################################################################
# 2. Six data-generating settings
#    Same settings as Section 5.1
##################################################################
case_par <- data.frame(case = 1:6,
                      a = c(2, 2, 2, 1, 1, 1),
                      sigma_min = c(0.1, 0.1, 0.4, 0.1, 0.1, 0.4),
                      sigma_max = c(0.4, 0.8, 0.8, 0.4, 0.8, 0.8))

##################################################################
# 3. Data-generating function
##################################################################
generate_piecewise_series <- function(n,
                                      true_cp,
                                      a,
                                      sigma_min,
                                      sigma_max) {
  
  true_cp <- sort(unique(as.integer(true_cp)))
  boundaries <- c(0, true_cp, n)
  num_segment <- length(true_cp) + 1
  
  # Segment-specific means and standard deviations are generated
  # independently for each series and each replication.
  mu_vec <- runif(num_segment, min = -a, max = a)
  sigma_vec <- runif(num_segment,min = sigma_min,max = sigma_max)
  x <- rep(NA_real_, n)
  
  for (k in seq_len(num_segment)) {
    idx <- (boundaries[k] + 1):boundaries[k + 1]
    x[idx] <- rnorm(length(idx),mean = mu_vec[k],sd = sigma_vec[k])
  }
  x
}

##################################################################
# 4. Change-point configurations
##################################################################
get_true_change_points <- function(num_change,
                                   delta) {
  
  if (num_change == 1) {
    true_X <- 50
    true_Y <- 50 + delta
    
  } else if (num_change == 2) {
    true_X <- c(40, 60)
    true_Y <- c(40 + delta, 60 + delta)
    
  } else {
    stop("num_change must be 1 or 2.")
  }
  
  list(X = true_X,Y = true_Y)
}

##################################################################
# 5. Interval utilities
##################################################################
normalize_intervals <- function(intervals,
                                n) {
  
  if (is.null(intervals) || length(intervals) == 0) {
    return(list())
  }
  
  # Allow a single numeric interval c(left, right)
  if (is.numeric(intervals) && length(intervals) == 2) {
    intervals <- list(intervals)
  }
  
  out <- list()
  for (z in intervals) {
    if (is.null(z) || length(z) < 2 ||
        any(!is.finite(z[1:2]))) {
      next
    }
    
    left <- max(1L, as.integer(ceiling(min(z[1:2]))))
    right <- min(n, as.integer(floor(max(z[1:2]))))
    if (left <= right) {
      out[[length(out) + 1L]] <- c(left, right)
    }
  }
  if (length(out) == 0) {
    return(list())
  }
  
  # Merge overlapping or adjacent intervals
  ord <- order(vapply(out, `[`, numeric(1), 1))
  out <- out[ord]
  
  merged <- list(out[[1]])
  if (length(out) >= 2) {
    for (j in 2:length(out)) {
      last <- merged[[length(merged)]]
      cur <- out[[j]]
      
      if (cur[1] <= last[2] + 1L) {
        merged[[length(merged)]] <- c(last[1], max(last[2], cur[2]))
      } else {
        merged[[length(merged) + 1L]] <- cur
      }
    }
  }
  merged
}

points_to_intervals <- function(x, n) {
  if (is.null(x) || length(x) == 0) {
    return(list())
  }
  
  x <- sort(unique(as.integer(x)))
  x <- x[is.finite(x) & x >= 1 & x <= n]
  
  if (length(x) == 0) {
    return(list())
  }
  
  breaks <- c(0L, which(diff(x) > 1L), length(x))
  intervals <- lapply(seq_len(length(breaks) - 1L), function(j) {c(x[breaks[j] + 1L],x[breaks[j + 1L]])})
  normalize_intervals(intervals, n)
}

union_interval_sets <- function(..., n) {
  sets <- list(...)
  all_intervals <- list()
  
  for (s in sets) {
    s <- normalize_intervals(s, n)
    if (length(s) > 0) {
      all_intervals <- c(all_intervals, s)
    }
  }
  normalize_intervals(all_intervals, n)
}

intersect_interval_sets <- function(A, B, n) {
  
  A <- normalize_intervals(A, n)
  B <- normalize_intervals(B, n)
  
  if (length(A) == 0 || length(B) == 0) {
    return(list())
  }
  
  intersections <- list()
  for (a in A) {
    for (b in B) {
      left <- max(a[1], b[1])
      right <- min(a[2], b[2])
      
      if (left <= right) {
        intersections[[length(intersections) + 1L]] <- c(left, right)
      }
    }
  }
  normalize_intervals(intersections, n)
}

##################################################################
# 6. Fixed Bi-MOSUM threshold calibration
#    The Monte Carlo threshold is calibrated at the same G used
#    by the corresponding MOSUM procedure. 
##################################################################
max_bi_mosum_null <- function(n, G) {
  WE <- cumsum(rnorm(n))
  WV <- cumsum(rnorm(n))
  
  LE <- (WE[(2 * G + 1):n] - 2 * WE[(G + 1):(n - G)] + WE[1:(n - 2 * G)]) / sqrt(2 * G)
  LV <- (WV[(2 * G + 1):n] - 2 * WV[(G + 1):(n - G)] + WV[1:(n - 2 * G)]) / sqrt(2 * G)
  max(sqrt(LE^2 + LV^2))
}

set.seed(20261012)
bi_thresholds <- vapply(G_values, function(G) {unname(quantile(
        replicate(B_threshold,max_bi_mosum_null(n = n, G = G)),
        probs = 1 - alpha,
        type = 1))
  },
  numeric(1)
)

names(bi_thresholds) <- as.character(G_values)
cat("Bi-MOSUM Monte Carlo thresholds:\n")
print(bi_thresholds)

##################################################################
# 7. Refinement function with G-specific Bi-MOSUM threshold
#    This local definition overrides refinements() from Library.R
#    after the library is sourced on each worker.
##################################################################
refinements_fixed <- function(test_stat,
                              type,
                              n,
                              Gl,
                              Gr,
                              alpha,
                              epsilon,
                              eta) {
  
  if (Gl != Gr) {
    stop("This simulation requires Gl = Gr.")
  }
  
  G <- Gl
  if (type == "univ") {
    an <- sqrt(2 * log(n / G))
    bn <- 2 * log(n / G) + 0.5 * log(log(n / G)) + log(3 / 2) - 0.5 * log(pi)
    c_alpha <- -log(log(1 / sqrt(1 - alpha)))
    Dn <- (bn + c_alpha) / an
    
  } else if (type == "bi") {
    key <- as.character(G)
    
    if (!key %in% names(bi_thresholds)) {
      stop("No precomputed Bi-MOSUM threshold for G = ", G)
    }
    Dn <- bi_thresholds[[key]]
    
  } else {
    stop("type must be 'univ' or 'bi'.")
  }
  
  idx_cand <- which(is.finite(test_stat) & abs(test_stat) >= Dn)
  
  # Epsilon criterion
  if (length(idx_cand) > 0) {
    grp <- cumsum(c(TRUE, diff(idx_cand) != 1))
    runs <- split(idx_cand, grp)
    cpt_eps <- integer(0)
    
    for (idx in runs) {
      left <- min(idx) - 1
      right <- max(idx) + 1
      
      if ((right - left) >= epsilon * G) {
        score <- abs(test_stat[idx])
        score[!is.finite(score)] <- -Inf
        
        if (any(is.finite(score))) { cpt_eps <- c(cpt_eps, idx[which.max(score)])
        }
      }
    }
    cpt_eps <- sort(unique(cpt_eps))
    if (length(cpt_eps) == 0) {
      cpt_eps <- NULL
    }
    
  } else {
    cpt_eps <- NULL
  }
  
  # Eta criterion
  if (length(idx_cand) > 0) {
    band_l <- floor(eta * Gl)
    band_r <- floor(eta * Gr)
    cpt_eta <- integer(0)
    
    for (idx in idx_cand) {
      idx_band <- max(1, idx - band_l): min(n, idx + band_r)
      score <- abs(test_stat[idx_band])
      score[!is.finite(score)] <- -Inf
      center_score <- abs(test_stat[idx])
      
      if (is.finite(center_score) &&
          center_score >= max(score)) {
        cpt_eta <- c(cpt_eta, idx)
      }
    }
    cpt_eta <- sort(unique(cpt_eta))
    if (length(cpt_eta) == 0) {
      cpt_eta <- NULL
    }
    
  } else {
    cpt_eta <- NULL
  }
  
  list(Dn = Dn, idx_cand = idx_cand, cpt_eta = cpt_eta, cpt_eps = cpt_eps)
}

##################################################################
# 8. Pointwise confidence intervals
#    Only pointwise CIs are needed for the hotspot rule in Section
#    5.2, so these functions avoid computing unused uniform CIs.
#    The radius is the empirical (1-alpha) quantile of the absolute
#    bootstrap localization errors, matching the manuscript
#    definition. quantile(type = 1) keeps interval endpoints integer.
##################################################################
clean_cpts <- function(x, n) {
  if (is.null(x) || length(x) == 0) {
    return(integer(0))
  }
  
  x <- as.integer(round(x))
  x <- x[is.finite(x)]
  x <- x[x > 0 & x < n]
  sort(unique(x))
}

local_argmax <- function(stat, center, G, n) {
  idx <- max(1L, center - G + 1L): min(n, center + G)
  score <- abs(stat[idx])
  score[!is.finite(score)] <- -Inf
  if (all(score == -Inf)) {
    return(as.integer(center))
  }
  as.integer(idx[which.max(score)])
}

joint_mosum_pointwise_ci <- function(series, refinement, n, G, B = 1000, alpha = 0.05) {
  
  cps <- clean_cpts(refinement$cpt_eta, n)
  if (length(cps) == 0) {
    return(list())
  }
  
  k_index <- c(0L, cps, n)
  boot_cps <- matrix(NA_integer_, nrow = B, ncol = length(cps))
  for (b in seq_len(B)) {
    samples <- integer(0)
    for (k in seq_len(length(k_index) - 1L)) {
      segment_idx <- (k_index[k] + 1L): k_index[k + 1L]
      samples <- c(samples,sample(segment_idx,length(segment_idx),replace = TRUE))
    }
    series_b <- series[samples]
    track_b <- single_tracker(series_b, n = n, Gl = G, Gr = G)
    for (j in seq_along(cps)) {
      boot_cps[b, j] <- local_argmax(track_b$Dist, center = cps[j], G = G, n = n)
    }
  }
  
  out <- vector("list", length(cps))
  for (j in seq_along(cps)) {
    deviation <- abs(boot_cps[, j] - cps[j])
    radius <- as.integer(quantile(deviation, probs = 1 - alpha, type = 1, na.rm = TRUE))
    out[[j]] <- c(max(1L, cps[j] - radius), min(n, cps[j] + radius))
  }
  normalize_intervals(out, n)
}

bivariate_mosum_pointwise_ci <- function(series_Y, series_X, bi_refinement,
                                         moment, n, G, B = 1000, alpha = 0.05) {
  
  moment <- match.arg(moment,c("mean_mean","mean_var", "var_mean", "var_var"))
  cps <- clean_cpts(bi_refinement$cpt_eta, n)
  
  if (length(cps) == 0) {
    return(list())
  }
  
  k_index <- c(0L, cps, n)
  boot_cps <- matrix(NA_integer_, nrow = B, ncol = length(cps))
  for (b in seq_len(B)) {
    
    samples <- integer(0)
    for (k in seq_len(length(k_index) - 1L)) {
      segment_idx <- (k_index[k] + 1L): k_index[k + 1L]
      samples <- c(samples, sample(segment_idx, length(segment_idx), replace = TRUE))
    }
    
    # Use the same resampled time indices for Y and X to preserve
    # their within-segment paired structure.
    Y_b <- series_Y[samples]
    X_b <- series_X[samples]
    
    Y_track_b <- single_tracker(Y_b, n = n, Gl = G, Gr = G)
    X_track_b <- single_tracker(X_b, n = n, Gl = G, Gr = G)
    bi_track_b <- bivariate_tracker(Y_track_b, X_track_b, n = n, Gl = G, Gr = G)
    stat_b <- switch(moment, mean_mean = bi_track_b$Dist_YX,
                     mean_var  = bi_track_b$Dist_YX2,
                     var_mean  = bi_track_b$Dist_Y2X,
                     var_var   = bi_track_b$Dist_Y2X2)
    
    for (j in seq_along(cps)) {
      boot_cps[b, j] <- local_argmax(stat_b, center = cps[j], G = G, n = n)
    }
  }
  
  out <- vector("list", length(cps))
  for (j in seq_along(cps)) {
    deviation <- abs(boot_cps[, j] - cps[j])
    
    radius <- as.integer(quantile(deviation, probs = 1 - alpha, type = 1, na.rm = TRUE))
    out[[j]] <- c(max(1L, cps[j] - radius), min(n, cps[j] + radius))
  }
  normalize_intervals(out, n)
}

##################################################################
# 9. Construct thresholding and CI hotspots for one window size
##################################################################
construct_hotspots <- function(series_X,
                               series_Y,
                               G,
                               n,
                               alpha,
                               epsilon,
                               eta,
                               B_ci) {
  
  fit <- bivariate_mosum(series_X, series_Y, n = n, Gl = G, Gr = G,
                         alpha = alpha, epsilon = epsilon, eta = eta)

  # Thresholding rule
  cross_threshold_points <- sort(unique(c(fit$YX$idx_cand,
                                          fit$YX2$idx_cand,
                                          fit$Y2X$idx_cand,
                                          fit$Y2X2$idx_cand)))
  Y_threshold_points <- sort(unique(fit$Y_both$idx_cand))
  threshold_points <- intersect(cross_threshold_points, Y_threshold_points)
  hotspot_threshold <- points_to_intervals(threshold_points, n = n)
  
  # CI rule: Y intervals
  CI_Y <- joint_mosum_pointwise_ci(series = series_Y, refinement = fit$Y_both, n = n,
                                   G = G, B = B_ci, alpha = alpha)
  
  # CI rule: four cross-series detector combinations
  CI_YX <- bivariate_mosum_pointwise_ci(series_Y = series_Y, series_X = series_X,
                                        bi_refinement = fit$YX, moment = "mean_mean",
                                        n = n,G = G, B = B_ci, alpha = alpha)
  CI_YX2 <- bivariate_mosum_pointwise_ci(series_Y = series_Y, series_X = series_X,
                                         bi_refinement = fit$YX2, moment = "mean_var",
                                         n = n, G = G, B = B_ci, alpha = alpha)
  CI_Y2X <- bivariate_mosum_pointwise_ci(series_Y = series_Y, series_X = series_X, 
                                         bi_refinement = fit$Y2X, moment = "var_mean",
                                         n = n, G = G, B = B_ci, alpha = alpha)
  CI_Y2X2 <- bivariate_mosum_pointwise_ci(series_Y = series_Y, series_X = series_X,
                                          bi_refinement = fit$Y2X2, moment = "var_var",
                                          n = n, G = G, B = B_ci, alpha = alpha)
  
  # First form the UNION of all four cross-series CI sets,
  # then intersect that union with the Y CI set.
  cross_CI_union <- union_interval_sets(CI_YX, CI_YX2, CI_Y2X, CI_Y2X2, n = n)
  hotspot_CI <- intersect_interval_sets(cross_CI_union, CI_Y, n = n)
  
  list(threshold = hotspot_threshold, CI = hotspot_CI)
}

##################################################################
# 10. Evaluate hotspot hit rate and hotspot interval length
##################################################################
evaluate_hotspot <- function(hotspots, true_Y, n) {
  
  hotspots <- normalize_intervals(hotspots, n)
  true_Y <- sort(unique(as.integer(true_Y)))
  K_true <- length(true_Y)
  
  hit <- logical(K_true)
  cp_length <- rep(n, K_true)
  if (length(hotspots) > 0) {
    for (j in seq_along(true_Y)) {
      
      containing <- which(vapply(hotspots, function(z) {true_Y[j] >= z[1] && true_Y[j] <= z[2]},
                                 logical(1)))
      if (length(containing) > 0) {
        hit[j] <- TRUE
        
        # Intervals have already been merged, so a true change point
        # should belong to at most one connected interval.
        z <- hotspots[[containing[1]]]
        cp_length[j] <- z[2] - z[1] + 1L
      }
    }
  }
  
  c(HR = mean(hit), HotspotLength = mean(cp_length))
}

##################################################################
# 11. Run one Monte Carlo replication
##################################################################
run_one_replication <- function(num_change,
                                case_id,
                                delta,
                                n,
                                alpha,
                                epsilon,
                                eta,
                                B_ci) {
  
  pars <- case_par[case_par$case == case_id, ]
  truth <- get_true_change_points(num_change = num_change,delta = delta)
  
  # X and Y are generated independently, with independent
  # segment-specific parameters.
  series_X <- generate_piecewise_series(n = n, true_cp = truth$X, a = pars$a, 
                                        sigma_min = pars$sigma_min, 
                                        sigma_max = pars$sigma_max)
  series_Y <- generate_piecewise_series(n = n, true_cp = truth$Y, a = pars$a,
                                        sigma_min = pars$sigma_min,
                                        sigma_max = pars$sigma_max)
  fit_G20 <- construct_hotspots(series_X = series_X, series_Y = series_Y,
                                G = 20, n = n, alpha = alpha,
                                epsilon = epsilon,
                                eta = eta, B_ci = B_ci)
  fit_G40 <- construct_hotspots(series_X = series_X, series_Y = series_Y,
                                G = 40, n = n, alpha = alpha,
                                epsilon = epsilon,
                                eta = eta, B_ci = B_ci)
  
  list(G20_Thrs = fit_G20$threshold,
       G40_Thrs = fit_G40$threshold,
       G20_CI = fit_G20$CI,
       G40_CI = fit_G40$CI)
}

##################################################################
# 12. Set up parallel backend
##################################################################
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Load Library.R once on each worker, then replace the old
# refinements() definition with the G-specific version above.
clusterExport(cl, varlist = c("bi_thresholds", "refinements_fixed",
                              "case_par", "generate_piecewise_series",
                              "get_true_change_points",
                              "normalize_intervals",
                              "points_to_intervals",
                              "union_interval_sets",
                              "intersect_interval_sets",
                              "clean_cpts", "local_argmax",
                              "joint_mosum_pointwise_ci",
                              "bivariate_mosum_pointwise_ci",
                              "construct_hotspots", "run_one_replication"),
              envir = environment())

clusterEvalQ(cl, {
    source("Library.R")
    refinements <- refinements_fixed
    NULL
  })

##################################################################
# 13. Run all hotspot simulation settings
##################################################################
summary_list <- list()
summary_count <- 0L
raw_file_list <- character(0)

for (num_change in c(1, 2)) {
  
  for (delta in delta_values) {
    
    for (case_id in 1:6) {
      
      cat("\n############################################################\n",
          "Running hotspot simulation: ",
          num_change, " change point(s), Case ", case_id, ", delta = ", delta, "\n",
          "############################################################\n", sep = "")
      
      # Scenario-specific seed gives reproducible parallel streams.
      scenario_seed <- 100000L + 10000L * num_change + 100L * case_id + as.integer(delta)
      
      list_iter <- foreach(
        iter = seq_len(iter_max),
        .errorhandling = "pass",
        .options.RNG = scenario_seed
      ) %dorng% {
        
        run_one_replication(num_change = num_change, case_id = case_id, 
                            delta = delta, n = n, alpha = alpha, 
                            epsilon = epsilon,
                            eta = eta,
                            B_ci = B_ci)
      }
      
      # Do not silently drop failed replications
      failed <- vapply(list_iter, function(x) inherits(x, "error"),
                       logical(1))
      if (any(failed)) {
        first_error <- list_iter[[which(failed)[1]]]
        
        stop(sum(failed), " of ", iter_max,
             " replications failed for num_change = ",
             num_change, ", case = ", case_id, ", delta = ",
             delta, ". First error: ",
          conditionMessage(first_error))
      }
      
      if (length(list_iter) != iter_max) {
        stop("Expected ", iter_max, " replications but received ", length(list_iter), ".")
      }
      
      # Save raw hotspot results for this setting
      truth <- get_true_change_points(num_change = num_change, delta = delta)
      
      raw_file <- paste0("hotspot_raw_change", num_change, "_case",
                         case_id, "_delta", delta, ".RData")
      save(list_iter, truth, num_change, case_id, delta, file = raw_file)
      raw_file_list <- c(raw_file_list, raw_file)
      
      # Evaluate the four rule/window combinations
      method_map <- list(Thrs_G20 = "G20_Thrs", Thrs_G40 = "G40_Thrs",
                         CI_G20 = "G20_CI", CI_G40 = "G40_CI")
      
      for (method_name in names(method_map)) {
        object_name <- method_map[[method_name]]
        
        metrics <- t(vapply(list_iter, function(res) {
              evaluate_hotspot(hotspots = res[[object_name]],true_Y = truth$Y, n = n)
            }, numeric(2)))
        
        rule <- if (grepl("^Thrs", method_name)) {
          "Thrs"
        } else {
          "CI"
        }
        
        G <- if (grepl("G20$", method_name)) {
          20
        } else {
          40
        }
        
        summary_count <- summary_count + 1L
        
        summary_list[[summary_count]] <- data.frame(num_change = num_change,
                                                    case = case_id,
                                                    delta = delta,
                                                    rule = rule, G = G,
                                                    HitRate = mean(metrics[, "HR"]),
                                                    HotspotLength = mean(metrics[, "HotspotLength"]))
      }
      cat("Done.\n")
    }
  }
}

##################################################################
# 14. Stop parallel backend
##################################################################
stopCluster(cl)

##################################################################
# 15. Combine and order results
##################################################################
hotspot_summary <- do.call(rbind, summary_list)
hotspot_summary$rule <- factor(hotspot_summary$rule, levels = c("Thrs", "CI"))
hotspot_summary <- hotspot_summary[order(hotspot_summary$num_change,
                                         hotspot_summary$delta,
                                         hotspot_summary$rule,
                                         hotspot_summary$G,
                                         hotspot_summary$case),]
row.names(hotspot_summary) <- NULL

##################################################################
# 16. Save final results
##################################################################
write.csv(hotspot_summary, file = "hotspot_simulation_summary.csv", row.names = FALSE)
save(hotspot_summary, case_par, bi_thresholds, raw_file_list,
     n, iter_max, alpha, epsilon, eta, B_ci, B_threshold,
     file = "hotspot_simulation_results.RData")

##################################################################
# 17. Display results
##################################################################
print(hotspot_summary)

cat("\nSaved:\n")
cat("  hotspot_simulation_summary.csv\n")
cat("  hotspot_simulation_results.RData\n")
cat("  hotspot_raw_change*_case*_delta*.RData\n")
