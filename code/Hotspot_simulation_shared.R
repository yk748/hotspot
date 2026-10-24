##################################################################
# Shared utilities for supplementary hotspot simulations
# Analytical functions copied from Simulations_hotspot.R without
# changing the main-simulation hotspot definitions or CI procedure.
# Run from a directory containing Library.R and this file.
##################################################################
source("Library.R")
library(doParallel)
library(doRNG)
library(foreach)

##################################################################
# Original DGP and interval utilities (Sections 2--5 of main script)
##################################################################

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

##################################################################
# Monte Carlo calibration: same null statistic and seed as main script
##################################################################
max_bi_mosum_null <- function(n, G) {
  WE <- cumsum(rnorm(n))
  WV <- cumsum(rnorm(n))
  LE <- (WE[(2 * G + 1):n] - 2 * WE[(G + 1):(n - G)] +
         WE[1:(n - 2 * G)]) / sqrt(2 * G)
  LV <- (WV[(2 * G + 1):n] - 2 * WV[(G + 1):(n - G)] +
         WV[1:(n - 2 * G)]) / sqrt(2 * G)
  max(sqrt(LE^2 + LV^2))
}

calibrate_bi_thresholds <- function(n, G_values, B_threshold, alpha,
                                    seed = 20261012L) {
  set.seed(seed)
  ans <- vapply(G_values, function(G) {
    unname(quantile(replicate(B_threshold, max_bi_mosum_null(n, G)),
                    probs = 1 - alpha, type = 1))
  }, numeric(1))
  names(ans) <- as.character(G_values)
  cat("Bi-MOSUM Monte Carlo thresholds:\n")
  print(ans)
  ans
}

##################################################################
# Original refinement, CI, hotspot and evaluation definitions
# (Sections 7--10 of main script, including B_ci=1000 when supplied)
##################################################################

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

##################################################################
# Additional shared helpers for sensitivity simulations
##################################################################

get_num_cores <- function() {
  num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
  if (is.na(num_cores) || num_cores < 1L) num_cores <- 1L
  num_cores
}

# Returns both G window-size outputs, as in main Simulations_hotspot.R.
fit_both_windows <- function(series_X, series_Y, n, G_values,
                             alpha, epsilon, eta, B_ci) {
  out <- list()
  for (G in G_values) {
    fit <- construct_hotspots(series_X, series_Y, G = G, n = n,
                              alpha = alpha, epsilon = epsilon,
                              eta = eta, B_ci = B_ci)
    out[[paste0("G", G, "_Thrs")]] <- fit$threshold
    out[[paste0("G", G, "_CI")]] <- fit$CI
  }
  out
}

method_map <- function(G_values) {
  out <- c()
  for (rule in c("Thrs", "CI")) {
    for (G in G_values) {
      out[paste0(rule, "_G", G)] <- paste0("G", G, "_", rule)
    }
  }
  out
}

# FHL counts observed time points in the UNION of all intervals.
# This is zero when no hotspots are detected.
hotspot_total_length <- function(hotspots, n) {
  z <- normalize_intervals(hotspots, n)
  if (length(z) == 0L) return(0L)
  sum(vapply(z, function(iv) iv[2L] - iv[1L] + 1L, numeric(1)))
}

# Stationary AR(1) innovations shared across segments; segment-specific
# means and standard deviations are drawn just as in the main DGP.
generate_piecewise_ar1_series <- function(n, true_cp, a,
                                           sigma_min, sigma_max, rho) {
  if (!is.finite(rho) || abs(rho) >= 1) stop("Require |rho| < 1.")
  true_cp <- sort(unique(as.integer(true_cp)))
  boundaries <- c(0L, true_cp, n)
  K <- length(true_cp) + 1L
  mu_vec <- runif(K, -a, a)
  sigma_vec <- runif(K, sigma_min, sigma_max)
  eps <- rnorm(n)
  errors <- numeric(n)
  errors[1L] <- eps[1L]
  if (n > 1L) {
    for (t in 2:n) {
      errors[t] <- rho * errors[t - 1L] + sqrt(1 - rho^2) * eps[t]
    }
  }
  x <- numeric(n)
  for (k in seq_len(K)) {
    idx <- (boundaries[k] + 1L):boundaries[k + 1L]
    x[idx] <- mu_vec[k] + sigma_vec[k] * errors[idx]
  }
  x
}

# Like Library.R::conv_smooth(), but takes an explicit vector of
# uniform draws. This allows all M transformations to be generated
# first and keeps the CI bootstrap random stream fixed across M.
transform_ordinal_with_uniforms <- function(y, uniforms) {
  n <- length(y)
  if (length(uniforms) != n || anyNA(y) || any(!y %in% 1:5)) {
    stop("Invalid ordinal values or uniforms.")
  }
  frequencies <- tabulate(y, nbins = 5L) / n
  bin <- c(0, cumsum(frequencies))
  u <- bin[y] + uniforms * frequencies[y]
  qnorm(u)
}

# Returns one Jaccard similarity for each pair of hotspot sets.
# J(empty,empty)=1, and J(empty,nonempty)=0.
hotspot_jaccard <- function(A, B, n) {
  mask_A <- rep(FALSE, n)
  mask_B <- rep(FALSE, n)
  for (iv in normalize_intervals(A, n)) mask_A[iv[1L]:iv[2L]] <- TRUE
  for (iv in normalize_intervals(B, n)) mask_B[iv[1L]:iv[2L]] <- TRUE
  union_size <- sum(mask_A | mask_B)
  if (union_size == 0L) return(1)
  sum(mask_A & mask_B) / union_size
}

average_pairwise_jaccard <- function(hotspot_sets, n) {
  M <- length(hotspot_sets)
  if (M < 2L) stop("At least two transformations required.")
  pairs <- utils::combn(M, 2L)
  mean(vapply(seq_len(ncol(pairs)), function(k) {
    hotspot_jaccard(hotspot_sets[[pairs[1L, k]]],
                    hotspot_sets[[pairs[2L, k]]], n)
  }, numeric(1)))
}

# Set up a PSOCK cluster exactly as the original foreach/doRNG flow.
# All common functions and any run-specific function are made available.
start_hotspot_cluster <- function(num_cores, additional_functions = character()) {
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  common <- c("bi_thresholds", "refinements_fixed", "case_par",
              "generate_piecewise_series", "get_true_change_points",
              "normalize_intervals", "points_to_intervals",
              "union_interval_sets", "intersect_interval_sets",
              "clean_cpts", "local_argmax", "joint_mosum_pointwise_ci",
              "bivariate_mosum_pointwise_ci", "construct_hotspots",
              "evaluate_hotspot", "fit_both_windows", "method_map",
              "hotspot_total_length", "generate_piecewise_ar1_series",
              "transform_ordinal_with_uniforms", "hotspot_jaccard",
              "average_pairwise_jaccard")
  parallel::clusterExport(cl, varlist = unique(c(common, additional_functions)),
                          envir = .GlobalEnv)
  parallel::clusterEvalQ(cl, {
    source("Library.R")
    refinements <- refinements_fixed
    NULL
  })
  cl
}

# Fail loudly if any replication fails, rather than silently dropping it.
check_replications <- function(list_iter, iter_max, description) {
  if (length(list_iter) != iter_max) {
    stop("Expected ", iter_max, " replications, got ", length(list_iter),
         ": ", description)
  }
  failed <- vapply(list_iter, inherits, logical(1), what = "error")
  if (any(failed)) {
    stop(sum(failed), " of ", iter_max, " replications failed: ",
         description, ". First error: ",
         conditionMessage(list_iter[[which(failed)[1L]]]))
  }
  invisible(TRUE)
}
