##################################################################
# Data Application: Temporal-shift negative control (LOCAL / SERIAL)
# Based on Applications_local_updated.R; THRESHOLDING ONLY.
#
# Required files in the working directory:
#   Library.R
#   hotspot_4subjects_data.RData
#   data_appl.RData   <-- from Applications_local_updated.R
#
# Output files (the original data_appl.RData is NEVER overwritten):
#   data_appl_shift_raw.csv             : all nonzero circular shifts
#   data_appl_shift_summary.csv         : 24 observed-vs-shifted results
#   data_appl_shift_table.tex           : manuscript tab:shift_control
#   data_appl_shift_table_with_coverage.tex : optional expanded table
#
# Resume: rerun this script after interruption; completed shifts are
#         read from data_appl_shift_raw.csv and are not recalculated.
##################################################################

source("Library.R")

# --------------------------------------------------------------- #
# 1. Settings
# --------------------------------------------------------------- #
run_shifts <- TRUE        # FALSE: recreate final tables using saved raw CSV
checkpoint_every <- 10L  # save progress after every 10 shifted alignments
G_values <- c(20L, 40L)
passive_names <- c("step_count", "conversation_percent", "sleep_duration")
passive_short <- c(step_count = "SC", conversation_percent = "CT",
                   sleep_duration = "SD")
subject_order <- c("sternrelief15", "sternrelief16", "sternrelief5",
                   "sternrelief6")
alpha <- 0.05

required <- c("hotspot_4subjects_data.RData", "data_appl.RData")
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing input file(s): ", paste(missing, collapse = ", "))

load("hotspot_4subjects_data.RData")
load("data_appl.RData")
if (!exists("df_list") || !exists("result_list") || length(result_list) != 24L) {
  stop("Expected df_list and 24 saved observed analyses in data_appl.RData.")
}

# --------------------------------------------------------------- #
# 2. Helpers: SAME thresholding hotspot rule as Applications.R.
#    Reuse Y's stored tracker and the observed (n,G) threshold.
#    No CI bootstrap, no Monte Carlo recalibration, no parallel R.
# --------------------------------------------------------------- #
get_observed <- function(subj, passive_name, G) {
  hit <- which(vapply(result_list, function(z) {
    identical(as.character(z$subj), subj) &&
      identical(as.character(z$passive), passive_name) &&
      identical(as.integer(z$G), as.integer(G))
  }, logical(1)))
  if (length(hit) != 1L) {
    stop("Cannot find exactly one observed fit: ", subj, ", ", passive_name, ", G=", G)
  }
  result_list[[hit]]
}

circular_shift <- function(x, shift) {
  n <- length(x)
  shift <- as.integer(shift) %% n
  if (shift == 0L) return(x)
  c(tail(x, shift), head(x, n - shift))
}

# Return a logical indicator of all discrete time points covered by
# the UNION of the estimated hotspot intervals (count once only).
hotspot_mask <- function(intervals, n) {
  covered <- rep(FALSE, n)
  if (is.null(intervals) || !length(intervals)) return(covered)
  if (is.matrix(intervals)) {
    if (ncol(intervals) != 2L) stop("Expected two-column hotspot interval matrix.")
    intervals <- lapply(seq_len(nrow(intervals)), function(i) intervals[i, ])
  }
  for (interval in intervals) {
    if (length(interval) != 2L || any(!is.finite(interval))) {
      stop("Unexpected hotspot interval format.")
    }
    left <- max(1L, as.integer(ceiling(min(interval))))
    right <- min(n, as.integer(floor(max(interval))))
    if (left <= right) covered[seq.int(left, right)] <- TRUE
  }
  covered
}

get_support <- function(mask, cp) {
  cp <- sort(unique(as.integer(cp)))
  cp <- cp[is.finite(cp)]
  if (!length(cp)) return(NA_real_)  # no stress-only CP -> undefined
  mean(mask[cp])
}

# Extract EXACTLY the same precomputed threshold as observed analysis.
# data_appl.RData from Applications_local_updated.R includes `thresholds`.
get_threshold <- function(fit, n, G) {
  key <- paste(n, G, sep = "_")
  if (exists("thresholds", inherits = TRUE) &&
      !is.null(thresholds[[key]])) {
    Dn <- as.numeric(thresholds[[key]])
  } else {
    # The threshold also survives in each saved MOSUM refinement.
    Dn <- as.numeric(fit$MOSUM$YX$Dn)
  }
  if (length(Dn) != 1L || !is.finite(Dn) || Dn <= 0) {
    stop("No valid saved threshold for n=", n, ", G=", G)
  }
  Dn
}

# Faster but equivalent to bivariate_mosum(passive, stress, ...)
# for the thresholding rule: Y_track is fixed under all shifts.
# All four bivariate detector distances are still recomputed.
fit_shift_thrs <- function(shifted_passive, saved_Y_track, saved_Y_idx, G, n, Dn) {
  X_track <- single_tracker(shifted_passive, n = n, Gl = G, Gr = G)
  btrack <- bivariate_tracker(saved_Y_track, X_track, n = n, Gl = G, Gr = G)
  cross_dist <- list(btrack$Dist_YX, btrack$Dist_YX2,
                     btrack$Dist_Y2X, btrack$Dist_Y2X2)
  cross_idx <- sort(unique(unlist(lapply(cross_dist, function(d) {
    which(is.finite(d) & abs(d) >= Dn)
  }), use.names = FALSE)))
  shared_idx <- intersect(cross_idx, saved_Y_idx)
  if (length(shared_idx)) get_intervals(sort(shared_idx)) else list()
}

# Save raw results periodically, so a long local serial run is resumable.
raw_file <- "data_appl_shift_raw.csv"
empty_raw <- function() {
  data.frame(subject = character(), G = integer(), passive = character(),
             n = integer(), shift = integer(), n_stress_cp = integer(),
             support = numeric(), hotspot_days = integer(),
             hotspot_fraction = numeric(), stringsAsFactors = FALSE)
}
shift_raw <- if (file.exists(raw_file)) {
  read.csv(raw_file, stringsAsFactors = FALSE)
} else empty_raw()
if (!identical(names(shift_raw), names(empty_raw()))) {
  stop("Existing raw CSV has unexpected columns: ", raw_file)
}

# --------------------------------------------------------------- #
# 3. All nonzero circular shifts: 4 participants x 3 variables x 2 G
#    Each shift s = 1,...,n-1 uses the SAME transformed Y, window
#    size, threshold, and four cross-series detector combinations.
# --------------------------------------------------------------- #
if (run_shifts) {
  for (subj in subject_order) {
    dat_ix <- which(vapply(df_list, function(z) {
      identical(as.character(z$subj), subj)
    }, logical(1)))
    if (length(dat_ix) != 1L) stop("Missing data for ", subj)
    d <- df_list[[dat_ix]]
    n <- as.integer(d$n)

    for (G in G_values) {
      if (2L * G >= n) stop("Invalid window: n=", n, ", G=", G)
      for (passive_name in passive_names) {
        fit <- get_observed(subj, passive_name, G)
        x <- as.numeric(d$dat[[passive_name]])
        if (length(x) != n || any(!is.finite(x))) {
          stop("Missing/invalid passive data for ", subj, ", ", passive_name)
        }
        if (is.null(fit$MOSUM$Y_track$series)) {
          stop("Saved observed fit lacks the stress tracker. Rerun Applications_local_updated.R.")
        }
        Y_track <- fit$MOSUM$Y_track
        Y_idx <- as.integer(fit$MOSUM$Y_both$idx_cand)
        Y_cp <- as.integer(fit$MOSUM$Y_both$cpt_eta)
        Y_cp <- sort(unique(Y_cp[is.finite(Y_cp)]))
        Dn <- get_threshold(fit, n, G)
        obs_mask <- hotspot_mask(fit$hotspot_thrs, n)

        # IMPORTANT validation: the optimized computation at shift=0 must
        # match the hotspots obtained in the ORIGINAL applications analysis.
        check0 <- fit_shift_thrs(x, Y_track, Y_idx, G, n, Dn)
        if (!identical(hotspot_mask(check0, n), obs_mask)) {
          stop("Shift=0 disagrees with saved observed Thrs hotspots: ",
               subj, ", ", passive_name, ", G=", G,
               ". Check whether the saved Library.R / thresholds differ.")
        }

        done <- shift_raw$shift[
          shift_raw$subject == subj & shift_raw$G == G &
          shift_raw$passive == passive_name]
        if (length(done) && (anyNA(done) || anyDuplicated(done) ||
                             any(done < 1L | done >= n))) {
          stop("Invalid checkpoint entries for ", subj, ", ", passive_name, ", G=", G)
        }
        remaining <- setdiff(seq_len(n - 1L), done)
        cat(sprintf("\n%s | G=%d | %s | remaining=%d/%d | S0=%s\n",
                    subj, G, passive_name, length(remaining), n - 1L,
                    if (length(Y_cp)) sprintf("%.3f", get_support(obs_mask, Y_cp))
                    else "NA (no stress-only CP)"))

        for (idx in seq_along(remaining)) {
          s <- remaining[[idx]]
          hot <- fit_shift_thrs(circular_shift(x, s), Y_track, Y_idx, G, n, Dn)
          mask <- hotspot_mask(hot, n)
          row <- data.frame(subject = subj, G = G, passive = passive_name,
                            n = n, shift = s, n_stress_cp = length(Y_cp),
                            support = get_support(mask, Y_cp),
                            hotspot_days = sum(mask),
                            hotspot_fraction = mean(mask),
                            stringsAsFactors = FALSE)
          shift_raw <- rbind(shift_raw, row)
          if (idx %% checkpoint_every == 0L || idx == length(remaining)) {
            write.csv(shift_raw, raw_file, row.names = FALSE, na = "NA")
            cat(sprintf("  saved %d/%d remaining shifts (%d raw rows total)\n",
                        idx, length(remaining), nrow(shift_raw)))
          }
        }
      }
    }
  }
}

# --------------------------------------------------------------- #
# 4. Summarize the observed versus all nonzero shifted alignments.
# --------------------------------------------------------------- #
if (file.exists(raw_file)) {
  shift_raw <- read.csv(raw_file, stringsAsFactors = FALSE)
} else stop("No shift results available. Set run_shifts <- TRUE.")

q <- function(x, p) unname(quantile(x, probs = p, names = FALSE, type = 7))
summary_rows <- list()
count <- 0L
for (subj in subject_order) {
  d <- df_list[[which(vapply(df_list, function(z) {
    identical(as.character(z$subj), subj)
  }, logical(1)))[[1L]]]]
  n <- as.integer(d$n)
  for (G in G_values) {
    for (passive_name in passive_names) {
      fit <- get_observed(subj, passive_name, G)
      Y_cp <- sort(unique(as.integer(fit$MOSUM$Y_both$cpt_eta)))
      Y_cp <- Y_cp[is.finite(Y_cp)]
      observed_mask <- hotspot_mask(fit$hotspot_thrs, n)
      S0 <- get_support(observed_mask, Y_cp)
      C0 <- mean(observed_mask)

      z <- shift_raw[shift_raw$subject == subj & shift_raw$G == G &
                       shift_raw$passive == passive_name, , drop = FALSE]
      if (nrow(z) != n - 1L || !setequal(as.integer(z$shift), seq_len(n - 1L)) ||
          anyDuplicated(z$shift)) {
        stop(sprintf("Incomplete controls: %s / G=%d / %s (%d of %d shifts).",
                     subj, G, passive_name, nrow(z), n - 1L))
      }
      if (any(z$n != n) || any(z$n_stress_cp != length(Y_cp))) {
        stop("Raw results appear to belong to a different observed dataset: ", subj)
      }
      if (length(Y_cp)) {
        z_support <- z$support
        if (anyNA(z_support)) stop("Missing support value for ", subj)
        S_med <- median(z_support)
        S_q25 <- q(z_support, 0.25)
        S_q75 <- q(z_support, 0.75)
        E_shift <- mean(z_support >= S0)
      } else {
        S_med <- S_q25 <- S_q75 <- E_shift <- NA_real_
      }
      coverage <- z$hotspot_fraction
      if (anyNA(coverage)) stop("Missing coverage value for ", subj)
      count <- count + 1L
      summary_rows[[count]] <- data.frame(
        Participants = subj, G = G,
        Variable = unname(passive_short[[passive_name]]),
        n = n, n_stress_cp = length(Y_cp), n_shifts = nrow(z),
        S0 = S0,
        Shifted_support_median = S_med,
        Shifted_support_q25 = S_q25,
        Shifted_support_q75 = S_q75,
        E_shift = E_shift,
        C0 = C0,
        Shifted_coverage_median = median(coverage),
        Shifted_coverage_q25 = q(coverage, 0.25),
        Shifted_coverage_q75 = q(coverage, 0.75),
        E_shift_coverage = mean(coverage >= C0),
        stringsAsFactors = FALSE)
    }
  }
}
shift_summary <- do.call(rbind, summary_rows)
stopifnot(nrow(shift_summary) == 24L)
write.csv(shift_summary, "data_appl_shift_summary.csv", row.names = FALSE)

# --------------------------------------------------------------- #
# 5. Generate LaTeX tab:shift_control (same structure as manuscript)
#    and an OPTIONAL expanded table with the hotspot coverage data.
# --------------------------------------------------------------- #
fmt_val <- function(z) if (is.na(z)) "--" else sprintf("%.3f", z)
fmt_median_iqr <- function(m, l, u) {
  if (anyNA(c(m, l, u))) "--" else
    sprintf("%.3f (%.3f--%.3f)", m, l, u)
}

write_shift_tex <- function(full_coverage = FALSE) {
  if (!full_coverage) {
    lines <- c(
      "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",
      "\\begin{table}[h]", "\\centering", "\\resizebox{0.6\\columnwidth}{!}{%",
      "\\begin{tabular}{cccccc}", "\\hline",
      "Participant & $G$ & Variable & Observed support & Shifted support & Exceedance \\\\",
      "& & & $S_0$ & Median (IQR) & $E_{\\mathrm{shift}}$ \\\\",
      "\\hline\\hline")
  } else {
    lines <- c(
      "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",
      "\\begin{table}[h]", "\\centering", "\\resizebox{0.9\\columnwidth}{!}{%",
      "\\begin{tabular}{ccccccccc}", "\\hline",
      "& & & \\multicolumn{3}{c}{CP support} & \\multicolumn{3}{c}{Hotspot coverage} \\\\",
      "Participant & $G$ & Variable & Observed $S_0$ & Shifted median (IQR) & $E_{\\mathrm{shift}}$ & Observed $C_0$ & Shifted median (IQR) & $E_{\\mathrm{cov}}$ \\\\",
      "\\hline\\hline")
  }

  for (subj in subject_order) {
    for (G in G_values) {
      for (j in seq_along(passive_names)) {
        v <- unname(passive_short[[passive_names[[j]]]])
        z <- shift_summary[shift_summary$Participants == subj &
                             shift_summary$G == G & shift_summary$Variable == v, ]
        stopifnot(nrow(z) == 1L)
        first <- if (G == G_values[[1L]] && j == 1L) {
          paste0("\\multirow{6}{*}{", subj, "}")
        } else ""
        gcell <- if (j == 1L) sprintf("\\multirow{3}{*}{%d}", G) else ""
        vals <- c(fmt_val(z$S0),
                  fmt_median_iqr(z$Shifted_support_median,
                                 z$Shifted_support_q25, z$Shifted_support_q75),
                  fmt_val(z$E_shift))
        if (full_coverage) {
          vals <- c(vals, fmt_val(z$C0),
                    fmt_median_iqr(z$Shifted_coverage_median,
                                   z$Shifted_coverage_q25, z$Shifted_coverage_q75),
                    fmt_val(z$E_shift_coverage))
        }
        lines <- c(lines, paste0(first, " & ", gcell, " & ", v, " & ",
                                paste(vals, collapse = " & "), " \\\\"))
      }
      if (G == G_values[[1L]]) lines <- c(lines, if (full_coverage) "\\cline{2-9}" else "\\cline{2-6}")
    }
    lines <- c(lines, if (subj == tail(subject_order, 1L))
               "\\hline" else "\\hline\\hline")
  }

  if (full_coverage) {
    caption <- paste0(
      "Temporal-shift negative-control analysis using all nonzero circular shifts ",
      "of the passive sensing series. $S_0$ and $C_0$ denote stress-only ",
      "change-point support and the proportion of time points covered by ",
      "thresholding hotspots under the observed alignment. Shifted columns ",
      "report median (first--third quartile). $E_{\\mathrm{shift}}$ and ",
      "$E_{\\mathrm{cov}}$ denote the proportions of shifted alignments ",
      "with support and coverage at least as large as observed, respectively. ",
      "SC, CT, and SD denote daily step count, conversation time, and sleep ",
      "duration. -- indicates undefined stress-only support when no ",
      "stress-only change points are detected. These exceedance measures ",
      "are descriptive, not calibrated $p$-values.")
    label <- "tab:shift_control_coverage"
  } else {
    caption <- paste0(
      "Results of the temporal-shift negative-control analysis using all ",
      "nonzero circular shifts of the passive sensing series. For each ",
      "participant, window size $G$, and passive sensing variable, $S_0$ ",
      "is the proportion of stress-only change points contained in the ",
      "thresholding hotspots under the observed alignment. Shifted support ",
      "is summarized as median (first--third quartile). ",
      "$E_{\\mathrm{shift}}$ is the proportion of shifted alignments with ",
      "support at least as large as observed. SC, CT, and SD denote daily ",
      "step count, conversation time, and sleep duration. -- indicates ",
      "undefined support when no stress-only change points are detected. ",
      "The exceedance proportion is descriptive, not a calibrated $p$-value.")
    label <- "tab:shift_control"
  }
  lines <- c(lines, "\\end{tabular}%", "}",
             paste0("\\caption{\\textcolor{red}{", caption, "}}"),
             paste0("\\label{", label, "}"), "\\end{table}",
             "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  filename <- if (full_coverage) "data_appl_shift_table_with_coverage.tex" else
    "data_appl_shift_table.tex"
  writeLines(lines, filename, useBytes = TRUE)
  cat("Saved ", filename, "\n", sep = "")
}

write_shift_tex(FALSE)
write_shift_tex(TRUE)
cat("\nFinished: ", nrow(shift_raw), " shifted alignments, ",
    nrow(shift_summary), " observed configurations.\n", sep = "")
cat("Main manuscript table: data_appl_shift_table.tex (tab:shift_control).\n")
cat("Coverage diagnostics: data_appl_shift_summary.csv or expanded LaTeX.\n")
