##################################################################
# Data Application (local, serial execution; NO temporal shifts)
# Based on the original Applications.R.
# Inputs : Library.R, hotspot_4subjects_data.RData
# Outputs: data_appl.RData, data_appl_table.csv,
#          data_appl_table.tex, figure/hotspot_<subject>.pdf
##################################################################

source("Library.R")
library(ggplot2)
library(patchwork)

# --------------------------------------------------------------- #
# 1. Settings
# --------------------------------------------------------------- #
run_analysis <- TRUE    # FALSE: regenerate figures/table from saved data_appl.RData
set.seed(1234)
alpha <- 0.05
epsilon <- 0.2
eta <- 0.2
N_boot <- 1000L
B_threshold <- 1000L
G_values <- c(20L, 40L)

# G values corresponding to the current manuscript figure captions.
main_figure_G <- c(sternrelief15 = 40L,
                   sternrelief16 = 20L,
                   sternrelief5  = 40L,
                   sternrelief6  = 40L)

load("hotspot_4subjects_data.RData")
stopifnot(exists("df_list"), length(df_list) == 4L)

passive_names <- c("step_count", "conversation_percent", "sleep_duration")
passive_short <- c(step_count = "SC", conversation_percent = "CT",
                   sleep_duration = "SD")
passive_ylabels <- c(step_count = "Daily step count",
                     conversation_percent = "Conversation time (%)",
                     sleep_duration = "Sleep duration (seconds)")
passive_colors <- c(step_count = "purple", conversation_percent = "orange",
                    sleep_duration = "navy")

# --------------------------------------------------------------- #
# 2. Fixed Monte Carlo calibration, once per (n,G).
#    This prevents repeated 1000-replicate calibrations inside every
#    bivariate_mosum() call. The detection rule and G remain unchanged.
#    Does NOT modify Library.R on disk.
# --------------------------------------------------------------- #
refinements_original <- refinements
thresholds <- list()

if (run_analysis) {
  for (n0 in sort(unique(vapply(df_list, function(d) as.integer(d$n), integer(1))))) {
    for (G0 in G_values) {
      if (2L * G0 >= n0) {
        stop(sprintf("G=%d invalid for n=%d (need 2*G < n).", G0, n0))
      }
      key <- paste(n0, G0, sep = "_")
      # Same deterministic calibration seeds as the prior UAHPC application.
      set.seed(20261012L + 1000L * n0 + G0)
      draws <- replicate(B_threshold, max_bi_mosum_null(n0, G0))
      thresholds[[key]] <- as.numeric(quantile(draws, 1 - alpha, type = 1))
    }
  }
  print(thresholds)
}

refinements <- function(test_stat, type, n, Gl, Gr, alpha, epsilon, eta) {
  if (type == "univ") {
    return(refinements_original(test_stat, type, n, Gl, Gr,
                                alpha, epsilon, eta))
  }
  if (type != "bi" || Gl != Gr) stop("Only equal-window Bi-MOSUM is supported.")
  key <- paste(n, Gl, sep = "_")
  Dn <- thresholds[[key]]
  if (is.null(Dn)) stop("Missing calibrated threshold for n,G: ", key)

  idx_cand <- which(is.finite(test_stat) & abs(test_stat) >= Dn)
  cpt_eps <- integer(0)
  cpt_eta <- integer(0)
  if (length(idx_cand)) {
    grp <- cumsum(c(TRUE, diff(idx_cand) != 1L))
    runs <- split(idx_cand, grp)
    for (idx in runs) {
      if (max(idx) - min(idx) + 2L >= epsilon * Gl) {
        cpt_eps <- c(cpt_eps, idx[which.max(abs(test_stat[idx]))])
      }
    }
    band_l <- floor(eta * Gl)
    band_r <- floor(eta * Gr)
    for (idx in idx_cand) {
      band <- seq.int(max(1L, idx - band_l), min(n, idx + band_r))
      if (abs(test_stat[idx]) >= max(abs(test_stat[band]), na.rm = TRUE)) {
        cpt_eta <- c(cpt_eta, idx)
      }
    }
  }
  list(Dn = Dn, idx_cand = idx_cand,
       cpt_eps = if (length(cpt_eps)) sort(unique(cpt_eps)) else NULL,
       cpt_eta = if (length(cpt_eta)) sort(unique(cpt_eta)) else NULL)
}

# --------------------------------------------------------------- #
# 3. Small helpers (original hotspot definitions retained)
# --------------------------------------------------------------- #
get_ptwise <- function(ci) {
  if (is.null(ci)) return(list())
  if (is.null(ci$ptwise)) return(list())
  ci$ptwise
}

# Works for the list-of-pairs output of the existing Library.R.
interval_matrix <- function(x) {
  out <- matrix(numeric(0), ncol = 2L)
  if (is.null(x) || length(x) == 0L) return(out)
  if (is.matrix(x)) {
    if (ncol(x) != 2L) stop("Intervals must have two columns.")
    return(matrix(as.numeric(x), ncol = 2L))
  }
  if (is.data.frame(x)) return(interval_matrix(as.matrix(x)))
  if (is.numeric(x) && length(x) == 2L) x <- list(x)
  if (!is.list(x)) stop("Unsupported interval representation.")
  pieces <- lapply(x, function(z) {
    if (is.null(z) || length(z) != 2L) return(NULL)
    as.numeric(z)
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) return(out)
  out <- do.call(rbind, pieces)
  matrix(as.numeric(out), ncol = 2L)
}

# Normalize intervals onto the observed day grid for tables and figures.
# A fractional interval [a,b] covers integer days ceiling(a):floor(b).
integer_intervals <- function(x, n) {
  x <- interval_matrix(x)
  if (!nrow(x)) return(x)
  start <- pmax(1L, as.integer(ceiling(pmin(x[, 1], x[, 2]))))
  end <- pmin(n, as.integer(floor(pmax(x[, 1], x[, 2]))))
  keep <- is.finite(start) & is.finite(end) & start <= end
  if (!any(keep)) return(matrix(numeric(0), ncol = 2L))
  intervals <- cbind(start[keep], end[keep])
  intervals <- intervals[order(intervals[, 1], intervals[, 2]), , drop = FALSE]
  merged <- list(intervals[1L, ])
  if (nrow(intervals) >= 2L) {
    for (j in seq.int(2L, nrow(intervals))) {
      last <- length(merged)
      if (intervals[j, 1] <= merged[[last]][2] + 1L) {
        merged[[last]][2] <- max(merged[[last]][2], intervals[j, 2])
      } else {
        merged[[length(merged) + 1L]] <- intervals[j, ]
      }
    }
  }
  matrix(as.integer(do.call(rbind, merged)), ncol = 2L)
}

stress_cp_support <- function(stress_cp, hotspots, n) {
  stress_cp <- sort(unique(as.integer(stress_cp)))
  stress_cp <- stress_cp[is.finite(stress_cp)]
  iv <- integer_intervals(hotspots, n)
  supported <- if (length(stress_cp) && nrow(iv)) {
    sum(vapply(stress_cp, function(cp) {
      any(cp >= iv[, 1] & cp <= iv[, 2])
    }, logical(1)))
  } else 0L
  total <- length(stress_cp)
  list(count = as.integer(supported), total = as.integer(total),
       proportion = if (total) supported / total else NA_real_)
}

as_stress <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  x <- as.numeric(x)
  if (anyNA(x) || !all(x %in% 1:5)) {
    stop("Expected imputed original stress responses coded 1,2,3,4,5.")
  }
  x
}

# --------------------------------------------------------------- #
# 4. Original analysis: 4 participants x 3 sensing variables x 2 G
#    Stress-only Y confidence intervals are computed once per (subject,G).
#    No temporal-shift controls are run.
# --------------------------------------------------------------- #
if (run_analysis) {
  result_list <- list()
  cnt <- 0L
  for (i in seq_along(df_list)) {
    subj <- as.character(df_list[[i]]$subj)
    dat <- df_list[[i]]$dat
    n <- as.integer(df_list[[i]]$n)
    stopifnot(length(dat$stress) == n)

    # Use the same subject-specific stress transformation seed as the
    # previous HPC application. This transformation is reused for all G/X.
    set.seed(100000L + i)
    stress <- conv_smooth(as_stress(dat$stress), n, scale = 1:5)
    ci_y_cache <- list()
    cp_y_cache <- list()

    for (j in seq_along(passive_names)) {
      passive_name <- passive_names[[j]]
      passive <- as.numeric(dat[[passive_name]])
      if (length(passive) != n || any(!is.finite(passive))) {
        stop("Invalid passive sensing data: ", subj, " / ", passive_name)
      }

      for (G in G_values) {
        cnt <- cnt + 1L
        # Reproducible bootstrap sampling for this subject/variable/G.
        set.seed(200000L + 1000L * i + 100L * j + G)
        MOSUM <- bivariate_mosum(passive, stress, n = n, Gl = G, Gr = G,
                                 alpha = alpha, epsilon = epsilon, eta = eta)

        # Y is Joint-MOSUM applied to the stress series alone.
        key <- as.character(G)
        cp_y <- sort(unique(as.integer(MOSUM$Y_both$cpt_eta)))
        if (!(key %in% names(ci_y_cache))) {
          CI_Y <- if (length(cp_y)) {
            joint_mosum_ci(stress, refinement = MOSUM$Y_both,
                           method = "eta", n = n, Gl = G, Gr = G,
                           N = N_boot, alpha = alpha,
                           epsilon = epsilon, eta = eta)
          } else NULL
          ci_y_cache[key] <- list(CI_Y)  # list(NULL) is retained, unlike [[<- NULL
          cp_y_cache[key] <- list(cp_y)
        } else {
          CI_Y <- ci_y_cache[[key]]
          if (!identical(cp_y, cp_y_cache[[key]])) {
            stop("Stress-only change points unexpectedly depend on passive series.")
          }
        }

        # Original four cross-series Bi-MOSUM CIs.
        moments <- c(YX = "mean_mean", YX2 = "mean_var",
                     Y2X = "var_mean", Y2X2 = "var_var")
        ci_list <- lapply(names(moments), function(z) {
          if (!length(MOSUM[[z]]$cpt_eta)) return(NULL)
          bivariate_mosum_ci(stress, passive,
                             bi_refinement = MOSUM[[z]],
                             moment = moments[[z]], method = "eta",
                             n = n, Gl = G, Gr = G,
                             N = N_boot, alpha = alpha,
                             epsilon = epsilon, eta = eta)
        })
        names(ci_list) <- names(moments)

        # Keep exactly the original hotspot constructions.
        paired_ptwise <- do.call(c, lapply(ci_list, get_ptwise))
        hotspot_ci <- if (length(paired_ptwise) && length(get_ptwise(CI_Y))) {
          get_union_of_intersects(paired_ptwise, get_ptwise(CI_Y),
                                  type = "ptwise")
        } else list()

        paired_threshold_points <- sort(unique(unlist(lapply(
          names(moments), function(z) MOSUM[[z]]$idx_cand), use.names = FALSE)))
        common <- intersect(paired_threshold_points, MOSUM$Y_both$idx_cand)
        hotspot_thrs <- if (length(common)) get_intervals(sort(common)) else list()

        result_list[[cnt]] <- list(
          subj = subj, passive = passive_name, G = G, n = n,
          MOSUM = MOSUM,
          CI_Y = CI_Y,
          CI_YX = ci_list$YX, CI_YX2 = ci_list$YX2,
          CI_Y2X = ci_list$Y2X, CI_Y2X2 = ci_list$Y2X2,
          hotspot_ci = hotspot_ci, hotspot_thrs = hotspot_thrs,
          stress_only_cp = cp_y,
          stress_support_thrs = stress_cp_support(cp_y, hotspot_thrs, n),
          stress_support_ci = stress_cp_support(cp_y, hotspot_ci, n)
        )
        cat(sprintf("[%02d/24] %s / %s / G=%d complete\n",
                    cnt, subj, passive_name, G))
      }
    }
    # Checkpoint: finished participants are retained if local R is interrupted.
    save(result_list, thresholds, file = "data_appl_partial.RData")
  }
  stopifnot(length(result_list) == 24L)
  save(result_list, thresholds, file = "data_appl.RData")
  cat("Saved data_appl.RData\n")
} else {
  load("data_appl.RData")
  stopifnot(exists("result_list"), length(result_list) == 24L)
}

# --------------------------------------------------------------- #
# 5. Table: original CP columns, hotspot intervals, AND stress-only
#    CP support; one row per (subject, G, passive variable).
# --------------------------------------------------------------- #
fmt_cps <- function(cp) {
  if (is.null(cp) || !length(cp)) return("--")
  cp <- sort(unique(as.integer(round(cp))))
  if (!length(cp)) return("--")
  paste0("$\\{", paste(cp, collapse = ","), "\\}$")
}
fmt_intervals <- function(hotspots, n) {
  iv <- integer_intervals(hotspots, n)
  if (!nrow(iv)) return("--")
  blocks <- apply(iv, 1L, function(z) sprintf("[%d,%d]", z[1], z[2]))
  paste0("$\\{", paste(blocks, collapse = ","), "\\}$")
}
fmt_support <- function(s) {
  if (s$total == 0L) return("--")
  sprintf("%d/%d", s$count, s$total)
}

rows <- lapply(result_list, function(z) {
  cp_y <- z$stress_only_cp
  if (is.null(cp_y)) cp_y <- z$MOSUM$Y_both$cpt_eta
  ss_thrs <- stress_cp_support(cp_y, z$hotspot_thrs, z$n)
  ss_ci <- stress_cp_support(cp_y, z$hotspot_ci, z$n)
  data.frame(
    Participants = as.character(z$subj), G = as.integer(z$G),
    Variable = unname(passive_short[[z$passive]]),
    Y = fmt_cps(cp_y),
    YX = fmt_cps(z$MOSUM$YX$cpt_eta),
    YX2 = fmt_cps(z$MOSUM$YX2$cpt_eta),
    Y2X = fmt_cps(z$MOSUM$Y2X$cpt_eta),
    Y2X2 = fmt_cps(z$MOSUM$Y2X2$cpt_eta),
    Thrs = fmt_intervals(z$hotspot_thrs, z$n),
    CI = fmt_intervals(z$hotspot_ci, z$n),
    Support_Thrs = fmt_support(ss_thrs),
    Support_CI = fmt_support(ss_ci),
    Support_Thrs_Proportion = ss_thrs$proportion,
    Support_CI_Proportion = ss_ci$proportion,
    stringsAsFactors = FALSE
  )
})
data_appl_table <- do.call(rbind, rows)
subject_order <- c("sternrelief15", "sternrelief16", "sternrelief5", "sternrelief6")
data_appl_table <- data_appl_table[order(
  match(data_appl_table$Participants, subject_order),
  data_appl_table$G, match(data_appl_table$Variable, c("SC", "CT", "SD"))), ]
rownames(data_appl_table) <- NULL
write.csv(data_appl_table, "data_appl_table.csv", row.names = FALSE)

# Generate the manuscript's full 24-row LaTeX table.
tex <- c(
  "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",
  "\\begin{table}[h]",
  "\\centering",
  "\\resizebox{0.8\\columnwidth}{!}{%",
  "\\begin{tabular}{cccccccccccc}",
  "\\hline",
  "& & & \\multicolumn{5}{c}{Change points} & \\multicolumn{2}{c}{Hotspots} & \\multicolumn{2}{c}{Stress-only CP support} \\\\",
  "\\cline{4-8}\\cline{9-10}\\cline{11-12}",
  "Participants & $G$ & Variable & $Y$ & $YX$ & $YX^2$ & $Y^2X$ & $Y^2X^2$ & Thrs & CI & Thrs & CI \\\\",
  "\\hline\\hline"
)
for (subject in subject_order) {
  for (G in G_values) {
    subset_rows <- data_appl_table[
      data_appl_table$Participants == subject & data_appl_table$G == G, , drop = FALSE]
    stopifnot(nrow(subset_rows) == 3L)
    for (j in seq_len(3L)) {
      z <- subset_rows[j, ]
      subject_cell <- if (G == G_values[1] && j == 1L) {
        paste0("\\multirow{6}{*}{", subject, "}")
      } else ""
      G_cell <- if (j == 1L) sprintf("\\multirow{3}{*}{$%d$}", G) else ""
      values <- unlist(z[c("Variable", "Y", "YX", "YX2", "Y2X", "Y2X2",
                            "Thrs", "CI", "Support_Thrs", "Support_CI")], use.names = FALSE)
      tex <- c(tex, paste0(subject_cell, " & ", G_cell, " & ",
                           paste(values, collapse = " & "), " \\\\"))
      if (j < 3L) tex <- c(tex, "\\cline{3-3}")
    }
    if (G == G_values[1]) tex <- c(tex, "\\cline{2-3}")
  }
  tex <- c(tex, if (subject == tail(subject_order, 1L)) "\\hline" else "\\hline\\hline")
}
tex <- c(tex,
  "\\end{tabular}%",
  "}",
  "\\caption{Estimated change points and identified hotspots in the real-data application. The five change-point columns report the stress-only Joint-MOSUM estimates ($Y$) and the four cross-series Bi-MOSUM estimates ($YX$, $YX^2$, $Y^2X$, and $Y^2X^2$). The hotspot columns report the thresholding (Thrs) and confidence-interval (CI) regions. Stress-only CP support is the number of stress-only change points contained in the paired hotspot divided by the total number of stress-only change points; -- indicates that no stress-only change points were detected, or no estimated points/intervals for other columns. SC, CT, and SD denote daily step count, conversation time, and sleep duration, respectively. The significance level is $\\alpha=0.05$.}",
  "\\label{tab:data_appl}",
  "\\end{table}",
  "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
writeLines(tex, "data_appl_table.tex", useBytes = TRUE)
cat("Saved data_appl_table.csv and data_appl_table.tex\n")

# --------------------------------------------------------------- #
# 6. Four manuscript figures with numeric y-axis tick labels/units.
#    Use ORIGINAL integer Likert responses, not cut(..., breaks=5).
# --------------------------------------------------------------- #
dir.create("figure", showWarnings = FALSE)

shade_df <- function(thrs, ci, n) {
  make <- function(x, type) {
    iv <- integer_intervals(x, n)
    if (!nrow(iv)) return(data.frame(xmin = numeric(), xmax = numeric(),
                                     type = character()))
    data.frame(xmin = iv[, 1] - 0.5, xmax = iv[, 2] + 0.5,
               type = rep(type, nrow(iv)))
  }
  out <- rbind(make(thrs, "Thrs"), make(ci, "CI"))
  out$type <- factor(out$type, levels = c("Thrs", "CI"))
  out
}

get_result <- function(subj, passive_name, G) {
  idx <- which(vapply(result_list, function(z) {
    identical(as.character(z$subj), as.character(subj)) &&
      identical(as.character(z$passive), as.character(passive_name)) &&
      identical(as.integer(z$G), as.integer(G))
  }, logical(1)))
  if (length(idx) != 1L) stop("Missing or duplicate result: ", subj,
                              " / ", passive_name, " / G=", G)
  result_list[[idx]]
}

mdpi_number_labels <- function(x) {
  labels <- scales::label_number(big.mark = "")(x)
  large <- is.finite(x) & abs(x) >= 10000
  labels[large] <- scales::comma(x[large])
  labels
}

plot_pair <- function(subject_data, passive_name, G, bottom = FALSE) {
  subj <- as.character(subject_data$subj)
  dat <- subject_data$dat
  n <- as.integer(subject_data$n)
  fit <- get_result(subj, passive_name, G)
  shade <- shade_df(fit$hotspot_thrs, fit$hotspot_ci, n)
  days <- seq_len(n)
  passive_value <- as.numeric(dat[[passive_name]])
  if (passive_name == "conversation_percent") passive_value <- 100 * passive_value
  stress_value <- as_stress(dat$stress)
  steps <- if (n <= 100L) 10L else 20L
  x_title <- if (bottom) "Day" else NULL

  shaded <- list(
    geom_rect(data = shade,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                  fill = type, alpha = type), inherit.aes = FALSE),
    scale_fill_manual(name = "Hotspots", values = c(Thrs = "red", CI = "blue"),
                      breaks = c("Thrs", "CI"), drop = FALSE),
    scale_alpha_manual(values = c(Thrs = 0.35, CI = 0.22), guide = "none"),
    guides(fill = guide_legend(override.aes = list(alpha = c(0.35, 0.22))))
  )
  base_theme <- theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(color = "black"),
          axis.ticks.y = element_line(color = "grey50"),
          legend.position = "bottom", panel.grid.minor = element_blank())
  xscale <- scale_x_continuous(breaks = seq(0, n, by = steps),
                               limits = c(0.5, n + 0.5), expand = c(0, 0))

  p1 <- ggplot(data.frame(day = days, value = passive_value), aes(day, value)) +
    shaded + geom_line(color = unname(passive_colors[[passive_name]]), size = 0.5) +
    xscale + scale_y_continuous(labels = mdpi_number_labels) +
    labs(x = x_title, y = unname(passive_ylabels[[passive_name]])) + base_theme

  p2 <- ggplot(data.frame(day = days, value = stress_value), aes(day, value)) +
    shaded + geom_step(color = "black", size = 0.5) +
    xscale + scale_y_continuous(breaks = 1:5, limits = c(1, 5)) +
    labs(x = x_title, y = "Self-reported stress (1–5)") + base_theme

  list(passive = p1, stress = p2)
}

for (d in df_list) {
  subj <- as.character(d$subj)
  G <- unname(main_figure_G[[subj]])
  if (is.null(G) || !length(G)) stop("No main G configured for ", subj)
  panels <- lapply(seq_along(passive_names), function(j) {
    plot_pair(d, passive_names[[j]], G, bottom = (j == 3L))
  })
  combined <- ((panels[[1]]$passive | panels[[1]]$stress) /
               (panels[[2]]$passive | panels[[2]]$stress) /
               (panels[[3]]$passive | panels[[3]]$stress)) +
    plot_layout(guides = "collect") + plot_annotation(title = subj) &
    theme(legend.position = "bottom")
  outfile <- file.path("figure", paste0("hotspot_", subj, ".pdf"))
  ggsave(outfile, plot = combined, width = 13, height = 11, units = "in")
  cat("Saved ", outfile, " (G=", G, ")\n", sep = "")
}

cat("Finished: data_appl.RData, table CSV/TEX, four PDFs. No temporal shifts.\n")
