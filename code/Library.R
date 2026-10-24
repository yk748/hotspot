max_bi_mosum_null <- function(n, G) {
  
  WE <- cumsum(rnorm(n))
  WV <- cumsum(rnorm(n))
  
  LE <- (
    WE[(2 * G + 1):n] -
      2 * WE[(G + 1):(n - G)] +
      WE[1:(n - 2 * G)]
  ) / sqrt(2 * G)
  
  LV <- (
    WV[(2 * G + 1):n] -
      2 * WV[(G + 1):(n - G)] +
      WV[1:(n - 2 * G)]
  ) / sqrt(2 * G)
  
  max(sqrt(LE^2 + LV^2))
}

normalize_intervals <- function(intervals, n) {
  
  if (is.null(intervals) || length(intervals) == 0) {
    return(list())
  }
  
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
  
  ord <- order(vapply(out, `[`, numeric(1), 1))
  out <- out[ord]
  
  merged <- list(out[[1]])
  
  if (length(out) >= 2) {
    for (j in 2:length(out)) {
      
      last <- merged[[length(merged)]]
      cur <- out[[j]]
      
      if (cur[1] <= last[2] + 1L) {
        merged[[length(merged)]] <-
          c(last[1], max(last[2], cur[2]))
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
  
  intervals <- lapply(
    seq_len(length(breaks) - 1L),
    function(j) {
      c(
        x[breaks[j] + 1L],
        x[breaks[j + 1L]]
      )
    }
  )
  
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
        intersections[[length(intersections) + 1L]] <-
          c(left, right)
      }
    }
  }
  
  normalize_intervals(intersections, n)
}


############################################
# Bivariate MOSUM (Bi-MOSUM) tracker:
############################################
bivariate_mosum <- function(xt,yt,n,Gl=20,Gr=20,alpha=0.05,epsilon=0.2,eta=0.4){
  
  # ---------------------------------------------------- #
  # Initial tracking:
  Y_track <- single_tracker(yt,n,Gl,Gr)
  X_track <- single_tracker(xt,n,Gl,Gr)
  bi_track <- bivariate_tracker(Y_track,X_track,n,Gl,Gr)
  
  # ---------------------------------------------------- #
  # Refinement:
  Y_mean <- refinements(Y_track$Test1,type="univ",n,Gl,Gr,alpha,epsilon,eta)
  Y_var <- refinements(Y_track$Test2,type="univ",n,Gl,Gr,alpha,epsilon,eta)
  Y_both <- refinements(Y_track$Dist,type="bi",n,Gl,Gr,alpha,epsilon,eta)
  
  X_mean <- refinements(X_track$Test1,type="univ",n,Gl,Gr,alpha,epsilon,eta)
  X_var <- refinements(X_track$Test2,type="univ",n,Gl,Gr,alpha,epsilon,eta)
  X_both <- refinements(X_track$Dist,type="bi",n,Gl,Gr,alpha,epsilon,eta)
  
  Y_mean_X_mean <- refinements(bi_track$Dist_YX,type="bi",n,Gl,Gr,alpha,epsilon,eta)
  Y_mean_X_var <- refinements(bi_track$Dist_YX2,type="bi",n,Gl,Gr,alpha,epsilon,eta)
  Y_var_X_mean <- refinements(bi_track$Dist_Y2X,type="bi",n,Gl,Gr,alpha,epsilon,eta)
  Y_var_X_var <- refinements(bi_track$Dist_Y2X2,type="bi",n,Gl,Gr,alpha,epsilon,eta)
  
  # Return:
  output <- list()
  output$X_track <- X_track
  output$Y_track <- Y_track
  output$bi_track <- bi_track
  
  output$Y_mean <- Y_mean
  output$Y_var <- Y_var
  output$Y_both <- Y_both
  output$X_mean <- X_mean
  output$X_var <- X_var
  output$X_both <- X_both
  
  output$YX <- Y_mean_X_mean
  output$YX2 <- Y_mean_X_var
  output$Y2X <- Y_var_X_mean
  output$Y2X2 <- Y_var_X_var
  return(output)
}

############################################
# P-value tracker in Bi-MOSUM:
############################################
bivariate_tracker <- function(series1_track,series2_track,n,Gl,Gr){
  
  # ------------------------------------------------------ #
  Cov_YX_bar <- Cov_YX2_bar <- Cov_Y2X_bar <- Cov_Y2X2_bar <- rep(NA,n)
  rho_YX <- rho_YX2 <- rho_Y2X <- rho_Y2X2 <- rep(NA,n)
  Dist_YX <- Dist_YX2 <- Dist_Y2X <- Dist_Y2X2 <- rep(NA,n)
  pval_YX <- pval_YX2 <- pval_Y2X <- pval_Y2X2 <- rep(NA,n)
  
  yt <- series1_track$series
  S2_Y_bar <- series1_track$S2_bar
  V2_Y_bar <- series1_track$V2_bar
  T1_Y <- series1_track$Test1
  T2_Y <- series1_track$Test2
  
  xt <- series2_track$series 
  S2_X_bar <- series2_track$S2_bar
  V2_X_bar <- series2_track$V2_bar
  T1_X <- series2_track$Test1
  T2_X <- series2_track$Test2
  
  for (k in 1:n){
    if (k <= Gl-1){
      
      Y_bar <- mean(yt[1:(Gl+Gr)])
      X_bar <- mean(xt[1:(Gl+Gr)])
      
      # YX
      Cov_YX_bar[k] <- mean( (yt[1:(Gl+Gr)] - Y_bar)*(xt[1:(Gl+Gr)] - X_bar) )
      rho_YX[k] <- Cov_YX_bar[k]/sqrt(S2_Y_bar[k] * S2_X_bar[k])
      Gamma_YX <- diag(1,2); Gamma_YX[1,2] <- Gamma_YX[2,1] <- rho_YX[k]
      Dist_YX[k] <- sqrt(t(c(T1_Y[k],T1_X[k])) %*% solve(Gamma_YX) %*% c(T1_Y[k],T1_X[k]))
      pval_YX[k] <- 1-pchisq(Dist_YX[k]^2,df=2)
      
      # Y2X
      Cov_Y2X_bar[k] <- mean( ((yt[1:(Gl+Gr)] - Y_bar)^2-S2_Y_bar[k])*(xt[1:(Gl+Gr)] - X_bar) )
      rho_Y2X[k] <- Cov_Y2X_bar[k]/sqrt(V2_Y_bar[k] * S2_X_bar[k])
      Gamma_Y2X <- diag(1,2); Gamma_Y2X[1,2] <- Gamma_Y2X[2,1] <- rho_Y2X[k]
      Dist_Y2X[k] <- sqrt(t(c(T2_Y[k],T1_X[k])) %*% solve(Gamma_Y2X) %*% c(T2_Y[k],T1_X[k]))
      pval_Y2X[k] <- 1-pchisq(Dist_Y2X[k]^2,df=2)
      
      # YX2
      Cov_YX2_bar[k] <- mean( (yt[1:(Gl+Gr)] - Y_bar)*((xt[1:(Gl+Gr)] - X_bar)^2-S2_X_bar[k]) )
      rho_YX2[k] <- Cov_YX2_bar[k]/sqrt(S2_Y_bar[k] * V2_X_bar[k])
      Gamma_YX2 <- diag(1,2); Gamma_YX2[1,2] <- Gamma_YX2[2,1] <- rho_YX2[k]
      Dist_YX2[k] <- sqrt(t(c(T1_Y[k],T2_X[k])) %*% solve(Gamma_YX2) %*% c(T1_Y[k],T2_X[k]))
      pval_YX2[k] <- 1-pchisq(Dist_YX2[k]^2,df=2)
      
      # Y2X2
      Cov_Y2X2_bar[k] <- mean( ((yt[1:(Gl+Gr)] - Y_bar)^2-S2_Y_bar[k])*((xt[1:(Gl+Gr)] - X_bar)^2-S2_X_bar[k]) )
      rho_Y2X2[k] <- Cov_Y2X2_bar[k]/sqrt(V2_Y_bar[k] * V2_X_bar[k])
      Gamma_Y2X2 <- diag(1,2); Gamma_Y2X2[1,2] <- Gamma_Y2X2[2,1] <- rho_Y2X2[k]
      Dist_Y2X2[k] <- sqrt(t(c(T2_Y[k],T2_X[k])) %*% solve(Gamma_Y2X2) %*% c(T2_Y[k],T2_X[k]))
      pval_Y2X2[k] <- 1-pchisq(Dist_Y2X2[k]^2,df=2)
      
    }else if (k >= n-Gr+1){
      
      Y_bar <- mean(yt[(n-(Gl+Gr)+1):n])
      X_bar <- mean(xt[(n-(Gl+Gr)+1):n])
      
      # YX
      Cov_YX_bar[k] <- mean( (yt[(n-(Gl+Gr)+1):n] - Y_bar)*(xt[(n-(Gl+Gr)+1):n] - X_bar) )
      rho_YX[k] <- Cov_YX_bar[k]/sqrt(S2_Y_bar[k] * S2_X_bar[k])
      Gamma_YX <- diag(1,2); Gamma_YX[1,2] <- Gamma_YX[2,1] <- rho_YX[k]
      Dist_YX[k] <- sqrt(t(c(T1_Y[k],T1_X[k])) %*% solve(Gamma_YX) %*% c(T1_Y[k],T1_X[k]))
      pval_YX[k] <- 1-pchisq(Dist_YX[k]^2,df=2)
      
      # Y2X
      Cov_Y2X_bar[k] <- mean( ((yt[(n-(Gl+Gr)+1):n] - Y_bar)^2-S2_Y_bar[k])*(xt[(n-(Gl+Gr)+1):n] - X_bar) )
      rho_Y2X[k] <- Cov_Y2X_bar[k]/sqrt(V2_Y_bar[k] * S2_X_bar[k])
      Gamma_Y2X <- diag(1,2); Gamma_Y2X[1,2] <- Gamma_Y2X[2,1] <- rho_Y2X[k]
      Dist_Y2X[k] <- sqrt(t(c(T2_Y[k],T1_X[k])) %*% solve(Gamma_Y2X) %*% c(T2_Y[k],T1_X[k]))
      pval_Y2X[k] <- 1-pchisq(Dist_Y2X[k]^2,df=2)
      
      # YX2
      Cov_YX2_bar[k] <- mean( (yt[(n-(Gl+Gr)+1):n] - Y_bar)*((xt[(n-(Gl+Gr)+1):n] - X_bar)^2-S2_X_bar[k]) )
      rho_YX2[k] <- Cov_YX2_bar[k]/sqrt(S2_Y_bar[k] * V2_X_bar[k])
      Gamma_YX2 <- diag(1,2); Gamma_YX2[1,2] <- Gamma_YX2[2,1] <- rho_YX2[k]
      Dist_YX2[k] <- sqrt(t(c(T1_Y[k],T2_X[k])) %*% solve(Gamma_YX2) %*% c(T1_Y[k],T2_X[k]))
      pval_YX2[k] <- 1-pchisq(Dist_YX2[k]^2,df=2)
      
      # Y2X2
      Cov_Y2X2_bar[k] <- mean( ((yt[(n-(Gl+Gr)+1):n] - Y_bar)^2-S2_Y_bar[k])*((xt[(n-(Gl+Gr)+1):n] - X_bar)^2-S2_X_bar[k]) )
      rho_Y2X2[k] <- Cov_Y2X2_bar[k]/sqrt(V2_Y_bar[k] * V2_X_bar[k])
      Gamma_Y2X2 <- diag(1,2); Gamma_Y2X2[1,2] <- Gamma_Y2X2[2,1] <- rho_Y2X2[k]
      Dist_Y2X2[k] <- sqrt(t(c(T2_Y[k],T2_X[k])) %*% solve(Gamma_Y2X2) %*% c(T2_Y[k],T2_X[k]))
      pval_Y2X2[k] <- 1-pchisq(Dist_Y2X2[k]^2,df=2)
      
    }else{
      
      Y_bar_prv <- mean(yt[(k-Gl+1):k])
      Y_bar_pst <- mean(yt[(k+1):(k+Gr)])
      X_bar_prv <- mean(xt[(k-Gl+1):k])
      X_bar_pst <- mean(xt[(k+1):(k+Gr)])
      
      S2_Y_prv <- mean( (yt[(k-Gl+1):k] - Y_bar_prv)^2 )
      S2_Y_pst <- mean( (yt[(k+1):(k+Gr)] - Y_bar_pst)^2 )
      S2_X_prv <- mean( (xt[(k-Gl+1):k] - X_bar_prv)^2 )
      S2_X_pst <- mean( (xt[(k+1):(k+Gr)] - X_bar_pst)^2 )
      
      # YX
      Cov_YX_prv <- mean((yt[(k-Gl+1):k] - Y_bar_prv)*(xt[(k-Gl+1):k] - X_bar_prv))
      Cov_YX_pst <- mean((yt[(k+1):(k+Gr)] - Y_bar_pst)*(xt[(k+1):(k+Gr)] - X_bar_pst))
      Cov_YX_bar[k] <- 1/2*(Cov_YX_pst + Cov_YX_prv)
      rho_YX[k] <- Cov_YX_bar[k]/sqrt(S2_Y_bar[k] * S2_X_bar[k])
      Gamma_YX <- diag(1,2); Gamma_YX[1,2] <- Gamma_YX[2,1] <- rho_YX[k]
      Dist_YX[k] <- sqrt(t(c(T1_Y[k],T1_X[k])) %*% solve(Gamma_YX) %*% c(T1_Y[k],T1_X[k]))
      pval_YX[k] <- 1-pchisq(Dist_YX[k]^2,df=2)
      
      # Y2X
      Cov_Y2X_prv <- mean( ((yt[(k-Gl+1):k] - Y_bar_prv)^2-S2_Y_prv)*(xt[(k-Gl+1):k] - X_bar_prv))
      Cov_Y2X_pst <- mean( ((yt[(k+1):(k+Gr)] - Y_bar_pst)^2-S2_Y_pst)*(xt[(k+1):(k+Gr)] - X_bar_pst))
      Cov_Y2X_bar[k] <- 1/2*(Cov_Y2X_pst + Cov_Y2X_prv)
      rho_Y2X[k] <- Cov_Y2X_bar[k]/sqrt(V2_Y_bar[k] * S2_X_bar[k])
      Gamma_Y2X <- diag(1,2); Gamma_Y2X[1,2] <- Gamma_Y2X[2,1] <- rho_Y2X[k]
      Dist_Y2X[k] <- sqrt(t(c(T2_Y[k],T1_X[k])) %*% solve(Gamma_Y2X) %*% c(T2_Y[k],T1_X[k]))
      pval_Y2X[k] <- 1-pchisq(Dist_Y2X[k]^2,df=2)
      
      # YX2
      Cov_YX2_prv <- mean((yt[(k-Gl+1):k] - Y_bar_prv)*((xt[(k-Gl+1):k] - X_bar_prv)^2-S2_X_prv) )
      Cov_YX2_pst <- mean((yt[(k+1):(k+Gr)] - Y_bar_pst)*((xt[(k+1):(k+Gr)] - X_bar_pst)^2-S2_X_pst) )
      Cov_YX2_bar[k] <- 1/2*(Cov_YX2_pst + Cov_YX2_prv)
      rho_YX2[k] <- Cov_YX2_bar[k]/sqrt(S2_Y_bar[k] * V2_X_bar[k])
      Gamma_YX2 <- diag(1,2); Gamma_YX2[1,2] <- Gamma_YX2[2,1] <- rho_YX2[k]
      Dist_YX2[k] <- sqrt(t(c(T1_Y[k],T2_X[k])) %*% solve(Gamma_YX2) %*% c(T1_Y[k],T2_X[k]))
      pval_YX2[k] <- 1-pchisq(Dist_YX2[k]^2,df=2)
      
      # Y2X2
      Cov_Y2X2_prv <- mean( ((yt[(k-Gl+1):k] - Y_bar_prv)^2-S2_Y_prv)*((xt[(k-Gl+1):k] - X_bar_prv)^2-S2_X_prv) )
      Cov_Y2X2_pst <- mean( ((yt[(k+1):(k+Gr)] - Y_bar_pst)^2-S2_Y_pst)*((xt[(k+1):(k+Gr)] - X_bar_pst)^2-S2_X_pst) )
      Cov_Y2X2_bar[k] <- 1/2*(Cov_Y2X2_pst + Cov_Y2X2_prv)
      rho_Y2X2[k] <- Cov_Y2X2_bar[k]/sqrt(V2_Y_bar[k] * V2_X_bar[k])
      Gamma_Y2X2 <- diag(1,2); Gamma_Y2X2[1,2] <- Gamma_Y2X2[2,1] <- rho_Y2X2[k]
      Dist_Y2X2[k] <- sqrt(t(c(T2_Y[k],T2_X[k])) %*% solve(Gamma_Y2X2) %*% c(T2_Y[k],T2_X[k]))
      pval_Y2X2[k] <- 1-pchisq(Dist_Y2X2[k]^2,df=2)
    }
  }
  
  # ---------------------------------------------------- #
  # Return:
  output <- list()
  output$Cov_YX_bar <- Cov_YX_bar
  output$Cov_YX2_bar <- Cov_YX2_bar
  output$Cov_Y2X_bar <- Cov_Y2X_bar
  output$Cov_Y2X2_bar <- Cov_Y2X2_bar
  
  output$rho_YX <- rho_YX
  output$rho_YX2 <- rho_YX2
  output$rho_Y2X <- rho_Y2X
  output$rho_Y2X2 <- rho_Y2X2
  
  output$Dist_YX <- Dist_YX
  output$Dist_YX2 <- Dist_YX2
  output$Dist_Y2X <- Dist_Y2X
  output$Dist_Y2X2 <- Dist_Y2X2
  
  output$pval_YX <- pval_YX
  output$pval_YX2 <- pval_YX2
  output$pval_Y2X <- pval_Y2X
  output$pval_Y2X2 <- pval_Y2X2
  return(output)
}

single_tracker <- function(series,n,Gl,Gr){
  
  # ------------------------------------------------------ #
  d_bar <- S2_bar <- Test1 <- pval_mu <- rep(NA,n)
  dS2 <- V2_bar <- Test2 <- pval_sigma2 <- rep(NA,n)
  K_bar <- rho <- Dist <- pval <- rep(NA,n)
  
  # ---------------------------------------------------- #
  for (k in 1:n){
    if (k <= Gl-1){
      
      # 1st order:
      series_bar <- mean(series[1:(Gl+Gr)])
      d_bar[k] <- 2/sqrt(k*(2*(Gl+Gr)-k))*sum(series_bar - series[1:k])
      S2_bar[k] <- mean( (series[1:(Gl+Gr)] - series_bar)^2 )
      Test1[k] <- sqrt(Gl*Gr/(Gl+Gr))*d_bar[k]/sqrt(S2_bar[k])
      pval_mu[k] <- 2*(1-pnorm(abs(Test1[k])))
      
      # 2nd order:
      series_bar2 <- mean( (series[1:(Gl+Gr)] - series_bar)^2 )
      dS2[k] <- 2/sqrt(k*(Gl+Gr-k))*sum(series_bar2 - (series[1:k]-series_bar)^2)
      V2_bar[k] <- mean(((series[1:(Gl+Gr)] - series_bar)^2 - S2_bar[k])^2)
      Test2[k] <- sqrt(Gl*Gr/(Gl+Gr))*dS2[k]/sqrt(V2_bar[k])
      pval_sigma2[k] <- 2*(1-pnorm(abs(Test2[k])))
      
      # Correlations:
      K_bar[k] <- mean( (series[1:(Gl+Gr)] - series_bar)^3 )
      rho[k] <- K_bar[k]/sqrt(S2_bar[k]*V2_bar[k])
      Gamma <- diag(1,2); Gamma[1,2] <- Gamma[2,1] <- rho[k]
      Dist[k] <- sqrt(t(c(Test1[k],Test2[k])) %*% solve(Gamma) %*% c(Test1[k],Test2[k]))
      pval[k] <- 1-pchisq(Dist[k]^2,df=2)
      
    }else if (k >= n-Gr+1){
      
      # 1st order:
      series_bar <- mean(series[(n-(Gl+Gr)+1):n])
      d_bar[k] <- 2/sqrt((n+1-k)*(k-n+Gl+Gr))*sum(series_bar - series[k:n])
      S2_bar[k] <- mean( (series[(n-(Gl+Gr)+1):n] - series_bar)^2 )
      Test1[k] <- sqrt(Gl*Gr/(Gl+Gr))*d_bar[k]/sqrt(S2_bar[k])
      pval_mu[k] <- 2*(1-pnorm(abs(Test1[k])))
      
      # 2nd order:
      series_bar2 <- mean( (series[(n-(Gl+Gr)+1):n] - series_bar)^2 )
      dS2[k] <- 2/sqrt((n+1-k)*(k-n+(Gl+Gr)))*sum(series_bar2 - (series[k:n]-series_bar)^2)
      V2_bar[k] <- mean(((series[(n-(Gl+Gr)+1):n] - series_bar)^2 - S2_bar[k])^2)
      Test2[k] <- sqrt(Gl*Gr/(Gl+Gr))*dS2[k]/sqrt(V2_bar[k])
      pval_sigma2[k] <- 2*(1-pnorm(abs(Test2[k])))
      
      # Correlations:
      K_bar[k] <- mean( (series[(n-(Gl+Gr)+1):n] - series_bar)^3 )
      rho[k] <- K_bar[k]/sqrt(S2_bar[k]*V2_bar[k])
      Gamma <- diag(1,2); Gamma[1,2] <- Gamma[2,1] <- rho[k]
      Dist[k] <- sqrt(t(c(Test1[k],Test2[k])) %*% solve(Gamma) %*% c(Test1[k],Test2[k]))
      pval[k] <- 1-pchisq(Dist[k]^2,df=2)
      
    }else{
      
      # 1st order:
      series_bar_prv <- mean(series[(k-Gl+1):k])
      series_bar_pst <- mean(series[(k+1):(k+Gr)])
      d_bar[k] <- series_bar_pst - series_bar_prv
      
      series_bar2_prv <- mean( (series[(k-Gl+1):k] - series_bar_prv)^2 )
      series_bar2_pst <- mean( (series[(k+1):(k+Gr)] - series_bar_pst)^2 )
      
      S2_bar[k] <- 1/2*(series_bar2_pst + series_bar2_prv)
      Test1[k] <- sqrt(Gl*Gr/(Gl+Gr))*d_bar[k]/sqrt(S2_bar[k])
      pval_mu[k] <-  2*(1-pnorm(abs(Test1[k])))
      
      # 2nd order:
      dS2[k] <- series_bar2_pst - series_bar2_prv
      series_bar4_prv <- mean( ((series[(k-Gl+1):k] - series_bar_prv)^2 - series_bar2_prv)^2 )
      series_bar4_pst <- mean( ((series[(k+1):(k+Gr)] - series_bar_pst)^2 - series_bar2_pst)^2 )
      
      V2_bar[k] <- 1/2*(series_bar4_pst + series_bar4_prv)
      Test2[k] <- sqrt(Gl*Gr/(Gl+Gr))*dS2[k]/sqrt(V2_bar[k])
      pval_sigma2[k] <-  2*(1-pnorm(abs(Test2[k])))
      
      # Correlations:
      series_bar3_prv <- mean((series[(k-Gl+1):k] - series_bar_prv)^3)
      series_bar3_pst <- mean((series[(k+1):(k+Gr)] - series_bar_pst)^3)
      K_bar[k] <- 1/2*(series_bar3_pst + series_bar3_prv)
      rho[k] <- K_bar[k]/sqrt(S2_bar[k] * V2_bar[k])
      Gamma <- diag(1,2); Gamma[1,2] <- Gamma[2,1] <- rho[k]
      Dist[k] <- sqrt(t(c(Test1[k],Test2[k])) %*% solve(Gamma) %*% c(Test1[k],Test2[k]))
      pval[k] <- 1-pchisq(Dist[k]^2,df=2)
    }
  }
  
  # ---------------------------------------------------- #
  # Return:
  output <- list()
  output$series <- series
  output$d_bar <- d_bar
  output$S2_bar <- S2_bar
  output$dS2 <- dS2
  output$V2_bar <- V2_bar
  output$K_bar <- K_bar
  output$rho <- rho
  
  output$Test1 <- Test1
  output$Test2 <- Test2
  output$Dist <- Dist
  
  output$pval_mu <- pval_mu
  output$pval_sigma2 <- pval_sigma2
  output$pval <- pval
  return(output)
}

############################################
# refinement function for MOSUM/Bi-MOSUM
############################################
refinements <- function(test_stat,type,n,Gl,Gr,alpha,epsilon,eta){
  
  # ------------------------------------------------------ #
  # Initial candidates:
  if (type == "univ"){
    G <- mean(c(Gl,Gr))
    an <- sqrt(2*log(n/G)) 
    bn <- 2*log(n/G) + 0.5*log(log(n/G)) + log(3/2) - 0.5*log(pi)
    c_alpha <- -log(log(1/sqrt(1-alpha)))
    Dn <- (bn + c_alpha)/an
    
  }else{
    # ------------------------------------------------------ #
    # Threshold:
    if (Gl != Gr){
      stop("Monte Carlo threshold calibration currently assumes Gl = Gr.")
    }
    
    G <- Gl
    
    max_G <- function(g, WE, WV) {
      LE <- (WE[(2*g+1):n] - 2*WE[(1+g):(n-g)] + WE[1:(n-2*g)])/sqrt(2*g)
      LV <- (WV[(2*g+1):n] - 2*WV[(1+g):(n-g)] + WV[1:(n-2*g)])/sqrt(2*g)
      max(sqrt(LE^2 + LV^2))
    }
    sim_global_max <- function(n, G) {
      WE <- cumsum(rnorm(n))
      WV <- cumsum(rnorm(n))
      max_G(G, WE, WV)
    }
    Dn <- quantile(replicate(1000, sim_global_max(n=n, G=G)),1-alpha)
  }
  idx_cand <- which(abs(test_stat) >= Dn)
  
  # ---------------------------------------------------- #
  # Epsilon-criterion: significant neighborhoods
  if (length(idx_cand)!=0){
    # Get runs:
    grp <- cumsum(c(TRUE, diff(idx_cand) != 1))
    runs <- split(idx_cand, grp)
    
    cpt_eps <- list()
    for (i in 1:length(runs)){
      idx <- runs[[i]]
      l <- min(idx)-1
      r <- max(idx)+1
      if (r-l >= epsilon*G){
        max_val <- max(abs(test_stat[idx]))
        cpt_eps[[i]] <- which(test_stat == max_val)
      }else{
        next
      }
    }
    cpt_eps <- unlist(cpt_eps) 
  }else{
    cpt_eps <- NULL
  }
  
  # ---------------------------------------------------- #
  # Eta-criterion: narrower search
  if (length(idx_cand)!=0){
    band_l <- floor(eta*Gl)
    band_r <- floor(eta*Gr)
    cpt_eta <- list()
    for (i in 1:length(idx_cand)){
      idx <- idx_cand[i]
      if ( max(idx-band_l)<1 || min(idx+band_r)>n ){
        next
        
      }else{
        idx_band <- (idx-band_l):(idx+band_r)
        idx_band <- idx_band[which(idx_band>=1)]
        idx_band <- idx_band[which(idx_band<=n)]
        
        if (max(abs(test_stat[idx_band])) == abs(test_stat[idx])){
          cpt_eta[[i]] <- idx
        }else{
          next
        }
      }
    }
    cpt_eta <- unlist(cpt_eta)
  }else{
    cpt_eta <- NULL
  }
  
  # ---------------------------------------------------- #
  # Return:
  output <- list()
  output$Dn <- Dn
  output$idx_cand <- idx_cand
  output$cpt_eta <- cpt_eta
  output$cpt_eps <- cpt_eps
  return(output)
}


############################################
# Confidence interval for Bi-MOSUM
############################################
bivariate_mosum_ci <- function(series_Y,series_X,bi_refinement,
                               moment=c("mean_mean","mean_var","var_mean","var_var"),
                               method=c("epsilon","eta"),n,Gl=30,Gr=30,N=1000,
                               alpha=0.05,epsilon=0.2,eta=0.4){
  
  # ---------------------------------------------------- #
  # Set up:
  if (method == "epsilon"){
    cps <- bi_refinement$cpt_eps
    if (length(cps)==0){
      stop("No change point detected using epsilon method") 
    }
    
  }else{ # method == "eta"
    cps <- bi_refinement$cpt_eta
    if (length(cps)==0){
      stop("No change point detected using eta method") 
    }
  }
  
  time_idx <- 1:n
  k_index <- c(0,cps,n)
  
  # ---------------------------------------------------- #
  # Bootstrapping:
  series_X_rep <- series_Y_rep <- list()
  single_X_rep <- single_Y_rep <- list()
  initial_rep <- list()
  cps_rep <- list()
  for (iter in 1:N){
    samples <- c()
    for (k in 1:(length(k_index)-1)){
      inteval_k <- (k_index[k]+1):k_index[k+1]
      length_k <- length(inteval_k)
      samples <- c(samples,sample(inteval_k,length_k,replace=TRUE))
    }
    series_X_rep[[iter]] <- series_X[samples]
    series_Y_rep[[iter]] <- series_Y[samples]
    single_X_rep[[iter]] <- single_tracker(series_X_rep[[iter]],n,Gl,Gr)
    single_Y_rep[[iter]] <- single_tracker(series_Y_rep[[iter]],n,Gl,Gr)
    
    initial_rep[[iter]] <- bivariate_tracker(single_Y_rep[[iter]],single_X_rep[[iter]],n,Gl,Gr)
    
    G <- min(Gl,Gr)
    sample_cps <- rep(NA,length(cps))
    for (k in 2:(length(k_index)-1)){
      nbh_k <- max((k_index[k]-G+1),1):min((k_index[k]+G),n)
      if (moment == "mean_mean"){
        sample_cps[k-1] <- which.max(abs(initial_rep[[iter]]$Dist_YX[nbh_k])) + (k_index[k]-G)
      }else if (moment == "mean_var"){
        sample_cps[k-1] <- which.max(abs(initial_rep[[iter]]$Dist_YX2[nbh_k])) + (k_index[k]-G)
      }else if (moment == "var_mean"){
        sample_cps[k-1] <- which.max(abs(initial_rep[[iter]]$Dist_Y2X[nbh_k])) + (k_index[k]-G)
      }else{
        sample_cps[k-1] <- which.max(abs(initial_rep[[iter]]$Dist_Y2X2[nbh_k])) + (k_index[k]-G)
      }
    }
    cps_rep[[iter]] <- sample_cps
  }
  
  # ---------------------------------------------------- #
  # Pointwise C.I:
  M_j_alpha <- rep(NA,length(cps))
  K_pw <- list()
  for (j in 1:length(cps)){
    rep_j <- abs(cps[j] - mapply(iter=1:N,function(iter)cps_rep[[iter]][j]))
    M_j_alpha[j] <- quantile(rep_j,(1-alpha/2))
    K_pw[[j]] <- c( (cps[j]-M_j_alpha[j]),(cps[j]+M_j_alpha[j]) )
  }
  
  # ---------------------------------------------------- #
  # Uniform C.I:
  K_unif <- list()
  scales <- matrix(NA,N,length(cps))
  coef_hat <- rep(NA,length(cps))
  for (j in 1:length(cps)){
    prev_j <- (k_index[j]+1):k_index[j+1]
    post_j <- (k_index[j+1]+1):k_index[j+2]
    
    Y_bar_prv <- mean(series_Y[prev_j])
    Y_bar_pst <- mean(series_Y[post_j])
    X_bar_prv <- mean(series_X[prev_j])
    X_bar_pst <- mean(series_X[post_j])
    
    S2_Y_prv <- mean( (series_Y[prev_j] - Y_bar_prv)^2 )
    S2_Y_pst <- mean( (series_Y[post_j] - Y_bar_pst)^2 )
    S2_X_prv <- mean( (series_X[prev_j] - X_bar_prv)^2 )
    S2_X_pst <- mean( (series_X[post_j] - X_bar_pst)^2 )
    
    S2_Y_bar <- 1/2*(S2_Y_prv + S2_Y_pst)
    S2_X_bar <- 1/2*(S2_X_prv + S2_X_pst)
    
    series_Y_bar4_prv <- mean( ((series_Y[prev_j] - S2_Y_prv)^2 - S2_Y_prv)^2 )
    series_Y_bar4_pst <- mean( ((series_Y[post_j] - S2_Y_pst)^2 - S2_Y_pst)^2 )
    series_X_bar4_prv <- mean( ((series_X[prev_j] - S2_X_prv)^2 - S2_X_prv)^2 )
    series_X_bar4_pst <- mean( ((series_X[post_j] - S2_X_pst)^2 - S2_X_pst)^2 )
    
    V2_Y_bar <- 1/2*(series_Y_bar4_prv + series_Y_bar4_pst)
    V2_X_bar <- 1/2*(series_X_bar4_prv + series_X_bar4_pst)
    
    if (moment == "mean_mean"){
      Cov_YX_prv <- mean((series_Y[prev_j] - Y_bar_prv)*(series_X[prev_j] - X_bar_prv))
      Cov_YX_pst <- mean((series_Y[post_j] - Y_bar_pst)*(series_X[post_j] - X_bar_pst))
      Cov_YX_bar <- 1/2*(Cov_YX_pst + Cov_YX_prv)
      rho_YX <- Cov_YX_bar/sqrt(S2_Y_bar * S2_X_bar)
      Gamma_YX <- diag(1,2); Gamma_YX[1,2] <- Gamma_YX[2,1] <- rho_YX
      
      T1_Y <- (Y_bar_pst - Y_bar_prv)/sqrt(S2_Y_bar)
      T1_X <- (X_bar_pst - X_bar_prv)/sqrt(S2_X_bar)
      coef_hat[j] <- t(c(T1_Y,T1_X)) %*% solve(Gamma_YX) %*% c(T1_Y,T1_X)
      
    }else if (moment == "mean_var"){
      Cov_YX2_prv <- mean((series_Y[prev_j] - Y_bar_prv)*((series_X[prev_j] - X_bar_prv)^2-S2_X_prv) )
      Cov_YX2_pst <- mean((series_Y[post_j] - Y_bar_pst)*((series_X[post_j] - X_bar_pst)^2-S2_X_pst) )
      Cov_YX2_bar <- 1/2*(Cov_YX2_pst + Cov_YX2_prv)
      rho_YX2 <- Cov_YX2_bar/sqrt(S2_Y_bar * V2_X_bar)
      Gamma_YX2 <- diag(1,2); Gamma_YX2[1,2] <- Gamma_YX2[2,1] <- rho_YX2
      
      T1_Y <- (Y_bar_pst - Y_bar_prv)/sqrt(S2_Y_bar)
      T2_X <- (S2_X_pst - S2_X_prv)/sqrt(V2_X_bar)
      coef_hat[j] <- t(c(T1_Y,T2_X)) %*% solve(Gamma_YX2) %*% c(T1_Y,T2_X)
      
    }else if (moment == "var_mean"){
      Cov_Y2X_prv <- mean( ((series_Y[prev_j] - Y_bar_prv)^2-S2_Y_prv)*(series_X[prev_j] - X_bar_prv))
      Cov_Y2X_pst <- mean( ((series_Y[post_j] - Y_bar_pst)^2-S2_Y_pst)*(series_X[post_j] - X_bar_pst))
      Cov_Y2X_bar <- 1/2*(Cov_Y2X_pst + Cov_Y2X_prv)
      rho_Y2X <- Cov_Y2X_bar/sqrt(V2_Y_bar * S2_X_bar)
      Gamma_Y2X <- diag(1,2); Gamma_Y2X[1,2] <- Gamma_Y2X[2,1] <- rho_Y2X
      
      T2_Y <- (S2_Y_pst - S2_Y_prv)/sqrt(V2_Y_bar)
      T1_X <- (X_bar_pst - X_bar_prv)/sqrt(S2_X_bar)
      coef_hat[j] <- t(c(T2_Y,T1_X)) %*% solve(Gamma_Y2X) %*% c(T2_Y,T1_X)
      
    }else if (moment == "var_var"){
      Cov_Y2X2_prv <- mean( ((series_Y[prev_j] - Y_bar_prv)^2-S2_Y_prv)*((series_X[prev_j] - X_bar_prv)^2-S2_X_prv) )
      Cov_Y2X2_pst <- mean( ((series_Y[post_j] - Y_bar_pst)^2-S2_Y_pst)*((series_X[post_j] - X_bar_pst)^2-S2_X_pst) )
      Cov_Y2X2_bar <- 1/2*(Cov_Y2X2_pst + Cov_Y2X2_prv)
      rho_Y2X2 <- Cov_Y2X2_bar/sqrt(V2_Y_bar * V2_X_bar)
      Gamma_Y2X2 <- diag(1,2); Gamma_Y2X2[1,2] <- Gamma_Y2X2[2,1] <- rho_Y2X2
      
      T2_Y <- (S2_Y_pst - S2_Y_prv)/sqrt(V2_Y_bar)
      T2_X <- (S2_X_pst - S2_X_prv)/sqrt(V2_X_bar)
      coef_hat[j] <- t(c(T2_Y,T2_X)) %*% solve(Gamma_Y2X2) %*% c(T2_Y,T2_X)
    }
    
    for (iter in 1:N){
      scales[iter,j] <- coef_hat[j]*(cps_rep[[iter]][j] - cps[j])
    }
  }
  
  const <- 0; 
  max_1m_alpha <- max(mapply(j=1:length(cps),function(j)sum(scales[,j]<=const)/N))
  if (max_1m_alpha >= 1-alpha){
    M_alpha <- const
  }else{
    while (max_1m_alpha < 1-alpha){
      const <- const + 0.1
      max_1m_alpha <- max(mapply(j=1:length(cps),function(j)sum(scales[,j]<=const)/N))
    }
    M_alpha <- const
  }
  
  for (j in 1:length(cps)){
    K_unif[[j]] <- c( (cps[j]-M_alpha/coef_hat[j]),(cps[j]+M_alpha/coef_hat[j]) )
  }
  
  # ---------------------------------------------------- #
  output <- list()
  output$ptwise <- K_pw
  output$uniform <- K_unif
  return(output)
}

############################################
# Confidence interval for Joint-MOSUM
############################################
joint_mosum_ci <- function(series,refinement,
                           method=c("epsilon","eta"),n,Gl=30,Gr=30,N=1000,
                           alpha=0.05,epsilon=0.2,eta=0.4){
  
  # ---------------------------------------------------- #
  # Set up:
  if (method == "epsilon"){
    cps <- refinement$cpt_eps
    if (length(cps)==0){
      stop("No change point detected using epsilon method") 
    }
    
  }else{ # method == "eta"
    cps <- refinement$cpt_eta
    if (length(cps)==0){
      stop("No change point detected using eta method") 
    }
  }
  
  time_idx <- 1:n
  k_index <- c(0,cps,n)
  
  # ---------------------------------------------------- #
  # Bootstrapping:
  series_rep <- single_rep <- list()
  initial_rep <- list()
  cps_rep <- list()
  for (iter in 1:N){
    samples <- c()
    for (k in 1:(length(k_index)-1)){
      inteval_k <- (k_index[k]+1):k_index[k+1]
      length_k <- length(inteval_k)
      samples <- c(samples,sample(inteval_k,length_k,replace=TRUE))
    }
    series_rep[[iter]] <- series[samples]
    single_rep[[iter]] <- single_tracker(series_rep[[iter]],n,Gl,Gr)
    
    G <- min(Gl,Gr)
    sample_cps <- rep(NA,length(cps))
    for (k in 2:(length(k_index)-1)){
      nbh_k <- max((k_index[k]-G+1),1):min((k_index[k]+G),n)
      sample_cps[k-1] <- which.max(abs(single_rep[[iter]]$Dist[nbh_k])) + (k_index[k]-G)
    }
    cps_rep[[iter]] <- sample_cps
  }
  
  # ---------------------------------------------------- #
  # Pointwise C.I:
  M_j_alpha <- rep(NA,length(cps))
  K_pw <- list()
  for (j in 1:length(cps)){
    rep_j <- abs(cps[j] - mapply(iter=1:N,function(iter)cps_rep[[iter]][j]))
    M_j_alpha[j] <- quantile(rep_j,(1-alpha/2))
    K_pw[[j]] <- c( (cps[j]-M_j_alpha[j]),(cps[j]+M_j_alpha[j]) )
  }
  
  # ---------------------------------------------------- #
  # Uniform C.I:
  K_unif <- list()
  scales <- matrix(NA,N,length(cps))
  coef_hat <- rep(NA,length(cps))
  for (j in 1:length(cps)){
    prev_j <- (k_index[j]+1):k_index[j+1]
    post_j <- (k_index[j+1]+1):k_index[j+2]
    
    series_bar_prv <- mean(series[prev_j])
    series_bar_pst <- mean(series[post_j])
    
    S2_prv <- mean( (series[prev_j] - series_bar_prv)^2 )
    S2_pst <- mean( (series[post_j] - series_bar_pst)^2 )
    
    S2_bar <- 1/2*(S2_prv + S2_pst)
    
    series_bar4_prv <- mean(((series[prev_j] - series_bar_prv)^2 - S2_prv)^2)
    series_bar4_pst <- mean(((series[post_j] - series_bar_pst)^2 - S2_pst)^2)
    
    V2_bar <- 1/2*(series_bar4_prv + series_bar4_pst)
    
    series_bar3_prv <- mean((series[prev_j] - series_bar_prv)^3)
    series_bar3_pst <- mean((series[post_j] - series_bar_pst)^3)
    K_bar <- 1/2*(series_bar3_pst + series_bar3_prv)
    rho <- K_bar/sqrt(S2_bar * V2_bar)
    Gamma <- diag(1,2); Gamma[1,2] <- Gamma[2,1] <- rho
    
    T1 <- (series_bar_pst - series_bar_prv)/sqrt(S2_bar)
    T2 <- (S2_pst - S2_prv)/sqrt(V2_bar)
    coef_hat[j] <- t(c(T1,T2)) %*% solve(Gamma) %*% c(T1,T2)
    
    for (iter in 1:N){
      scales[iter,j] <- coef_hat[j]*(cps_rep[[iter]][j] - cps[j])
    }
  }
  
  const <- 0; 
  max_1m_alpha <- max(mapply(j=1:length(cps),function(j)sum(scales[,j]<=const)/N))
  if (max_1m_alpha >= 1-alpha){
    M_alpha <- const
  }else{
    while (max_1m_alpha < 1-alpha){
      const <- const + 0.1
      max_1m_alpha <- max(mapply(j=1:length(cps),function(j)sum(scales[,j]<=const)/N))
    }
    M_alpha <- const
  }
  
  for (j in 1:length(cps)){
    K_unif[[j]] <- c( (cps[j]-M_alpha/coef_hat[j]),(cps[j]+M_alpha/coef_hat[j]) )
  }
  
  # ---------------------------------------------------- #
  output <- list()
  output$ptwise <- K_pw
  output$uniform <- K_unif
  return(output)
}

############################################
# Confidence interval for MOSUM
############################################
univariate_mosum_ci <- function(series,single_refinement,moment=c("mean","var"),
                                method=c("epsilon","eta"),n,Gl=30,Gr=30,N=1000,
                                alpha=0.05,epsilon=0.2,eta=0.4){
  
  # ---------------------------------------------------- #
  # Set up:
  if (method == "epsilon"){
    cps <- single_refinement$cpt_eps
    if (length(cps)==0){
      stop("No change point detected using epsilon method") 
    }
    
  }else{ # method == "eta"
    cps <- single_refinement$cpt_eta
    if (length(cps)==0){
      stop("No change point detected using eta method") 
    }
  }
  time_idx <- 1:n
  k_index <- c(0,cps,n)
  
  # ---------------------------------------------------- #
  # Bootstrapping:
  series_rep <- list()
  initial_rep <- list()
  cps_rep <- list()
  for (iter in 1:N){
    samples <- c()
    for (k in 1:(length(k_index)-1)){
      inteval_k <- (k_index[k]+1):k_index[k+1]
      length_k <- length(inteval_k)
      samples <- c(samples,sample(inteval_k,length_k,replace=TRUE))
    }
    series_rep[[iter]] <- series[samples]
    initial_rep[[iter]] <- single_tracker(series_rep[[iter]],n,Gl,Gr)
    
    G <- min(Gl,Gr)
    sample_cps <- rep(NA,length(cps))
    for (k in 2:(length(k_index)-1)){
      nbh_k <- max((k_index[k]-G+1),1):min((k_index[k]+G),n)
      if (moment == "mean"){
        sample_cps[k-1] <- which.max(abs(initial_rep[[iter]]$Test1[nbh_k])) + (k_index[k]-G)
      }else{
        sample_cps[k-1] <- which.max(abs(initial_rep[[iter]]$Test2[nbh_k])) + (k_index[k]-G)
      }
    }
    cps_rep[[iter]] <- sample_cps
  }
  
  # ---------------------------------------------------- #
  # Pointwise C.I:
  M_j_alpha <- rep(NA,length(cps))
  K_pw <- list()
  for (j in 1:length(cps)){
    rep_j <- abs(cps[j] - mapply(iter=1:N,function(iter)cps_rep[[iter]][j]))
    M_j_alpha[j] <- quantile(rep_j,(1-alpha/2))
    K_pw[[j]] <- c( (cps[j]-M_j_alpha[j]),(cps[j]+M_j_alpha[j]) )
  }
  
  # ---------------------------------------------------- #
  # Uniform C.I:
  K_unif <- list()
  scales <- matrix(NA,N,length(cps))
  coef_hat <- rep(NA,length(cps))
  for (j in 1:length(cps)){
    prev_j <- (k_index[j]+1):k_index[j+1]
    post_j <- (k_index[j+1]+1):k_index[j+2]
    mean_prev_j <- mean(series[prev_j])
    mean_post_j <- mean(series[post_j])
    var_prev_j <- mean(series[prev_j]^2)
    var_post_j <- mean(series[post_j]^2)
    
    if (moment == "mean"){
      dbar_j_hat2 <- (mean_post_j - mean_prev_j)^2
      S2_j_hat2 <- (sum((series[prev_j] - mean_prev_j)^2) 
                    + sum((series[post_j] - mean_post_j)^2)) / (k_index[j+2] - k_index[j] - 1)
      coef_hat[j] <- dbar_j_hat2/S2_j_hat2
      
    }else if (moment == "var"){ 
      dS2_j_hat2 <- (var_post_j - var_prev_j)^2
      V2_j_hat2 <- (sum( ((series[prev_j] - mean_prev_j)^2 - var_prev_j)^2 ) 
                    + sum( ((series[post_j] - mean_post_j)^2 - var_post_j)^2 ) ) / (k_index[j+2] - k_index[j] - 1)
      coef_hat[j] <- dS2_j_hat2/V2_j_hat2
      
    }
    
    for (iter in 1:N){
      scales[iter,j] <- coef_hat[j]*(cps_rep[[iter]][j] - cps[j])
    }
  }
  
  const <- 0; 
  max_1m_alpha <- max(mapply(j=1:length(cps),function(j)sum(scales[,j]<=const)/N))
  if (max_1m_alpha >= 1-alpha){
    M_alpha <- const
  }else{
    while (max_1m_alpha < 1-alpha){
      const <- const + 0.1
      max_1m_alpha <- max(mapply(j=1:length(cps),function(j)sum(scales[,j]<=const)/N))
    }
    M_alpha <- const
  }
  
  for (j in 1:length(cps)){
    K_unif[[j]] <- c( (cps[j]-M_alpha/coef_hat[j]),(cps[j]+M_alpha/coef_hat[j]) )
  }
  
  # ---------------------------------------------------- #
  output <- list()
  output$ptwise <- K_pw
  output$uniform <- K_unif
  return(output)
}


############################################
# convert ordinal responses to continuous variables
############################################
conv_smooth <- function(series,n,scale){
  
  prob_est <- table(factor(series, levels = scale))/(n)
  bin <- c(0,cumsum(prob_est))
  Ut <- rep(NA,n)
  for (k in 1:n){
    Ut[k] <- bin[series[k]] + runif(1,0,1)*prob_est[series[k]] 
  }
  return(qnorm(Ut))
}

############################################
# List intervals merger:
############################################
get_union_of_intersects <- function(a, b, type = c("uniform", "ptwise"), tol = 0) {
  type <- match.arg(type)
  
  # Helper: extract the interval list (works whether a/b are full objects or sublists)
  extract_intervals <- function(x) {
    if (is.list(x) && !is.null(x[[type]])) {
      intervals <- x[[type]]
    } else {
      intervals <- x
    }
    if (is.numeric(intervals) && length(intervals) == 2) intervals <- list(intervals)
    if (is.null(intervals)) intervals <- list()
    intervals
  }
  
  intervals1 <- extract_intervals(a)
  intervals2 <- extract_intervals(b)
  
  # Intersection between two intervals
  interval_intersect <- function(x, y) {
    low  <- max(x[1], y[1])
    high <- min(x[2], y[2])
    if (low <= high + tol) return(c(low, high))
    NULL
  }
  
  # Get all intersections
  intersections <- list()
  for (x in intervals1) {
    for (y in intervals2) {
      inter <- interval_intersect(x, y)
      if (!is.null(inter)) intersections <- append(intersections, list(inter))
    }
  }
  
  # If no intersections found
  if (length(intersections) == 0) return(list())
  
  # Merge overlapping intervals
  merge_intervals <- function(intervals) {
    if (length(intervals) <= 1) return(intervals)
    starts <- sapply(intervals, `[`, 1)
    intervals <- intervals[order(starts)]
    merged <- list(intervals[[1]])
    for (i in seq(2, length(intervals))) {
      last <- merged[[length(merged)]]
      cur  <- intervals[[i]]
      if (cur[1] <= last[2] + tol) {
        merged[[length(merged)]] <- c(last[1], max(last[2], cur[2]))
      } else {
        merged <- append(merged, list(cur))
      }
    }
    merged
  }
  
  merge_intervals(intersections)
}


############################################
# List merger (for two names):
############################################
merge_lists <- function(x, y) {
  all_names <- union(names(x), names(y))
  merged <- list()
  
  for (nm in all_names) {
    x_part <- if (!is.null(x[[nm]])) x[[nm]] else list()
    y_part <- if (!is.null(y[[nm]])) y[[nm]] else list()
    merged[[nm]] <- c(x_part, y_part)
  }
  
  merged
}

############################################
# Get intervals:
############################################
get_intervals <- function(x) {
  x <- sort(unique(x))  # ensure sorted and unique values
  
  # Find breaks between non-consecutive numbers
  breaks <- c(0, which(diff(x) > 1), length(x))
  
  # Build and return list of (left, right) pairs
  lapply(seq_along(breaks[-1]), function(i) {
    c(x[breaks[i] + 1], x[breaks[i + 1]])
  })
}