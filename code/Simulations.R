##################################################################
# Numerical Experiment 1 (univariate):
##################################################################
get_sim <- function(i,j,iter_max){
  
  library(doParallel)
  require(dqrng)
  require(doRNG)
  
  cl <- makeCluster(50)
  registerDoParallel(cl)
  
  n <- 100
  if (i == 1){
    num_jump <- 1
    start <- c(1,51)
    end <- c(50,n)
  }else if (i == 2){
    num_jump <- 2
    start <- c(1,41,61)
    end <- c(40,60,n)
  }else if (i == 3){
    num_jump <- 3
    start <- c(1,26,51,76)
    end <- c(25,50,75,n)
  }
  
  if (j == 1){
    level_range <- 2
    sigma_min <- 0.1
    sigma_max <- 0.4
  }else if (j == 2){
    level_range <- 2
    sigma_min <- 0.1
    sigma_max <- 0.8
  }else if (j == 3){
    level_range <- 2
    sigma_min <- 0.4
    sigma_max <- 0.8
  }else if (j == 4){
    level_range <- 1
    sigma_min <- 0.1
    sigma_max <- 0.4
  }else if (j == 5){
    level_range <- 1
    sigma_min <- 0.1
    sigma_max <- 0.8
  }else if (j == 6){
    level_range <- 1
    sigma_min <- 0.4
    sigma_max <- 0.8
  }
  
  set.seed(123)
  list_iter <- foreach(iter = 1:iter_max,
                       .errorhandling = 'remove',
                       .packages = c("mvtnorm","stringr",
                                     "breakfast","changepoint","stepR","mosum")) %dopar% {
                                       
                                       source("Library.R")
                                       
                                       # ------------------------------------------- #
                                       n <- 100
                                       series <- series_noisy <- rep(NA,n)
                                       levels_vec <- runif((num_jump+1), min=-level_range, max=level_range)
                                       sigma_vec <- runif((num_jump+1), min=sigma_min, max=sigma_max)
                                       for (k in 1:(num_jump+1)) {
                                         series[start[k]:end[k]] <- levels_vec[k]
                                         series_noisy[start[k]:end[k]] <- series[start[k]:end[k]] + rnorm(end[k]-start[k]+1, mean=0, sd=sigma_vec[k])
                                       }
                                       
                                       # ------------------------------------------- #
                                       # Joint-MOSUM:
                                       G20 <- single_tracker(series_noisy,n,Gl=20,Gr=20)
                                       G20_both <- refinements(G20$Dist,type="both",n,
                                                               Gl=20,Gr=20,alpha=0.05,epsilon=0.2,eta=0.2)
                                       G20_both <- G20_both$cpt_eta
                                       
                                       G40 <- single_tracker(series_noisy,n,Gl=30,Gr=30)
                                       G40_both <- refinements(G40$Dist,type="both",n,
                                                               Gl=40,Gr=40,alpha=0.05,epsilon=0.2,eta=0.2)
                                       G40_both <- G40_both$cpt_eta
                                       
                                       
                                       # ------------------------------------------- #
                                       # PELT:
                                       PELT <- cpt.meanvar(series_noisy,penalty="MBIC",method="AMOC",
                                                           Q=5,test.stat="Normal",class=TRUE,
                                                           param.estimates=TRUE)
                                       PELT <- cpts(PELT)
                                       
                                       # ------------------------------------------- #
                                       # stepR:
                                       stepR <- stepFit(series_noisy,alpha=0.05,jumpint=TRUE,confband=TRUE)
                                       stepR <- setdiff(stepR$leftEnd,c(1,n))
                                       
                                       # ------------------------------------------- #
                                       # WBS:
                                       WBS <- breakfast(series_noisy,solution.path="wbs")
                                       WBS <- WBS$cptmodel.list[[1]]$cpts
                                       
                                       # ------------------------------------------- #
                                       # MOSUM:
                                       G20_bm <- mosum(series_noisy,G=20)$cpts.info$cpts
                                       G40_bm <- mosum(series_noisy,G=40)$cpts.info$cpts
                                       
                                       # --------------------------------------------- #
                                       output <- list()
                                       output$G20_both <- G20_both
                                       output$G40_both <- G40_both
                                       output$PELT <- PELT
                                       output$stepR <- stepR
                                       output$WBS <- WBS
                                       output$G20_bm <- G20_bm
                                       output$G40_bm <- G40_bm
                                       output
                                     }
  stopCluster(cl)
  print(length(list_iter))
  if (length(list_iter) == 0) {
    warning("No successful iterations were returned.")
  }
  save(list_iter,file=paste0("univariate_jump",i,"_snr",j,".RData"))
}

# --------------------------------------------------- #
# Run iterations
for (i in 1:3){
  for (j in 1:6){
    cat("The current iteration is i:",i,"j:",j,"\n")
    get_sim(i,j,iter_max=500)
    cat("Done.\n") 
  }
}

# --------------------------------------------------- #
# Results:
max_iter <- 500
for (i in 1:3){
  for (j in 1:6){
    
    load(paste0("univariate_jump",i,"_snr",j,".RData"))
    eta <- 5
    
    # --------------------------------------------------- #
    G20_both <- mapply(iter=1:max_iter,function(iter)list_iter[[iter]]$G20_both)
    G40_both <- mapply(iter=1:max_iter,function(iter)list_iter[[iter]]$G40_both)
    PELT <- mapply(iter=1:max_iter,function(iter)list_iter[[iter]]$PELT)
    stepR <- mapply(iter=1:max_iter,function(iter)list_iter[[iter]]$stepR)
    wbs <-mapply(iter=1:max_iter,function(iter)list_iter[[iter]]$WBS)
    G20_bm <- mapply(iter=1:max_iter,function(iter)list_iter[[iter]]$G20_bm)
    G40_bm <- mapply(iter=1:max_iter,function(iter)list_iter[[iter]]$G40_bm)
    
    if (i == 1){
      # --------------------------------------------------- #
      # Power:
      power_G20_both <- sum(sapply(G20_both, function(vec) {
        in_50 <- any(abs(vec-50) <= eta); 
        in_50 }))/max_iter
      power_G40_both <- sum(sapply(G40_both, function(vec) {
        in_50 <- any(abs(vec-50) <= eta); 
        in_50 }))/max_iter
      power_PELT <- sum(sapply(PELT, function(vec) {
        in_50 <- any(abs(vec-50) <= eta); 
        in_50 }))/max_iter
      power_stepR <- sum(sapply(stepR, function(vec) {
        in_50 <- any(abs(vec-50) <= eta); 
        in_50 }))/max_iter
      power_wbs <- sum(sapply(wbs, function(vec) {
        in_50 <- any(abs(vec-50) <= eta); 
        in_50 }))/max_iter
      power_G20_bm <- sum(sapply(G20_bm, function(vec) {
        in_50 <- any(abs(vec-50) <= eta); 
        in_50 }))/max_iter
      power_G40_bm <- sum(sapply(G40_bm, function(vec) {
        in_50 <- any(abs(vec-50) <= eta); 
        in_50 }))/max_iter
      
      # --------------------------------------------------- #
      # FDR:
      FDR_G20_both <- mean(sapply(G20_both, function(vec) {
        if (length(vec) < 1) { 0 } 
        else if (all( abs(vec-50) > eta)) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_G40_both <- mean(sapply(G40_both, function(vec) {
        if (length(vec) < 1) { 0 } 
        else if (all( abs(vec-50) > eta)) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_PELT <- mean(sapply(PELT, function(vec) {
        if (length(vec) < 1) { 0 } 
        else if (all( abs(vec-50) > eta)) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_stepR <- mean(sapply(stepR, function(vec) {
        if (length(vec) < 1) { 0 } 
        else if (all( abs(vec-50) > eta)) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_wbs <- mean(sapply(wbs, function(vec) {
        if (length(vec) < 1) { 0 } 
        else if (all( abs(vec-50) > eta)) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_G20_bm <- mean(sapply(G20_bm, function(vec) {
        if (length(vec) < 1) { 0 } 
        else if (all( abs(vec-50) > eta)) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_G40_bm <- mean(sapply(G40_bm, function(vec) {
        if (length(vec) < 1) { 0 } 
        else if (all( abs(vec-50) > eta)) { 1 } 
        else { 0 } }), na.rm = TRUE)
      
      
    }else if (i == 2){
      # --------------------------------------------------- #
      # Power:
      power_G20_both <- sum(sapply(G20_both, function(vec) {
        in_40 <- any(abs(vec-40) <= eta); 
        in_60 <- any(abs(vec-60) <= eta);
        in_40 & in_60 }))/max_iter
      power_G40_both <- sum(sapply(G40_both, function(vec) {
        in_40 <- any(abs(vec-40) <= eta); 
        in_60 <- any(abs(vec-60) <= eta);
        in_40 & in_60 }))/max_iter
      power_PELT <- sum(sapply(PELT, function(vec) {
        in_40 <- any(abs(vec-40) <= eta); 
        in_60 <- any(abs(vec-60) <= eta);
        in_40 & in_60 }))/max_iter
      power_stepR <- sum(sapply(stepR, function(vec) {
        in_40 <- any(abs(vec-40) <= eta); 
        in_60 <- any(abs(vec-60) <= eta);
        in_40 & in_60 }))/max_iter
      power_wbs <- sum(sapply(wbs, function(vec) {
        in_40 <- any(abs(vec-40) <= eta); 
        in_60 <- any(abs(vec-60) <= eta);
        in_40 & in_60 }))/max_iter
      power_G20_bm <- sum(sapply(G20_bm, function(vec) {
        in_40 <- any(abs(vec-40) <= eta); 
        in_60 <- any(abs(vec-60) <= eta);
        in_40 & in_60 }))/max_iter
      power_G40_bm <- sum(sapply(G40_bm, function(vec) {
        in_40 <- any(abs(vec-40) <= eta); 
        in_60 <- any(abs(vec-60) <= eta);
        in_40 & in_60 }))/max_iter
      
      # --------------------------------------------------- #
      # FDR:
      FDR_G20_both <- mean(sapply(G20_both, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-40) > eta))&(all( abs(vec-60) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_G40_both <- mean(sapply(G40_both, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-40) > eta))&(all( abs(vec-60) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_PELT <- mean(sapply(PELT, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-40) > eta))&(all( abs(vec-60) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_stepR <- mean(sapply(stepR, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-40) > eta))&(all( abs(vec-60) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_wbs <- mean(sapply(wbs, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-40) > eta))&(all( abs(vec-60) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_G20_bm <- mean(sapply(G20_bm, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-40) > eta))&(all( abs(vec-60) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_G40_bm <- mean(sapply(G40_bm, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-40) > eta))&(all( abs(vec-60) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      
      
    }else if (i == 3){
      # --------------------------------------------------- #
      # Power:
      power_G20_both <- sum(sapply(G20_both, function(vec) {
        in_25 <- any(abs(vec-25) <= eta); 
        in_50 <- any(abs(vec-50) <= eta);
        in_75 <- any(abs(vec-75) <= eta);
        in_25 & in_50 & in_75 }))/max_iter
      power_G40_both <- sum(sapply(G40_both, function(vec) {
        in_25 <- any(abs(vec-25) <= eta); 
        in_50 <- any(abs(vec-50) <= eta);
        in_75 <- any(abs(vec-75) <= eta);
        in_25 & in_50 & in_75 }))/max_iter
      power_PELT <- sum(sapply(PELT, function(vec) {
        in_25 <- any(abs(vec-25) <= eta); 
        in_50 <- any(abs(vec-50) <= eta);
        in_75 <- any(abs(vec-75) <= eta);
        in_25 & in_50 & in_75 }))/max_iter
      power_stepR <- sum(sapply(stepR, function(vec) {
        in_25 <- any(abs(vec-25) <= eta); 
        in_50 <- any(abs(vec-50) <= eta);
        in_75 <- any(abs(vec-75) <= eta);
        in_25 & in_50 & in_75 }))/max_iter
      power_wbs <- sum(sapply(wbs, function(vec) {
        in_25 <- any(abs(vec-25) <= eta); 
        in_50 <- any(abs(vec-50) <= eta);
        in_75 <- any(abs(vec-75) <= eta);
        in_25 & in_50 & in_75 }))/max_iter
      power_G20_bm <- sum(sapply(G20_bm, function(vec) {
        in_25 <- any(abs(vec-25) <= eta); 
        in_50 <- any(abs(vec-50) <= eta);
        in_75 <- any(abs(vec-75) <= eta);
        in_25 & in_50 & in_75 }))/max_iter
      power_G40_bm <- sum(sapply(G40_bm, function(vec) {
        in_25 <- any(abs(vec-25) <= eta); 
        in_50 <- any(abs(vec-50) <= eta);
        in_75 <- any(abs(vec-75) <= eta);
        in_25 & in_50 & in_75 }))/max_iter
      
      # --------------------------------------------------- #
      # FDR:
      FDR_G20_both <- mean(sapply(G20_both, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-25) > eta))&(all( abs(vec-50) > eta))&(all( abs(vec-75) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_G40_both <- mean(sapply(G40_both, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-25) > eta))&(all( abs(vec-50) > eta))&(all( abs(vec-75) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_PELT <- mean(sapply(PELT, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-25) > eta))&(all( abs(vec-50) > eta))&(all( abs(vec-75) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_stepR <- mean(sapply(stepR, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-25) > eta))&(all( abs(vec-50) > eta))&(all( abs(vec-75) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_wbs <- mean(sapply(wbs, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-25) > eta))&(all( abs(vec-50) > eta))&(all( abs(vec-75) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_G20_bm <- mean(sapply(G20_bm, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-25) > eta))&(all( abs(vec-50) > eta))&(all( abs(vec-75) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_G40_bm <- mean(sapply(G40_bm, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-25) > eta))&(all( abs(vec-50) > eta))&(all( abs(vec-75) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      
    }
    
    # display:
    cat("num of jumps",i,"Case",j,"produced.\n")
    cat("Power:\n")
    cat("Joint-MOSUM (G=20),Joint-MOSUM (G=40),PELT,stepR,wbs,MOSUM (G=20),MOSUM (G=40) \n")
    print(round(c(power_G20_both,power_G40_both,
                  power_PELT,power_stepR,power_wbs,
                  power_G20_bm,power_G40_bm),3))
    cat("FDR:\n")
    cat("Joint-MOSUM (G=20),Joint-MOSUM (G=40),PELT,stepR,wbs,MOSUM (G=20),MOSUM (G=40) \n")
    print(round(c(FDR_G20_both,FDR_G40_both,
                  FDR_PELT,FDR_stepR,FDR_wbs,
                  FDR_G20_bm,FDR_G40_bm),3))
    cat("#---------------------------------------------- #\n")
  }
}



##################################################################
# Numerical Experiment 1 (bivariate):
##################################################################
get_sim <- function(i,j,iter_max){
  
  n <- 100
  if (i == 1){
    num_jump1 <- 1
    start1 <- c(1,61)
    end1 <- c(60,n)
    
    num_jump2 <- 1
    start2 <- c(1,41)
    end2 <- c(40,n)
  }else if (i == 2){
    num_jump1 <- 1
    start1 <- c(1,51)
    end1 <- c(50,n)
    
    num_jump2 <- 2
    start2 <- c(1,41,61)
    end2 <- c(40,60,n)
  }else if (i == 3){
    num_jump1 <- 2
    start1 <- c(1,41,61)
    end1 <- c(40,60,n)
    
    num_jump2 <- 2
    start2 <- c(1,31,71)
    end2 <- c(30,70,n)
  }
  
  if (j == 1){
    level_range <- 2
    sigma_min <- 0.1
    sigma_max <- 0.4
  }else if (j == 2){
    level_range <- 2
    sigma_min <- 0.1
    sigma_max <- 0.8
  }else if (j == 3){
    level_range <- 2
    sigma_min <- 0.4
    sigma_max <- 0.8
  }else if (j == 4){
    level_range <- 1
    sigma_min <- 0.1
    sigma_max <- 0.4
  }else if (j == 5){
    level_range <- 1
    sigma_min <- 0.1
    sigma_max <- 0.8
  }else if (j == 6){
    level_range <- 1
    sigma_min <- 0.4
    sigma_max <- 0.8
  }
  
  library(doParallel)
  require(dqrng)
  require(doRNG)
  
  # cl <- makeCluster(50)
  cl <- makeCluster(7)
  registerDoParallel(cl)
  
  set.seed(123)
  list_iter <- foreach(iter = 1:iter_max,
                       .errorhandling = 'remove',
                       .packages = c("mvtnorm","stringr",
                                     "ecp","hdbinseg","changepoints")) %dopar% {
                                       
                                       source("Library.R")
                                       
                                       # ------------------------------------------- #
                                       n <- 100
                                       series1 <- series_noisy1 <- rep(NA,n)
                                       levels_vec1 <- runif((num_jump1+1), min=-level_range, max=level_range)
                                       sigma_vec1 <- runif((num_jump1+1), min=sigma_min, max=sigma_max)
                                       for (k in 1:(num_jump1+1)) {
                                         series1[start1[k]:end1[k]] <- levels_vec1[k]
                                         series_noisy1[start1[k]:end1[k]] <- series1[start1[k]:end1[k]] + rnorm(end1[k]-start1[k]+1, mean=0, sd=sigma_vec1[k])
                                       }
                                       
                                       series2 <- series_noisy2 <- rep(NA,n)
                                       levels_vec2 <- runif((num_jump2+1), min=-level_range, max=level_range)
                                       sigma_vec2 <- runif((num_jump2+1), min=sigma_min, max=sigma_max)
                                       for (k in 1:(num_jump2+1)) {
                                         series2[start2[k]:end2[k]] <- levels_vec2[k]
                                         series_noisy2[start2[k]:end2[k]] <- series2[start2[k]:end2[k]] + rnorm(end2[k]-start2[k]+1, mean=0, sd=sigma_vec2[k])
                                       }
                                       
                                       # --------------------------------------------- #
                                       # Bi-MOSUM:
                                       Bi_MOSUM_G20 <- bivariate_mosum(series_noisy1,series_noisy2,
                                                                       n=n,Gl=20,Gr=20,alpha=0.05,eta=0.2)
                                       Bi_MOSUM_G20 <-sort(unique(c(Bi_MOSUM_G20$YX$cpt_eta,
                                                                    Bi_MOSUM_G20$YX2$cpt_eta,
                                                                    Bi_MOSUM_G20$Y2X$cpt_eta,
                                                                    Bi_MOSUM_G20$Y2X2$cpt_eta)))
                                       Bi_MOSUM_G40 <- bivariate_mosum(series_noisy1,series_noisy2,
                                                                       n=n,Gl=40,Gr=40,alpha=0.05,eta=0.2)
                                       Bi_MOSUM_G40 <-sort(unique(c(Bi_MOSUM_G40$YX$cpt_eta,
                                                                    Bi_MOSUM_G40$YX2$cpt_eta,
                                                                    Bi_MOSUM_G40$Y2X$cpt_eta,
                                                                    Bi_MOSUM_G40$Y2X2$cpt_eta)))
                                       
                                       # --------------------------------------------- #
                                       # ecp: 
                                       ECP_alp1 <- setdiff(e.divisive(cbind(series_noisy1,series_noisy2),sig.lvl=0.05,
                                                                      R=100,alpha=1)$estimates,c(1,(n+1)))
                                       ECP_alp2 <- setdiff(e.divisive(cbind(series_noisy1,series_noisy2),sig.lvl=0.05,
                                                                      R=100,alpha=2)$estimates,c(1,(n+1)))
                                       ECP <- sort(unique(c(ECP_alp1,ECP_alp2)))
                                       
                                       # --------------------------------------------- #
                                       # SBS:
                                       SBS_alp1 <- sbs.alg(rbind(series_noisy1,series_noisy2),cp.type=1)$ecp
                                       SBS_alp2 <- sbs.alg(rbind(series_noisy1,series_noisy2),cp.type=2)$ecp
                                       SBS <- sort(unique(c(SBS_alp1,SBS_alp2)))
                                       
                                       # --------------------------------------------- #
                                       # DCBS:
                                       DCBS_alp1 <- dcbs.alg(rbind(series_noisy1,series_noisy2),cp.type=1)$ecp
                                       DCBS_alp2 <- dcbs.alg(rbind(series_noisy1,series_noisy2),cp.type=2)$ecp
                                       DCBS <- sort(unique(c(DCBS_alp1,DCBS_alp2)))
                                       
                                       # --------------------------------------------- #
                                       # MNP:
                                       intervals <- WBS.intervals(M=10,lower=1,upper=n) 
                                       temp <- WBS.multi.nonpar(Y=rbind(series_noisy1,series_noisy2), 
                                                                W=rbind(series_noisy1,series_noisy2), 1, n, 
                                                                intervals$Alpha, intervals$Beta, 
                                                                h=5*(30*log(n)/n)^{1/2}, delta = 10)
                                       MNP <- sort(unique(thresholdBS(temp, median(temp$Dval))$cpt_hat$location))
                                       
                                       # --------------------------------------------- #
                                       output <- list()
                                       output$Bi_MOSUM_G20 <- Bi_MOSUM_G20
                                       output$Bi_MOSUM_G40 <- Bi_MOSUM_G40
                                       output$SBS <- SBS
                                       output$DCBS <- DCBS
                                       output$ECP <- ECP
                                       output$MNP <- MNP
                                       output
                                     }
  stopCluster(cl)
  print(length(list_iter))
  if (length(list_iter) == 0) {
    warning("No successful iterations were returned.")
  }
  save(list_iter,file=paste0("bivariate_jump",i,"_snr",j,".RData"))
}


# --------------------------------------------------- #
# Run iterations
for (i in 1:3){
  for (j in 1:6){
    
    cat("The current iteration is i:",i,"j:",j,"\n")
    get_sim(i,j,iter_max=500)
    cat("Done.\n") 
  }
}


# --------------------------------------------------- #
# Results:
max_iter <- 500
for (i in 1:3){
  for (j in 1:6){
    
    load(paste0("bivariate_jump",i,"_snr",j,".RData"))
    max_iter <- length(list_iter)
    
    # --------------------------------------------------- #
    G20_both <- mapply(iter=1:max_iter,function(iter)list_iter[[iter]]$Bi_MOSUM_G20)
    G40_both <- mapply(iter=1:max_iter,function(iter)list_iter[[iter]]$Bi_MOSUM_G40)
    SBS <- mapply(iter=1:max_iter,function(iter)list_iter[[iter]]$SBS)
    DCBS <- mapply(iter=1:max_iter,function(iter)list_iter[[iter]]$DCBS)
    ECP <- mapply(iter=1:max_iter,function(iter)list_iter[[iter]]$ECP)
    MNP <- mapply(iter=1:max_iter,function(iter)list_iter[[iter]]$MNP)
    
    if (i == 1){
      # --------------------------------------------------- #
      # Power:
      eta <- 5
      power_G20_both <- sum(sapply(G20_both, function(vec) {
        in_40 <- any(abs(vec-40) <= eta); 
        in_60 <- any(abs(vec-60) <= eta);
        in_40 & in_60 }))/max_iter
      power_G40_both <- sum(sapply(G40_both, function(vec) {
        in_40 <- any(abs(vec-40) <= eta); 
        in_60 <- any(abs(vec-60) <= eta);
        in_40 & in_60 }))/max_iter
      power_SBS <- sum(sapply(SBS, function(vec) {
        in_40 <- any(abs(vec-40) <= eta); 
        in_60 <- any(abs(vec-60) <= eta);
        in_40 & in_60 }))/max_iter
      power_DCBS <- sum(sapply(DCBS, function(vec) {
        in_40 <- any(abs(vec-40) <= eta); 
        in_60 <- any(abs(vec-60) <= eta);
        in_40 & in_60 }))/max_iter
      power_ECP <- sum(sapply(ECP, function(vec) {
        in_40 <- any(abs(vec-40) <= eta); 
        in_60 <- any(abs(vec-60) <= eta);
        in_40 & in_60 }))/max_iter
      power_MNP <- sum(sapply(MNP, function(vec) {
        in_40 <- any(abs(vec-40) <= eta); 
        in_60 <- any(abs(vec-60) <= eta);
        in_40 & in_60 }))/max_iter
      
      # --------------------------------------------------- #
      # FDR:
      eta <- 0
      FDR_G20_both <- mean(sapply(G20_both, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-40) > eta))&(all( abs(vec-60) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_G40_both <- mean(sapply(G40_both, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-40) > eta))&(all( abs(vec-60) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_SBS <- mean(sapply(SBS, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-40) > eta))&(all( abs(vec-60) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_DCBS <- mean(sapply(DCBS, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-40) > eta))&(all( abs(vec-60) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_ECP <- mean(sapply(ECP, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-40) > eta))&(all( abs(vec-60) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_MNP <- mean(sapply(MNP, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-40) > eta))&(all( abs(vec-60) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      
      
    }else if (i == 2){
      # --------------------------------------------------- #
      # Power:
      eta <- 5
      power_G20_both <- sum(sapply(G20_both, function(vec) {
        in_40 <- any(abs(vec-40) <= eta); 
        in_50 <- any(abs(vec-50) <= eta);
        in_60 <- any(abs(vec-60) <= eta);
        in_40 & in_50 & in_60 }))/max_iter
      power_G40_both <- sum(sapply(G40_both, function(vec) {
        in_40 <- any(abs(vec-40) <= eta); 
        in_50 <- any(abs(vec-50) <= eta);
        in_60 <- any(abs(vec-60) <= eta);
        in_40 & in_50 & in_60 }))/max_iter
      power_SBS <- sum(sapply(SBS, function(vec) {
        in_40 <- any(abs(vec-40) <= eta); 
        in_50 <- any(abs(vec-50) <= eta);
        in_60 <- any(abs(vec-60) <= eta);
        in_40 & in_50 & in_60 }))/max_iter
      power_DCBS <- sum(sapply(DCBS, function(vec) {
        in_40 <- any(abs(vec-40) <= eta); 
        in_50 <- any(abs(vec-50) <= eta);
        in_60 <- any(abs(vec-60) <= eta);
        in_40 & in_50 & in_60 }))/max_iter
      power_ECP <- sum(sapply(ECP, function(vec) {
        in_40 <- any(abs(vec-40) <= eta); 
        in_50 <- any(abs(vec-50) <= eta);
        in_60 <- any(abs(vec-60) <= eta);
        in_40 & in_50 & in_60 }))/max_iter
      power_MNP <- sum(sapply(MNP, function(vec) {
        in_40 <- any(abs(vec-40) <= eta); 
        in_50 <- any(abs(vec-50) <= eta);
        in_60 <- any(abs(vec-60) <= eta);
        in_40 & in_50 & in_60 }))/max_iter
      
      # --------------------------------------------------- #
      # FDR:
      eta <- 0
      FDR_G20_both <- mean(sapply(G20_both, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-40) > eta))&(all( abs(vec-50) > eta))&(all( abs(vec-60) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_G40_both <- mean(sapply(G40_both, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-40) > eta))&(all( abs(vec-50) > eta))&(all( abs(vec-60) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_SBS <- mean(sapply(SBS, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-40) > eta))&(all( abs(vec-50) > eta))&(all( abs(vec-60) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_DCBS <- mean(sapply(DCBS, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-40) > eta))&(all( abs(vec-50) > eta))&(all( abs(vec-60) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_ECP <- mean(sapply(ECP, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-40) > eta))&(all( abs(vec-50) > eta))&(all( abs(vec-60) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_MNP <- mean(sapply(MNP, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-40) > eta))&(all( abs(vec-50) > eta))&(all( abs(vec-60) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      
    }else if (i == 3){
      # --------------------------------------------------- #
      # Power:
      eta <- 5
      power_G20_both <- sum(sapply(G20_both, function(vec) {
        in_30 <- any(abs(vec-30) <= eta);
        in_40 <- any(abs(vec-40) <= eta); 
        in_60 <- any(abs(vec-60) <= eta);
        in_70 <- any(abs(vec-70) <= eta);
        in_30 & in_40 & in_60 & in_70 }))/max_iter
      power_G40_both <- sum(sapply(G40_both, function(vec) {
        in_30 <- any(abs(vec-30) <= eta);
        in_40 <- any(abs(vec-40) <= eta); 
        in_60 <- any(abs(vec-60) <= eta);
        in_70 <- any(abs(vec-70) <= eta);
        in_30 & in_40 & in_60 & in_70 }))/max_iter
      power_SBS <- sum(sapply(SBS, function(vec) {
        in_30 <- any(abs(vec-30) <= eta);
        in_40 <- any(abs(vec-40) <= eta); 
        in_60 <- any(abs(vec-60) <= eta);
        in_70 <- any(abs(vec-70) <= eta);
        in_30 & in_40 & in_60 & in_70 }))/max_iter
      power_DCBS <- sum(sapply(DCBS, function(vec) {
        in_30 <- any(abs(vec-30) <= eta);
        in_40 <- any(abs(vec-40) <= eta); 
        in_60 <- any(abs(vec-60) <= eta);
        in_70 <- any(abs(vec-70) <= eta);
        in_30 & in_40 & in_60 & in_70 }))/max_iter
      power_ECP <- sum(sapply(ECP, function(vec) {
        in_30 <- any(abs(vec-30) <= eta);
        in_40 <- any(abs(vec-40) <= eta); 
        in_60 <- any(abs(vec-60) <= eta);
        in_70 <- any(abs(vec-70) <= eta);
        in_30 & in_40 & in_60 & in_70 }))/max_iter
      power_MNP <- sum(sapply(MNP, function(vec) {
        in_30 <- any(abs(vec-30) <= eta);
        in_40 <- any(abs(vec-40) <= eta); 
        in_60 <- any(abs(vec-60) <= eta);
        in_70 <- any(abs(vec-70) <= eta);
        in_30 & in_40 & in_60 & in_70 }))/max_iter
      
      # --------------------------------------------------- #
      # FDR:
      eta <- 0
      FDR_G20_both <- mean(sapply(G20_both, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-30) > eta))&(all( abs(vec-40) > eta))&(all( abs(vec-60) > eta))&(all( abs(vec-70) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_G40_both <- mean(sapply(G40_both, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-30) > eta))&(all( abs(vec-40) > eta))&(all( abs(vec-60) > eta))&(all( abs(vec-70) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_SBS <- mean(sapply(SBS, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-30) > eta))&(all( abs(vec-40) > eta))&(all( abs(vec-60) > eta))&(all( abs(vec-70) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_DCBS <- mean(sapply(DCBS, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-30) > eta))&(all( abs(vec-40) > eta))&(all( abs(vec-60) > eta))&(all( abs(vec-70) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_ECP <- mean(sapply(ECP, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-30) > eta))&(all( abs(vec-40) > eta))&(all( abs(vec-60) > eta))&(all( abs(vec-70) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      FDR_MNP <- mean(sapply(MNP, function(vec) {
        if (length(vec) < 1) { NA } 
        else if ((all( abs(vec-30) > eta))&(all( abs(vec-40) > eta))&(all( abs(vec-60) > eta))&(all( abs(vec-70) > eta))) { 1 } 
        else { 0 } }), na.rm = TRUE)
      
    }
    
    # display:
    cat("num of jumps",(i+1),"Case",j,"produced.\n")
    cat("Power:\n")
    cat("Bi-MOSUM (G=20),Bi-MOSUM (G=40),SBS,DCBS,ECP,MNP \n")
    print(round(c(power_G20_both,power_G40_both,
                  power_SBS,power_DCBS,power_ECP,power_MNP),3))
    cat("FDR:\n")
    cat("Bi-MOSUM (G=20),Bi-MOSUM (G=40),SBS,DCBS,ECP,MNP \n")
    print(round(c(FDR_G20_both,FDR_G40_both,
                  FDR_SBS,FDR_DCBS,FDR_ECP,FDR_MNP),3))
    cat("#---------------------------------------------- #\n")
  }
}


##################################################################
# Numerical Experiment 2 (hotspots):
##################################################################

get_sim <- function(i,j,delta,iter_max){
  
  n <- 100
  alpha <- 0.05
  if (i == 1){
    # X:
    num_jump1 <- 1
    start1 <- c(1,51)
    end1 <- c(50,n)
    
    # Y:
    num_jump2 <- 1
    start2 <- c(1,(51+delta))
    end2 <- c((50+delta),n)
  }else if (i == 2){
    # X:
    num_jump1 <- 2
    start1 <- c(1,41,61)
    end1 <- c(40,60,n)
    
    # Y:
    num_jump2 <- 2
    start2 <- c(1,(41+delta),(61+delta))
    end2 <- c((40+delta),(60+delta),n)
  }
  
  if (j == 1){
    level_range <- 4
    sigma_min <- 0.1
    sigma_max <- 0.4
  }else if (j == 2){
    level_range <- 4
    sigma_min <- 0.1
    sigma_max <- 0.8
  }else if (j == 3){
    level_range <- 4
    sigma_min <- 0.4
    sigma_max <- 0.8
  }else if (j == 4){
    level_range <- 2
    sigma_min <- 0.1
    sigma_max <- 0.4
  }else if (j == 5){
    level_range <- 2
    sigma_min <- 0.1
    sigma_max <- 0.8
  }else if (j == 6){
    level_range <- 2
    sigma_min <- 0.4
    sigma_max <- 0.8
  }
  
  
  library(doParallel)
  require(dqrng)
  require(doRNG)
  
  cl <- makeCluster(50)
  registerDoParallel(cl)
  
  set.seed(123)
  list_iter <- foreach(iter = 1:iter_max,
                       .errorhandling = 'remove',
                       .packages = c("mvtnorm","stringr")) %dopar% {
                         
                         source("Library.R")
                         # ------------------------------------------- #
                         n <- 100
                         series1 <- series_noisy1 <- rep(NA,n)
                         levels_vec1 <- runif((num_jump1+1), min=-level_range, max=level_range)
                         sigma_vec1 <- runif((num_jump1+1), min=sigma_min, max=sigma_max)
                         for (k in 1:(num_jump1+1)) {
                           series1[start1[k]:end1[k]] <- levels_vec1[k]
                           series_noisy1[start1[k]:end1[k]] <- series1[start1[k]:end1[k]] + rnorm(end1[k]-start1[k]+1, mean=0, sd=sigma_vec1[k])
                         }
                         
                         series2 <- series_noisy2 <- rep(NA,n)
                         levels_vec2 <- runif((num_jump2+1), min=-level_range, max=level_range)
                         sigma_vec2 <- runif((num_jump2+1), min=sigma_min, max=sigma_max)
                         for (k in 1:(num_jump2+1)) {
                           series2[start2[k]:end2[k]] <- levels_vec2[k]
                           series_noisy2[start2[k]:end2[k]] <- series2[start2[k]:end2[k]] + rnorm(end2[k]-start2[k]+1, mean=0, sd=sigma_vec2[k])
                         }
                         
                         
                         # --------------------------------------------- #
                         # Bi-MOSUM:
                         Bi_MOSUM_G20 <- bivariate_mosum(series_noisy1,series_noisy2,
                                                         n=n,Gl=20,Gr=20,alpha=0.05,eta=0.2)
                         Bi_MOSUM_G40 <- bivariate_mosum(series_noisy1,series_noisy2,
                                                         n=n,Gl=40,Gr=40,alpha=0.05,eta=0.2)
                         
                         # --------------------------------------------- #
                         # Threshold:
                         hotspot_G20_thres <- get_intervals(sort(intersect(unique(c(Bi_MOSUM_G20$YX$idx_cand,
                                                                                    Bi_MOSUM_G20$YX2$idx_cand,
                                                                                    Bi_MOSUM_G20$Y2X$idx_cand,
                                                                                    Bi_MOSUM_G20$Y2X2$idx_cand)),
                                                                           Bi_MOSUM_G20$Y_both$idx_cand)))
                         hotspot_G40_thres <- get_intervals(sort(intersect(unique(c(Bi_MOSUM_G40$YX$idx_cand,
                                                                                    Bi_MOSUM_G40$YX2$idx_cand,
                                                                                    Bi_MOSUM_G40$Y2X$idx_cand,
                                                                                    Bi_MOSUM_G40$Y2X2$idx_cand)),
                                                                           Bi_MOSUM_G40$Y_both$idx_cand)))
                         
                         # --------------------------------------------- #
                         # CI:
                         # On Y:
                         if (!is.null(Bi_MOSUM_G20$Y_both$cpt_eta)){
                           CI_Y_G20 <- joint_mosum_ci(series_noisy2,refinement=Bi_MOSUM_G20$Y_both,
                                                      method="eta",n,Gl=20,Gr=20,N=1000,
                                                      alpha=0.05,epsilon=0.2,eta=0.2)
                         }else {
                           CI_Y_G20 <- NULL
                         }
                         if (!is.null(Bi_MOSUM_G40$Y_both$cpt_eta)){
                           CI_Y_G40 <- joint_mosum_ci(series_noisy2,refinement=Bi_MOSUM_G40$Y_both,
                                                      method="eta",n,Gl=40,Gr=40,N=1000,
                                                      alpha=0.05,epsilon=0.2,eta=0.2)
                         }else {
                           CI_Y_G40 <- NULL
                         }
                         
                         # On YX:
                         if (!is.null(Bi_MOSUM_G20$YX$cpt_eta)){
                           CI_YX_G20 <- bivariate_mosum_ci(series_noisy2,series_noisy1,
                                                           bi_refinement=Bi_MOSUM_G20$YX,
                                                           moment="mean_mean",
                                                           method="eta",n,Gl=20,Gr=20,N=1000,
                                                           alpha=0.05,epsilon=0.2,eta=0.2)
                         }else{
                           CI_YX_G20 <- NULL
                         }
                         if (!is.null(Bi_MOSUM_G40$YX$cpt_eta)){
                           CI_YX_G40 <- bivariate_mosum_ci(series_noisy2,series_noisy1,
                                                           bi_refinement=Bi_MOSUM_G40$YX,
                                                           moment="mean_mean",
                                                           method="eta",n,Gl=40,Gr=40,N=1000,
                                                           alpha=0.05,epsilon=0.2,eta=0.2)
                         }else{
                           CI_YX_G40 <- NULL
                         }
                         
                         # On YX2:
                         if (!is.null(Bi_MOSUM_G20$YX2$cpt_eta)){
                           CI_YX2_G20 <- bivariate_mosum_ci(series_noisy2,series_noisy1,
                                                            bi_refinement=Bi_MOSUM_G20$YX2,
                                                            moment="mean_var",
                                                            method="eta",n,Gl=20,Gr=20,N=1000,
                                                            alpha=0.05,epsilon=0.2,eta=0.2)
                         }else{
                           CI_YX2_G20 <- NULL
                         }
                         if (!is.null(Bi_MOSUM_G40$YX2$cpt_eta)){
                           CI_YX2_G40 <- bivariate_mosum_ci(series_noisy2,series_noisy1,
                                                            bi_refinement=Bi_MOSUM_G40$YX2,
                                                            moment="mean_var",
                                                            method="eta",n,Gl=40,Gr=40,N=1000,
                                                            alpha=0.05,epsilon=0.2,eta=0.2)
                         }else{
                           CI_YX2_G40 <- NULL
                         }
                         
                         # On Y2X:
                         if (!is.null(Bi_MOSUM_G20$Y2X$cpt_eta)){
                           CI_Y2X_G20 <- bivariate_mosum_ci(series_noisy2,series_noisy1,
                                                            bi_refinement=Bi_MOSUM_G20$Y2X,
                                                            moment="var_mean",
                                                            method="eta",n,Gl=20,Gr=20,N=1000,
                                                            alpha=0.05,epsilon=0.2,eta=0.2)
                         }else{
                           CI_Y2X_G20 <- NULL
                         }
                         if (!is.null(Bi_MOSUM_G40$Y2X$cpt_eta)){
                           CI_Y2X_G40 <- bivariate_mosum_ci(series_noisy2,series_noisy1,
                                                            bi_refinement=Bi_MOSUM_G40$Y2X,
                                                            moment="var_mean",
                                                            method="eta",n,Gl=20,Gr=20,N=1000,
                                                            alpha=0.05,epsilon=0.2,eta=0.2)
                         }else{
                           CI_Y2X_G40 <- NULL
                         }
                         
                         # On Y2X2:
                         if (!is.null(Bi_MOSUM_G20$Y2X2$cpt_eta)){
                           CI_Y2X2_G20 <- bivariate_mosum_ci(series_noisy2,series_noisy1,
                                                             bi_refinement=Bi_MOSUM_G20$Y2X2,
                                                             moment="var_var",
                                                             method="eta",n,Gl=20,Gr=20,N=1000,
                                                             alpha=0.05,epsilon=0.2,eta=0.2)
                         }else{
                           CI_Y2X2_G20 <- NULL
                         }
                         if (!is.null(Bi_MOSUM_G40$Y2X2$cpt_eta)){
                           CI_Y2X2_G40 <- bivariate_mosum_ci(series_noisy2,series_noisy1,
                                                             bi_refinement=Bi_MOSUM_G40$Y2X2,
                                                             moment="var_var",
                                                             method="eta",n,Gl=40,Gr=40,N=1000,
                                                             alpha=0.05,epsilon=0.2,eta=0.2)
                         }else{
                           CI_Y2X2_G40 <- NULL
                         }
                         
                         
                         # --------------------------------------------- #
                         # G=20:
                         if (is.null(CI_Y_G20$ptwise)){
                           hotspot_G20_ptwise <- NULL
                         }else{
                           list_tmp <- CI_Y_G20$ptwise
                           flag_G20_ptwise <- FALSE
                           if (!is.null(CI_YX_G20$ptwise)){
                             list_tmp <- get_union_of_intersects(CI_YX_G20$ptwise,list_tmp,type="ptwise")
                             flag_G20_ptwise <- TRUE
                           }
                           if (!is.null(CI_YX2_G20$ptwise)){
                             list_tmp <- get_union_of_intersects(CI_YX2_G20$ptwise,list_tmp,type="ptwise")
                             flag_G20_ptwise <- TRUE
                           }
                           if (!is.null(CI_Y2X_G20$ptwise)){
                             list_tmp <- get_union_of_intersects(CI_Y2X_G20$ptwise,list_tmp,type="ptwise")
                             flag_G20_ptwise <- TRUE
                           }
                           if (!is.null(CI_Y2X2_G20$ptwise)){
                             list_tmp <- get_union_of_intersects(CI_Y2X2_G20$ptwise,list_tmp,type="ptwise")
                             flag_G20_ptwise <- TRUE
                           }
                           
                           if (flag_G20_ptwise == TRUE){
                             hotspot_G20_ptwise <- list_tmp
                           }else{
                             hotspot_G20_ptwise <- NULL
                           }
                         }
                         
                         # --------------------------------------------- #
                         # G=40:
                         if (is.null(CI_Y_G40$ptwise)){
                           hotspot_G40_ptwise <- NULL
                         }else{
                           list_tmp <- CI_Y_G40$ptwise
                           flag_G40_ptwise <- FALSE
                           if (!is.null(CI_YX_G40$ptwise)){
                             list_tmp <- get_union_of_intersects(CI_YX_G40$ptwise,list_tmp,type="ptwise")
                             flag_G40_ptwise <- TRUE
                           }
                           if (!is.null(CI_YX2_G40$ptwise)){
                             list_tmp <- get_union_of_intersects(CI_YX2_G40$ptwise,list_tmp,type="ptwise")
                             flag_G40_ptwise <- TRUE
                           }
                           if (!is.null(CI_Y2X_G40$ptwise)){
                             list_tmp <- get_union_of_intersects(CI_Y2X_G40$ptwise,list_tmp,type="ptwise")
                             flag_G40_ptwise <- TRUE
                           }
                           if (!is.null(CI_Y2X2_G40$ptwise)){
                             list_tmp <- get_union_of_intersects(CI_Y2X2_G40$ptwise,list_tmp,type="ptwise")
                             flag_G40_ptwise <- TRUE
                           }
                           
                           if (flag_G40_ptwise == TRUE){
                             hotspot_G40_ptwise <- list_tmp
                           }else{
                             hotspot_G40_ptwise <- NULL
                           }
                         }
                         
                         
                         # --------------------------------------------- #
                         output <- list()
                         output$Bi_MOSUM_G20 <- Bi_MOSUM_G20
                         output$Bi_MOSUM_G40 <- Bi_MOSUM_G40
                         
                         output$hotspot_G20_thres <- hotspot_G20_thres
                         output$hotspot_G40_thres <- hotspot_G40_thres
                         
                         output$CI_Y_G20 <- CI_Y_G20
                         output$CI_YX_G20 <- CI_YX_G20
                         output$CI_Y_G40 <- CI_Y_G40
                         output$CI_YX_G40 <- CI_YX_G40
                         
                         output$hotspot_G20_ptwise  <- hotspot_G20_ptwise
                         output$hotspot_G40_ptwise  <- hotspot_G40_ptwise
                         
                         output
                       }
  stopCluster(cl)
  print(length(list_iter))
  if (length(list_iter) == 0) {
    warning("No successful iterations were returned.")
  }
  save(list_iter,file=paste0("hotspots_jump",i,"_snr",j,"_delta",delta,".RData"))
}


# --------------------------------------------------- #
# Run iterations
delta_list <- c(0,5)
for (i in 1:2){
  for (j in 1:6){
    for (k in 1:2){
      
      delta <- delta_list[k]
      cat("The current iteration is i:",i,"j:",j,"delta",delta,"\n")
      get_sim(i,j,delta,iter_max=500)
      cat("Done.\n") 
    }
  }
}


# --------------------------------------------------- #
# Run iterations
eta <- 0
delta_list <- c(0,5)
for (i in 1:2){
  for (k in c(1,2)){
    delta <- delta_list[k]
    
    for (j in 1:6){
      
      load(paste0("hotspots_jump",i,"_snr",j,"_delta",delta,".RData"))
      max_iter <- length(list_iter)
      
      n <- 100
      if (i == 1){
        true_Y <- (50+delta)
      }else if (i == 2){
        true_Y <- c((40+delta),(60+delta))
      }
      
      
      # ------------------------------------ #
      # Thrs (G=20)
      ratio_G20_thrs <- rep(NA,max_iter)
      length_G20_thrs <- rep(0,max_iter)
      for (iter in 1:max_iter){
        tmp <- list_iter[[iter]]$hotspot_G20_thres
        if (is.null(tmp) || (is.list(tmp) && length(tmp) == 0)){
          ratio_G20_thrs[iter] <- 0
          length_G20_thrs[iter] <- NA
        }else if (is.null(tmp[[1]]) || length(tmp[[1]]) == 0 || (length(tmp[[1]]) == 1 && is.na(tmp[[1]]))){
          ratio_G20_thrs[iter] <- 0
          length_G20_thrs[iter] <- NA
        }else{
          hit_array <- array(NA,dim=c(length(tmp),length(true_Y)))
          for (h in 1:length(tmp)){
            for (l in 1:length(true_Y)){
              hit_array[h,l] <- true_Y[l] >= (tmp[[h]][1]-eta) & true_Y[l] <= (tmp[[h]][2]+eta)
            }
          }
          # Hit ratio
          if (any(rowSums(hit_array)/length(true_Y) == 1)) {
            ratio_G20_thrs[iter] <- 1 
          }else {
            ratio_G20_thrs[iter] <-max(rowSums(hit_array)/length(true_Y))
          } 
          
          # Length
          if (any(hit_array!=0)){
            num <- which(rowSums(hit_array)!=0)
            for (m in c(num)){
              length_G20_thrs[iter] <- length_G20_thrs[iter] + tmp[[m]][2] - tmp[[m]][[1]] +1
            }
          }else{
            length_G20_thrs[iter] <- 0
          }
        }
      }
      ratio_G20_thrs <- mean(ratio_G20_thrs)
      length_G20_thrs <-mean(length_G20_thrs,na.rm=TRUE)
      
      # ptwise (G=20)
      ratio_G20_ptwise <- rep(NA,max_iter)
      length_G20_ptwise <- rep(0,max_iter)
      for (iter in 1:max_iter){
        tmp <- list_iter[[iter]]$hotspot_G20_ptwise
        if (is.null(tmp) || (is.list(tmp) && length(tmp) == 0)){
          ratio_G20_ptwise[iter] <- 0
          length_G20_ptwise[iter] <- NA
        }else if (is.null(tmp[[1]]) || length(tmp[[1]]) == 0 || (length(tmp[[1]]) == 1 && is.na(tmp[[1]]))){
          ratio_G20_ptwise[iter] <- 0
          length_G20_ptwise[iter] <- NA
        }else{
          hit_array <- array(NA,dim=c(length(tmp),length(true_Y)))
          for (h in 1:length(tmp)){
            for (l in 1:length(true_Y)){
              hit_array[h,l] <- true_Y[l] >= (tmp[[h]][1]-eta) & true_Y[l] <= (tmp[[h]][2]+eta)
            }
          }
          # Hit ratio
          if (any(rowSums(hit_array)/length(true_Y) == 1)) {
            ratio_G20_ptwise[iter] <- 1 
          }else {
            ratio_G20_ptwise[iter] <-max(rowSums(hit_array)/length(true_Y))
          } 
          
          # Length
          if (any(hit_array!=0)){
            num <- which(rowSums(hit_array)!=0)
            for (m in c(num)){
              length_G20_ptwise[iter] <- length_G20_ptwise[iter] + tmp[[m]][2] - tmp[[m]][[1]] +1
            }
          }else{
            length_G20_ptwise[iter] <- 0
          }
        }
      }
      ratio_G20_ptwise <- mean(ratio_G20_ptwise)
      length_G20_ptwise <- mean(length_G20_ptwise,na.rm=TRUE)
      
      
      # ------------------------------------ #
      # Thrs (G=40)
      ratio_G40_thrs <- rep(NA,max_iter)
      length_G40_thrs <- rep(0,max_iter)
      for (iter in 1:max_iter){
        tmp <- list_iter[[iter]]$hotspot_G40_thres
        if (is.null(tmp) || (is.list(tmp) && length(tmp) == 0)){
          ratio_G40_thrs[iter] <- 0
          length_G40_thrs[iter] <- NA
        }else if (is.null(tmp[[1]]) || length(tmp[[1]]) == 0 || (length(tmp[[1]]) == 1 && is.na(tmp[[1]]))){
          ratio_G40_thrs[iter] <- 0
          length_G40_thrs[iter] <- NA
        }else{
          hit_array <- array(NA,dim=c(length(tmp),length(true_Y)))
          for (h in 1:length(tmp)){
            for (l in 1:length(true_Y)){
              hit_array[h,l] <- true_Y[l] >= (tmp[[h]][1]-eta) & true_Y[l] <= (tmp[[h]][2]+eta)
            }
          }
          # Hit ratio
          if (any(rowSums(hit_array)/length(true_Y) == 1)) {
            ratio_G40_thrs[iter] <- 1 
          }else {
            ratio_G40_thrs[iter] <-max(rowSums(hit_array)/length(true_Y))
          } 
          
          # Length
          if (any(hit_array!=0)){
            num <- which(rowSums(hit_array)!=0)
            for (m in c(num)){
              length_G40_thrs[iter] <- length_G40_thrs[iter] + tmp[[m]][2] - tmp[[m]][[1]] +1
            }
          }else{
            length_G40_thrs[iter] <- 0
          }
        }
      }
      ratio_G40_thrs <- mean(ratio_G40_thrs)
      length_G40_thrs <- mean(length_G40_thrs,na.rm=TRUE)
      
      # ptwise (G=40)
      ratio_G40_ptwise <- rep(NA,max_iter)
      length_G40_ptwise <- rep(0,max_iter)
      for (iter in 1:max_iter){
        tmp <- list_iter[[iter]]$hotspot_G40_ptwise
        if (is.null(tmp) || (is.list(tmp) && length(tmp) == 0)){
          ratio_G40_ptwise[iter] <- 0
          length_G40_ptwise[iter] <- NA
        }else if (is.null(tmp[[1]]) || length(tmp[[1]]) == 0 || (length(tmp[[1]]) == 1 && is.na(tmp[[1]]))){
          ratio_G40_ptwise[iter] <- 0
          length_G40_ptwise[iter] <- NA
        }else{
          hit_array <- array(NA,dim=c(length(tmp),length(true_Y)))
          for (h in 1:length(tmp)){
            for (l in 1:length(true_Y)){
              hit_array[h,l] <- true_Y[l] >= (tmp[[h]][1]-eta) & true_Y[l] <= (tmp[[h]][2]+eta)
            }
          }
          # Hit ratio
          if (any(rowSums(hit_array)/length(true_Y) == 1)) {
            ratio_G40_ptwise[iter] <- 1 
          }else {
            ratio_G40_ptwise[iter] <-max(rowSums(hit_array)/length(true_Y))
          } 
          
          # Length
          if (any(hit_array!=0)){
            num <- which(rowSums(hit_array)!=0)
            for (m in c(num)){
              length_G40_ptwise[iter] <- length_G40_ptwise[iter] + tmp[[m]][2] - tmp[[m]][[1]] +1
            }
          }else{
            length_G40_ptwise[iter] <- 0
          }
        }
      }
      ratio_G40_ptwise <- mean(ratio_G40_ptwise)
      length_G40_ptwise <-mean(length_G40_ptwise,na.rm=TRUE)
      
      
      # ------------------------------------ #
      # display:
      cat("num of Y jumps",i,"Case",j,"delta",delta,"produced.\n")
      cat("hit-miss:\n")
      cat("Thrs (G=20),Thrs (G=40),CI (G=20),CI (G=40)\n")
      print(round(c(ratio_G20_thrs,ratio_G40_thrs,
                    ratio_G20_ptwise,ratio_G40_ptwise),3))
      cat("length:\n")
      cat("Thrs (G=20),Thrs (G=40),CI (G=20),CI (G=40)\n")
      print(round(c(length_G20_thrs,length_G40_thrs,
                    length_G20_ptwise,length_G40_ptwise),3))
      cat("#---------------------------------------------- #\n")
      
    }
  }
}
