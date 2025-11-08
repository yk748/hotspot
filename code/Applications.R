##################################################################
# Data Application:
##################################################################
rm(list=ls())
source("Library.R")

library("ggplot2")
library("tidyr")
library("RColorBrewer")
library("ggpubr")
library("cowplot")
library("patchwork")

# ---------------------------------------------- #
load("hotspot_4subjects_data.RData")
set.seed(1234)
result_list <- list()
cnt <- 0
for (i in 1:4){
  subj <- df_list[[i]]$subj
  dat <- df_list[[i]]$dat
  n <- df_list[[i]]$n
  stress <- df_list[[i]]$dat$stress
  stress <- conv_smooth(stress,n,scale=1:5)
  
  # ---------------------------------------------- #
  for (j in 1:3){
    
    if (j == 1){
      passive <- dat$step_count 
      passive_label <- "step counts"
    }else if (j == 2){
      passive <- dat$conversation_percent
      passive_label <- "conversation_percent"
    }else if (j == 3){
      passive <- dat$sleep_duration
      passive_label <- "sleep_duration"
    }
    
    for (k in 1:2){
      if (k == 1){
        G <- 20
      }else if (k == 2){
        G <- 40
      }
      cnt <- cnt + 1
      
      # ---------------------------------------------- #
      MOSUM <- bivariate_mosum(passive,stress,n=n,Gl=G,Gr=G,alpha=0.05,eta=0.2)
      
      # ---------------------------------------------- #
      # Y:
      if (!is.null(MOSUM$Y_both$cpt_eta)){
        CI_Y <- joint_mosum_ci(stress,refinement=MOSUM$Y_both,
                               method="eta",n=n,
                               Gl=G,Gr=G,N=1000,alpha=0.05,epsilon=0.2,eta=0.2)
      }else{
        CI_Y <- NULL
      }
      
      # ---------------------------------------------- #
      # Y & X:
      if (!is.null(MOSUM$YX$cpt_eta)){
        CI_YX <- bivariate_mosum_ci(stress,passive,
                                    bi_refinement=MOSUM$YX,
                                    moment="mean_mean",
                                    method="eta",n=n,
                                    Gl=G,Gr=G,N=1000,alpha=0.05,epsilon=0.2,eta=0.2)
      }else{
        CI_YX <- NULL
      }
      
      # Y & X2:
      if (!is.null(MOSUM$YX2$cpt_eta)){
        CI_YX2 <- bivariate_mosum_ci(stress,passive,
                                     bi_refinement=MOSUM$YX2,
                                     moment="mean_var",
                                     method="eta",n=n,
                                     Gl=G,Gr=G,N=1000,alpha=0.05,epsilon=0.2,eta=0.2)
      }else{
        CI_YX2 <- NULL
      }
      
      # Y2 & X:
      if (!is.null(MOSUM$Y2X$cpt_eta)){
        CI_Y2X <- bivariate_mosum_ci(stress,passive,
                                     bi_refinement=MOSUM$Y2X,
                                     moment="var_mean",
                                     method="eta",n=n,
                                     Gl=G,Gr=G,N=1000,alpha=0.05,epsilon=0.2,eta=0.2)
      }else{
        CI_Y2X <- NULL
      }
      
      # Y2 & X2:
      if (!is.null(MOSUM$Y2X2$cpt_eta)){
        CI_Y2X2 <- bivariate_mosum_ci(stress,passive,
                                      bi_refinement=MOSUM$Y2X2,
                                      moment="var_var",
                                      method="eta",n=n,
                                      Gl=G,Gr=G,N=1000,alpha=0.05,epsilon=0.2,eta=0.2)
      }else{
        CI_Y2X2 <- NULL
      }
      
      # ---------------------------------------------- #
      # Hotspots:
      hotspot_ci <-get_union_of_intersects(c(CI_YX$ptwise,
                                                 CI_YX2$ptwise,
                                                 CI_Y2X$ptwise,
                                                 CI_Y2X2$ptwise),
                                               CI_Y$ptwise,type="ptwise")
      
      hotspot_thrs <- get_intervals(sort(intersect(unique(c(MOSUM$YX$idx_cand,
                                                            MOSUM$YX2$idx_cand,
                                                            MOSUM$Y2X$idx_cand,
                                                            MOSUM$Y2X2$idx_cand)),
                                                   MOSUM$Y_both$idx_cand)))
      # ------------------------------------ #
      # output:            
      output <- list()
      output$subj <- subj
      output$passive <- passive_label
      output$G <- G
      output$MOSUM <- MOSUM
      output$CI_Y <- CI_Y
      output$CI_YX <- CI_YX
      output$CI_YX2 <- CI_YX2
      output$CI_Y2X <- CI_Y2X
      output$CI_Y2X2 <- CI_Y2X2
      output$hotspot_ci <- hotspot_ci
      output$hotspot_thrs <- hotspot_thrs
      result_list[[cnt]] <- output
      
      cat(subj,"\n")
      cat(passive_label,"\n")
      cat(cnt,"th work done.\n")
    }
  }
}
save(result_list,file="data_appl.RData")

#---------------------------------------------- #

load("data_appl.RData")
for (iter in 1:24){
  # ------------------------------------ #
  # display:
  cat("Subject:",result_list[[iter]]$subj,"\n")
  cat("Passive sensing:",result_list[[iter]]$passive,"\n")
  cat("G:",result_list[[iter]]$G,"\n")
  cat("Change points:\n")
  cat("Y:",result_list[[iter]]$MOSUM$Y_both$cpt_eta,"\n")
  cat("YX:",result_list[[iter]]$MOSUM$YX$cpt_eta,"\n")
  cat("YX2:",result_list[[iter]]$MOSUM$YX2$cpt_eta,"\n")
  cat("Y2X:",result_list[[iter]]$MOSUM$Y2X$cpt_eta,"\n")
  cat("Y2X2:",result_list[[iter]]$MOSUM$Y2X2$cpt_eta,"\n")
  cat("Hotspots:\n")
  cat("Thrs\n")
  print(lapply(as.list(result_list[[iter]]$hotspot_thrs),round,3))
  cat("CI\n")
  print(lapply(result_list[[iter]]$hotspot_ci,round,3))
  cat("#---------------------------------------------- #\n")
}


#---------------------------------------------- #
load("hotspot_4subjects_data.RData")
load("data_appl.RData")
for (i in 1:4){
  
  subj <- df_list[[i]]$subj
  dat <- df_list[[i]]$dat
  n <- df_list[[i]]$n
  
  df <- data.frame(time=1:n,
                   value = c(dat$step_count,dat$conversation_percent,dat$sleep_duration),
                   variable = c(rep("SC",n),rep("CP",n),rep("SD",n)))
  stress <- cut(df_list[[i]]$dat$stress, breaks=5,labels=1:5)
  stress <- as.numeric(as.character(stress))
  df <- rbind(df, data.frame(time=1:n,value=stress,variable="Stress"))
  
  if (i == 1){
    cnt1 <- 2; cnt2 <- 4; cnt3 <- 6
  }else if (i == 2){
    cnt1 <- 7; cnt2 <- 9; cnt3 <- 11
  }else if (i == 3){
    cnt1 <- 14; cnt2 <- 16; cnt3 <- 18
  }else if (i == 4){
    cnt1 <- 20; cnt2 <- 22; cnt3 <- 24
  }
  
  # --------------------------------------------------------- #
  hotspot_thrs1 <- data.frame(do.call(rbind, result_list[[cnt1]]$hotspot_thrs))
  hotspot_ci1 <- data.frame(do.call(rbind, result_list[[cnt1]]$hotspot_ptwise))
  names(hotspot_thrs1) <- names(hotspot_ptwise1) <- c("xmin","xmax")
  hotspot_thrs1$type <- rep("Thrs",dim(hotspot_thrs1)[1])
  hotspot_ci1$type <- rep("CI",dim(hotspot_ci1)[1])
  hotspot1 <- rbind(hotspot_thrs1,hotspot_ci1)
  hotspot1$type <- factor(hotspot1$type,levels=c("Thrs","CI"),labels=c("Thrs","CI"))
  
  hotspot_SC1_jittered <- hotspot1
  hotspot_SC1_jittered$xmin <- hotspot1$xmin + runif(nrow(hotspot1), -0.4,-0.3)
  hotspot_SC1_jittered$xmax <- hotspot1$xmax + runif(nrow(hotspot1), 0.3,0.4)
  
  
  hotspot_thrs2 <- data.frame(do.call(rbind, result_list[[cnt2]]$hotspot_thrs))
  hotspot_ci2 <- data.frame(do.call(rbind, result_list[[cnt2]]$hotspot_ptwise))
  names(hotspot_thrs2) <- names(hotspot_ptwise2) <- c("xmin","xmax")
  hotspot_thrs2$type <- rep("Thrs",dim(hotspot_thrs2)[1])
  hotspot_ci2$type <- rep("CI",dim(hotspot_ci2)[1])
  hotspot2 <- rbind(hotspot_thrs2,hotspot_ci2)
  hotspot2$type <- factor(hotspot2$type,levels=c("Thrs","CI"),labels=c("Thrs","CI"))
  
  hotspot_SC2_jittered <- hotspot2
  hotspot_SC2_jittered$xmin <- hotspot2$xmin + runif(nrow(hotspot2), -0.4,-0.3)
  hotspot_SC2_jittered$xmax <- hotspot2$xmax + runif(nrow(hotspot2), 0.3,0.4)
  
  
  hotspot_thrs3 <- data.frame(do.call(rbind, result_list[[cnt3]]$hotspot_thrs))
  hotspot_ci3 <- data.frame(do.call(rbind, result_list[[cnt3]]$hotspot_ptwise))
  names(hotspot_thrs3) <- names(hotspot_ptwise3) <- c("xmin","xmax")
  hotspot_thrs3$type <- rep("Thrs",dim(hotspot_thrs3)[1])
  hotspot_ci3$type <- rep("CI",dim(hotspot_ci3)[1])
  hotspot3 <- rbind(hotspot_thrs3,hotspot_ci3)
  hotspot3$type <- factor(hotspot3$type,levels=c("Thrs","CI"),labels=c("Thrs","CI"))
  
  hotspot_SC3_jittered <- hotspot3
  hotspot_SC3_jittered$xmin <- hotspot3$xmin + runif(nrow(hotspot3), -0.4,-0.3)
  hotspot_SC3_jittered$xmax <- hotspot3$xmax + runif(nrow(hotspot3), 0.3,0.4)
  
  # --------------------------------------------------------- #
  pl_passive_SC <- ggplot(df[df$variable=="SC",], 
                          aes(x=time, y=value)) +
    geom_line(color = "purple") +
    geom_rect(data = hotspot_SC_jittered,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                  fill = type, alpha = type),
              inherit.aes = FALSE) +
    scale_fill_manual(name = "Hotspots",
                      values = c("Thrs"="red","CI"="blue")) +
    scale_alpha_manual(values = c("Thrs"=0.4,"CI"=0.2)) +
    guides(alpha = "none",
           fill = guide_legend(override.aes = list(
             alpha = c(0.4, 0.2)))) +
    scale_x_continuous(name="Time", breaks = seq(0, n, by = 10)) +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom")
  
  pl_stress_SC <- ggplot(df[df$variable=="Stress",], aes(x=time, y =value)) +
    geom_step(color="black") +
    geom_rect(data = hotspot_SC_jittered,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                  fill = type, alpha = type),
              inherit.aes = FALSE) +
    scale_fill_manual(name = "Hotspots",
                      values = c("Thrs"="red","CI"="blue")) +
    scale_alpha_manual(values = c("Thrs"=0.4,"CI"=0.2)) +
    guides(alpha = "none",
           fill = guide_legend(override.aes = list(
             alpha = c(0.4, 0.2)))) +
    scale_x_continuous(name="Time", breaks = seq(0, n, by = 10)) +
    scale_y_continuous(breaks = 1:5, limits = c(1,5)) +
    labs(x = "Time", y = "Stress") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom")
  
  
  pl_passive_CP <- ggplot(df[df$variable=="CP",], 
                          aes(x=time, y=value)) +
    geom_line(color = "orange") +
    geom_rect(data = hotspot_CP_jittered,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                  fill = type, alpha = type),
              inherit.aes = FALSE) +
    scale_fill_manual(name = "Hotspots",
                      values = c("Thrs"="red","CI"="blue")) +
    scale_alpha_manual(values = c("Thrs"=0.4,"CI"=0.2)) +
    guides(alpha = "none",
           fill = guide_legend(override.aes = list(
             alpha = c(0.4, 0.2)))) +
    scale_x_continuous(name="Time", breaks = seq(0, n, by = 10)) +
    labs(x = "", y = "Conversation Percentage", title = "") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom")
  
  pl_stress_CP <- ggplot(df[df$variable=="Stress",], aes(x=time, y =value)) +
    geom_step(color="black") +
    geom_rect(data = hotspot_CP_jittered,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                  fill = type, alpha = type),
              inherit.aes = FALSE) +
    scale_fill_manual(name = "Hotspots",
                      values = c("Thrs"="red","CI"="blue")) +
    scale_alpha_manual(values = c("Thrs"=0.4,"CI"=0.2)) +
    guides(alpha = "none",
           fill = guide_legend(override.aes = list(
             alpha = c(0.4, 0.2)))) +
    scale_x_continuous(name="Time", breaks = seq(0, n, by = 10)) +
    scale_y_continuous(breaks = 1:5, limits = c(1,5)) +
    labs(x = "Time", y = "Stress") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom")
  
  
  pl_passive_SD <- ggplot(df[df$variable=="SD",], 
                          aes(x=time, y=value)) +
    geom_line(color = "navy") +
    geom_rect(data = hotspot_SD_jittered,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                  fill = type, alpha = type),
              inherit.aes = FALSE) +
    scale_fill_manual(name = "Hotspots",
                      values = c("Thrs"="red","CI"="blue")) +
    scale_alpha_manual(values = c("Thrs"=0.4,"CI"=0.2)) +
    guides(alpha = "none",
           fill = guide_legend(override.aes = list(
             alpha = c(0.4, 0.2)))) +
    scale_x_continuous(name="Time", breaks = seq(0, n, by = 10)) +
    labs(x = "", y = "Sleep duration", title = "") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom")
  
  pl_stress_SD <- ggplot(df[df$variable=="Stress",], aes(x=time, y =value)) +
    geom_step(color="black") +
    geom_rect(data = hotspot_SD_jittered,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                  fill = type, alpha = type),
              inherit.aes = FALSE) +
    scale_fill_manual(name = "Hotspots",
                      values = c("Thrs"="red","CI"="blue")) +
    scale_alpha_manual(values = c("Thrs"=0.4,"CI"=0.2)) +
    guides(alpha = "none",
           fill = guide_legend(override.aes = list(
             alpha = c(0.4, 0.2)))) +
    scale_x_continuous(name="Time", breaks = seq(0, n, by = 10)) +
    scale_y_continuous(breaks = 1:5, limits = c(1,5)) +
    labs(x = "Time", y = "Stress") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom")
  
  if (i == 1){
    pl_passive_SC <- pl_passive_SC + labs(x = "", y = "Step Counts", title = "sternrelief15")
  }else if (i == 2){
    pl_passive_SC <- pl_passive_SC + labs(x = "", y = "Step Counts", title = "sternrelief16")
  }else if (i == 3){
    pl_passive_SC <- pl_passive_SC + labs(x = "", y = "Step Counts", title = "sternrelief5")
  }else if (i == 4){
    pl_passive_SC <- pl_passive_SC + labs(x = "", y = "Step Counts", title = "sternrelief6")
  }
  
  # --------------------------------------------------------- #
  combined_plot <- ((pl_passive_SC|pl_stress_SC)/
      (pl_passive_CP|pl_stress_CP)/
      (pl_passive_SD|pl_stress_SD)) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  pdf(paste0("hotspot_", subj, ".pdf"), width = 13, height = 11)
  print(combined_plot)
  dev.off()
}


