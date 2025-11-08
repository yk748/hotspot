##################################################################
# Illustration of the algorithm:
##################################################################
rm(list=ls())
source("Library.R")

library("mvtnorm")
library("stringr")

library("ggplot2")
library("tidyr")
library("RColorBrewer")
library("ggpubr")
library("cowplot")
library("patchwork")

# ---------------------------------------- #
# Example 1:
# ---------------------------------------- #
n <- 100
set.seed(123)
yt <- rep(NA,n)
for (t in 1:n){
  prob1 <- c(0.6,0.3,0.1,0.0,0.0)
  prob2 <- c(0.0,0.0,0.1,0.3,0.6)
  
  if (t < n/2){
    yt[t] <- sample(c(1:5),size=1,prob=prob1,replace=TRUE)
  }else{
    yt[t] <- sample(c(1:5),size=1,prob=prob2,replace=TRUE)
  }
}
zt <- conv_smooth(yt,n,scale=1:5)

xt <- rep(NA,n)
start <- c(1,46)
end <- c(45,n)
levels_vec <- runif(2, min=-2, max=2)
sigma_vec <- runif(2, min=0.4, max=0.8)
for (k in 1:2) {
  xt[start[k]:end[k]] <- levels_vec[k]
  xt[start[k]:end[k]] <- xt[start[k]:end[k]] + rnorm(end[k]-start[k]+1, mean=0, sd=sigma_vec[k])
}


# ---------------------------------------- #
# Plot:
df_series <- data.frame(Time = 1:n,
                        Stress = yt,
                        Continuous_Stress = zt,
                        Passive_Sensing = xt)
df_series <- pivot_longer(df_series,
                          cols = c("Stress", "Continuous_Stress", "Passive_Sensing"),
                          names_to = "Series",
                          values_to = "Value")
df_series$Series <- factor(df_series$Series,
                           levels = c("Stress", "Continuous_Stress", "Passive_Sensing"))


plot_Y <- ggplot(data.frame(Time = 1:n, Value = yt), aes(x = Time, y = Value)) +
  geom_step(color = "black") +
  geom_vline(xintercept = 50 + runif(1,-0.2,0.2),
             linetype = 1, color = "red") +
  scale_x_continuous(breaks = seq(0, n, by = 10)) +
  scale_y_continuous(breaks = 1:5) +
  labs(y = "Stress", x = NULL) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())
plot_Z <- ggplot(data.frame(Time = 1:n, Value = zt), aes(x = Time, y = Value)) +
  geom_line(color = "black") +
  geom_vline(xintercept = 50 + runif(1,-0.2,0.2),
             linetype = 1, color = "red") +
  scale_x_continuous(breaks = seq(0, n, by = 10)) +
  labs(y = "Continuous Stress", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())
plot_X <- ggplot(data.frame(Time = 1:n, Value = xt), aes(x = Time, y = Value)) +
  geom_line(color = "black") +
  geom_vline(xintercept = 45 + runif(1,-0.2,0.2),
             linetype = 1, color = "red") +
  scale_x_continuous(breaks = seq(0, n, by = 10)) +
  labs(y = "Passive Sensing", x = "Time") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())


pdf("example_samplepath1.pdf", width = 13, height = 11)
plot_Y/plot_Z/plot_X
dev.off()

# ---------------------------------------- #
# Approach:
Bi_MOSUM <- bivariate_mosum(xt,zt,n=n,Gl=20,Gr=20,alpha=0.05,epsilon=0.2,eta=0.2)

# Y both + X both:
dist_Y <- Bi_MOSUM$Y_track$Dist
dist_X <- Bi_MOSUM$X_track$Dist

# All boths:
dist_YX <- Bi_MOSUM$bi_track$Dist_YX
dist_YX2 <- Bi_MOSUM$bi_track$Dist_YX2

time_list_full <- list()
time_list_full$Y_both <- if (!is.null(Bi_MOSUM$Y_both$cpt_eta)) Bi_MOSUM$Y_both$idx_cand else integer(0)
time_list_full$X_both <- if (!is.null(Bi_MOSUM$X_both$cpt_eta)) Bi_MOSUM$X_both$idx_cand else integer(0)
time_list_full$YX  <- if (!is.null(Bi_MOSUM$YX$cpt_eta))  Bi_MOSUM$YX$idx_cand  else integer(0)
time_list_full$YX2 <- if (!is.null(Bi_MOSUM$YX2$cpt_eta)) Bi_MOSUM$YX2$idx_cand else integer(0)


CI_Y <- joint_mosum_ci(zt,refinement=Bi_MOSUM$Y_both,
                       method="eta",n,Gl=20,Gr=20,N=1000,
                       alpha=0.05,epsilon=0.2,eta=0.2)
CI_X <- joint_mosum_ci(xt,refinement=Bi_MOSUM$X_both,
                       method="eta",n,Gl=20,Gr=20,N=1000,
                       alpha=0.05,epsilon=0.2,eta=0.2)

CI_YX <- bivariate_mosum_ci(zt,xt,bi_refinement=Bi_MOSUM$YX,
                            moment="mean_mean",
                            method="eta",n,Gl=20,Gr=20,N=1000,
                            alpha=0.05,epsilon=0.2,eta=0.2)
CI_YX2 <- bivariate_mosum_ci(zt,xt,bi_refinement=Bi_MOSUM$YX2,
                             moment="mean_var",
                             method="eta",n,Gl=20,Gr=20,N=1000,
                             alpha=0.05,epsilon=0.2,eta=0.2)


# Hotspots:
hotspot_thres <- get_intervals(sort(unique(intersect(unique(c(Bi_MOSUM$YX$idx_cand,
                                                              Bi_MOSUM$YX2$idx_cand)),
                                                     Bi_MOSUM$Y_both$idx_cand))))
hotspot_ci <- get_union_of_intersects(
  c(CI_YX$ptwise,CI_YX2$ptwise),
  CI_Y$ptwise,type="ptwise")

hotspot_thres
hotspot_ci

# Thresholds:
max_G <- function(g, WE, WV) {
  Tt <- length(WE)
  LE <- (WE[(2*g+1):(n)] - 2*WE[(1+g):(n-g)] + WE[1:(n-2*g)])/sqrt(2*g)
  LV <- (WV[(2*g+1):(n)] - 2*WV[(1+g):(n-g)] + WV[1:(n-2*g)])/sqrt(2*g)
  max(sqrt(LE^2 + LV^2))
}
sim_global_max <- function(n, G) {
  WE <- cumsum(rnorm(n))
  WV <- cumsum(rnorm(n))
  maxL_H <- vapply(as.matrix(G), FUN = max_G, WE = WE, 
                   WV = WV, numeric(length(1)))
  max(maxL_H)
}

alpha <- 0.05
G <- 20
G_seq <- seq(G,min((n-1)/2,200),10)
Dn <- quantile(replicate(1000, sim_global_max(n=n, G=G_seq)), 1-alpha)
jit_max <- 0.2


# data frame
df_joint <- data.frame(Time = 1:n,
                       Y_both = dist_Y,
                       X_both = dist_X)
df_joint_long <- reshape(df_joint ,
                         varying = c("Y_both", "X_both"),
                         v.names = "Value",
                         times = c("Y both", "X both"),
                         timevar = "Series",
                         direction = "long")
df_bi <- data.frame(Time = 1:n,
                    YX = dist_YX,
                    YX2 = dist_YX2)
df_bi_long <- reshape(df_bi,
                      varying = c("YX", "YX2"),
                      v.names = "Value",
                      times = c("YX", "YX2"),
                      timevar = "Series",
                      direction = "long")


# Confidence intervals
df_ci_Y <- do.call(rbind, lapply(CI_Y$ptwise, function(x) data.frame(xleft = x[1], xright = x[2], Series = "Y both")))
df_ci_X <- do.call(rbind, lapply(CI_X$ptwise, function(x) data.frame(xleft = x[1], xright = x[2], Series = "X both")))
df_ci1 <- rbind(df_ci_Y, df_ci_X)
df_ci_YX <- do.call(rbind, lapply(CI_YX$ptwise, function(x) data.frame(xleft = x[1], xright = x[2], Series = "YX")))
df_ci_YX2 <- do.call(rbind, lapply(CI_YX2$ptwise, function(x) data.frame(xleft = x[1], xright = x[2], Series = "YX2")))
df_ci2 <- rbind(df_ci_YX, df_ci_YX2)


# ---------------------------------------- #
pl_joint <- ggplot(df_joint_long, aes(x = Time, y = Value, 
                                      color = Series, linetype = Series)) +
  geom_rect(data = df_ci1,
            inherit.aes = FALSE,
            aes(xmin = xleft, xmax = xright, ymin = -Inf, ymax = Inf, fill = Series),
            alpha = 0.2) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = Dn, color = "black") +
  geom_vline(xintercept = time_list_full$Y_both + runif(1, -jit_max, jit_max), color = "red", linetype = 1) +
  geom_vline(xintercept = time_list_full$X_both + runif(1, -jit_max, jit_max), color = "blue", linetype = 2) +
  geom_vline(xintercept = Bi_MOSUM$Y_both$cpt_eta + runif(1, -jit_max, jit_max), color = "red", linewidth = 1.2) +
  geom_vline(xintercept = Bi_MOSUM$X_both$cpt_eta + runif(1, -jit_max, jit_max), color = "blue", linewidth = 1.2) +
  scale_color_manual(name="Type",values = c("Y both" = "red", "X both" = "blue")) +
  scale_fill_manual(name="Type",values = c("Y both" = rgb(1, 0, 0, 0.2), 
                                           "X both" = rgb(0, 0, 1, 0.2))) +
  scale_linetype_manual(name="Type",values = c("Y both" = 1, "X both" = 2)) +
  scale_x_continuous(breaks = seq(0, n, by = 10)) +
  labs(y = "Dist", x = "Time") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

pl_bi <- ggplot(df_bi_long, aes(x = Time, y = Value, 
                                color = Series, linetype = Series)) +
  geom_rect(data = df_ci2,
            inherit.aes = FALSE,
            aes(xmin = xleft, xmax = xright, ymin = -Inf, ymax = Inf, fill = Series),
            alpha = 0.2) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = Dn, color = "black") +
  geom_vline(xintercept = time_list_full$YX + runif(1, -jit_max, jit_max), color = "purple", linetype = 3) +
  geom_vline(xintercept = time_list_full$YX2 + runif(1, -jit_max, jit_max), color = "green", linetype = 4) +
  geom_vline(xintercept = Bi_MOSUM$YX$cpt_eta + runif(1, -jit_max, jit_max), color = "purple", linewidth = 1.2) +
  geom_vline(xintercept = Bi_MOSUM$YX2$cpt_eta + runif(1, -jit_max, jit_max), color = "green", linewidth = 1.2) +
  scale_color_manual(name="Type",values = c("YX" = "purple", "YX2" = "green")) +
  scale_fill_manual(name="Type",values = c("YX" = rgb(0.6, 0, 0.8, 0.2), "YX2" = rgb(0, 1, 0, 0.2))) +
  scale_linetype_manual(name="Type",values = c("YX" = 3, "YX2" = 4)) +
  scale_x_continuous(breaks = seq(0, n, by = 10)) +
  labs(y = "Dist", x = "Time") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

pdf("example_result1.pdf", width = 13, height = 11)
pl_joint / pl_bi
dev.off()


# ---------------------------------------- #
# Example 2:
# ---------------------------------------- #
n <- 100
set.seed(123)
yt <- rep(NA,n)
for (t in 1:n){
  prob1 <- c(0.0,1.0,0.0,0.0,0.0)
  prob2 <- c(0.3,0.0,0.4,0.0,0.3)
  prob3 <- c(0.0,0.0,0.0,1.0,0.0)
  
  if (t <= 45){
    yt[t] <- sample(c(1:5),size=1,prob=prob1,replace=TRUE)
  }else if (t > 45 & t <= 65){
    yt[t] <- sample(c(1:5),size=1,prob=prob2,replace=TRUE)
  }else{
    yt[t] <- sample(c(1:5),size=1,prob=prob3,replace=TRUE)
  }
}
zt <- conv_smooth(yt,n,scale=1:5)

xt <- rep(NA,n)
start <- c(1,41,61)
end <- c(40,60,n)
levels_vec <- runif(3, min=-1, max=1)
sigma_vec <- runif(3, min=0.1, max=0.8)
for (k in 1:3) {
  xt[start[k]:end[k]] <- levels_vec[k]
  xt[start[k]:end[k]] <- xt[start[k]:end[k]] + rnorm(end[k]-start[k]+1, mean=0, sd=sigma_vec[k])
}

plot_Y <- ggplot(data.frame(Time = 1:n, Value = yt), aes(x = Time, y = Value)) +
  geom_step(color = "black") +
  geom_vline(xintercept = c(45,65) + runif(2,-0.2,0.2),
             linetype = c(1,2), color = "red") +
  scale_x_continuous(breaks = seq(0, n, by = 10)) +
  scale_y_continuous(breaks = 1:5) +
  labs(y = "Stress", x = NULL) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())
plot_Z <- ggplot(data.frame(Time = 1:n, Value = zt), aes(x = Time, y = Value)) +
  geom_line(color = "black") +
  geom_vline(xintercept = c(45,65) + runif(2,-0.2,0.2),
             linetype = c(1,2), color = "red") +
  scale_x_continuous(breaks = seq(0, n, by = 10)) +
  labs(y = "Continuous Stress", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())
plot_X <- ggplot(data.frame(Time = 1:n, Value = xt), aes(x = Time, y = Value)) +
  geom_line(color = "black") +
  geom_vline(xintercept = c(40,60) + runif(2,-0.2,0.2),
             linetype = c(1,2), color = "red") +
  scale_x_continuous(breaks = seq(0, n, by = 10)) +
  labs(y = "Passive Sensing", x = "Time") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())


pdf("example_samplepath2.pdf", width = 13, height = 11)
plot_Y/plot_Z/plot_X
dev.off()

# ---------------------------------------- #
# Approach:
Bi_MOSUM <- bivariate_mosum(xt,zt,n=n,Gl=20,Gr=20,alpha=0.05,epsilon=0.2,eta=0.2)

# Y both + X both:
dist_Y <- Bi_MOSUM$Y_track$Dist
dist_X <- Bi_MOSUM$X_track$Dist

# All boths:
dist_Y2X <- Bi_MOSUM$bi_track$Dist_YX
dist_Y2X2 <- Bi_MOSUM$bi_track$Dist_YX2

time_list_full <- list()
time_list_full$Y_both <- if (!is.null(Bi_MOSUM$Y_both$cpt_eta)) Bi_MOSUM$Y_both$idx_cand else integer(0)
time_list_full$X_both <- if (!is.null(Bi_MOSUM$X_both$cpt_eta)) Bi_MOSUM$X_both$idx_cand else integer(0)
time_list_full$Y2X  <- if (!is.null(Bi_MOSUM$Y2X$cpt_eta))  Bi_MOSUM$Y2X$idx_cand  else integer(0)
time_list_full$Y2X2 <- if (!is.null(Bi_MOSUM$Y2X2$cpt_eta)) Bi_MOSUM$Y2X2$idx_cand else integer(0)


CI_Y <- joint_mosum_ci(zt,refinement=Bi_MOSUM$Y_both,
                       method="eta",n,Gl=20,Gr=20,N=1000,
                       alpha=0.05,epsilon=0.2,eta=0.2)
CI_X <- joint_mosum_ci(xt,refinement=Bi_MOSUM$X_both,
                       method="eta",n,Gl=20,Gr=20,N=1000,
                       alpha=0.05,epsilon=0.2,eta=0.2)

CI_Y2X <- bivariate_mosum_ci(zt,xt,bi_refinement=Bi_MOSUM$Y2X,
                             moment="var_mean",
                             method="eta",n,Gl=20,Gr=20,N=1000,
                             alpha=0.05,epsilon=0.2,eta=0.2)
CI_Y2X2 <- bivariate_mosum_ci(zt,xt,bi_refinement=Bi_MOSUM$Y2X2,
                              moment="var_var",
                              method="eta",n,Gl=20,Gr=20,N=1000,
                              alpha=0.05,epsilon=0.2,eta=0.2)


# Hotspots:
hotspot_thres <- get_intervals(sort(unique(intersect(unique(c(Bi_MOSUM$Y2X$idx_cand,
                                                              Bi_MOSUM$Y2X2$idx_cand)),
                                                     Bi_MOSUM$Y_both$idx_cand))))
hotspot_ci <- get_union_of_intersects(
  c(CI_Y2X$ptwise,CI_Y2X2$ptwise),
  CI_Y$ptwise,type="ptwise")

hotspot_thres
hotspot_ci


# Thresholds:
max_G <- function(g, WE, WV) {
  Tt <- length(WE)
  LE <- (WE[(2*g+1):(n)] - 2*WE[(1+g):(n-g)] + WE[1:(n-2*g)])/sqrt(2*g)
  LV <- (WV[(2*g+1):(n)] - 2*WV[(1+g):(n-g)] + WV[1:(n-2*g)])/sqrt(2*g)
  max(sqrt(LE^2 + LV^2))
}
sim_global_max <- function(n, G) {
  WE <- cumsum(rnorm(n))
  WV <- cumsum(rnorm(n))
  maxL_H <- vapply(as.matrix(G), FUN = max_G, WE = WE, 
                   WV = WV, numeric(length(1)))
  max(maxL_H)
}

alpha <- 0.05
G <- 20
G_seq <- seq(G,min((n-1)/2,200),10)
Dn <- quantile(replicate(1000, sim_global_max(n=n, G=G_seq)), 1-alpha)
jit_max <- 0.2


# data frame
df_joint <- data.frame(Time = 1:n,
                       Y_both = dist_Y,
                       X_both = dist_X)
df_joint_long <- reshape(df_joint ,
                         varying = c("Y_both", "X_both"),
                         v.names = "Value",
                         times = c("Y both", "X both"),
                         timevar = "Series",
                         direction = "long")
df_bi <- data.frame(Time = 1:n,
                    Y2X = dist_Y2X,
                    Y2X2 = dist_Y2X2)
df_bi_long <- reshape(df_bi,
                      varying = c("Y2X", "Y2X2"),
                      v.names = "Value",
                      times = c("Y2X", "Y2X2"),
                      timevar = "Series",
                      direction = "long")

# Confidence intervals
df_ci_Y <- do.call(rbind, lapply(CI_Y$ptwise, function(x) data.frame(xleft = x[1], xright = x[2], Series = "Y both")))
df_ci_X <- do.call(rbind, lapply(CI_X$ptwise, function(x) data.frame(xleft = x[1], xright = x[2], Series = "X both")))
df_ci1 <- rbind(df_ci_Y, df_ci_X)
df_ci_Y2X <- do.call(rbind, lapply(CI_Y2X$ptwise, function(x) data.frame(xleft = x[1], xright = x[2], Series = "Y2X")))
df_ci_Y2X2 <- do.call(rbind, lapply(CI_Y2X2$ptwise, function(x) data.frame(xleft = x[1], xright = x[2], Series = "Y2X2")))
df_ci2 <- rbind(df_ci_Y2X, df_ci_Y2X2)

# ---------------------------------------- #
pl_joint <- ggplot(df_joint_long, aes(x = Time, y = Value, 
                                      color = Series, linetype = Series)) +
  geom_rect(data = df_ci1,
            inherit.aes = FALSE,
            aes(xmin = xleft, xmax = xright, ymin = -Inf, ymax = Inf, fill = Series),
            alpha = 0.2) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = Dn, color = "black") +
  geom_vline(xintercept = time_list_full$Y_both + runif(1, -jit_max, jit_max), color = "red", linetype = 1) +
  geom_vline(xintercept = time_list_full$X_both + runif(1, -jit_max, jit_max), color = "blue", linetype = 2) +
  geom_vline(xintercept = Bi_MOSUM$Y_both$cpt_eta + runif(1, -jit_max, jit_max), color = "red", linewidth = 1.2) +
  geom_vline(xintercept = Bi_MOSUM$X_both$cpt_eta + runif(1, -jit_max, jit_max), color = "blue", linewidth = 1.2) +
  scale_color_manual(name="Type",values = c("Y both" = "red", "X both" = "blue")) +
  scale_fill_manual(name="Type",values = c("Y both" = rgb(1, 0, 0, 0.2), 
                                           "X both" = rgb(0, 0, 1, 0.2))) +
  scale_linetype_manual(name="Type",values = c("Y both" = 1, "X both" = 2)) +
  scale_x_continuous(breaks = seq(0, n, by = 10)) +
  labs(y = "Dist", x = "Time") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

pl_bi <- ggplot(df_bi_long, aes(x = Time, y = Value, 
                                color = Series, linetype = Series)) +
  geom_rect(data = df_ci2,
            inherit.aes = FALSE,
            aes(xmin = xleft, xmax = xright, ymin = -Inf, ymax = Inf, fill = Series),
            alpha = 0.2) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = Dn, color = "black") +
  geom_vline(xintercept = time_list_full$Y2X + runif(1, -jit_max, jit_max), color = "purple", linetype = 3) +
  geom_vline(xintercept = time_list_full$Y2X2 + runif(1, -jit_max, jit_max), color = "green", linetype = 4) +
  geom_vline(xintercept = Bi_MOSUM$Y2X$cpt_eta + runif(1, -jit_max, jit_max), color = "purple", linewidth = 1.2) +
  geom_vline(xintercept = Bi_MOSUM$Y2X2$cpt_eta + runif(1, -jit_max, jit_max), color = "green", linewidth = 1.2) +
  scale_color_manual(name="Type",values = c("Y2X" = "purple", "Y2X2" = "green")) +
  scale_fill_manual(name="Type",values = c("Y2X" = rgb(0.6, 0, 0.8, 0.2), "Y2X2" = rgb(0, 1, 0, 0.2))) +
  scale_linetype_manual(name="Type",values = c("Y2X" = 3, "Y2X2" = 4)) +
  scale_x_continuous(breaks = seq(0, n, by = 10)) +
  labs(y = "Dist", x = "Time") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())


pdf("example_result2.pdf", width = 13, height = 11)
pl_joint / pl_bi
dev.off()
