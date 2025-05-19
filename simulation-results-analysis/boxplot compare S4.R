library(ggplot2)
library(tidyr)

results_path <- "simulation-output/proposed_S4/"
setwd(results_path)

#1 p10 rho0.3 N5000
tran_g_dat=NULL
a_dat=NULL
b_dat=NULL
mu_dat=NULL
sigma_dat=NULL
cor_all_N5000_r5_p10_J15=NULL

N=20000
p=20
J=30
rho=2
for (seed in 80:89){
  if(file.exists(paste0('seed',seed,'_tran_gamma_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))){
    temp_p10_rho0.3_N5000 = read.csv(paste0('seed',seed,'_tran_gamma_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    tran_g_dat = rbind(tran_g_dat, temp_p10_rho0.3_N5000)
    
    temp_p10_rho0.3_N5000 = read.csv(paste0('seed',seed,'_alpha_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    a_dat = rbind(a_dat, temp_p10_rho0.3_N5000)
    
    temp_p10_rho0.3_N5000 = read.csv(paste0('seed',seed,'_tran_d_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    b_dat = rbind(b_dat, temp_p10_rho0.3_N5000)
    
    temp_p10_rho0.3_N5000 = read.csv(paste0('seed',seed,'_tran_mu_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    mu_dat = rbind(mu_dat, temp_p10_rho0.3_N5000)
    
    temp_p10_rho0.3_N5000 = read.csv(paste0('seed',seed,'_sigma_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    sigma_dat = rbind(sigma_dat, temp_p10_rho0.3_N5000)
    
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    cor_all_N5000_r5_p10_J15 = rbind(cor_all_N5000_r5_p10_J15, temp)
    
    
   
  }
}

results_path <- "simulation-output/baseline_S4/"
setwd(results_path)

#1 p10 rho0.3 N5000
a_dat_bcmk=NULL
b_dat_bcmk=NULL
mu_dat_bcmk=NULL
sigma_dat_bcmk=NULL
cor_all_N5000_r5_p10_J15_bcmk=NULL

N=20000
p=20
J=30
rho=2
for (seed in 80:89){
  if(file.exists(paste0('seed',seed,'_alpha_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    cor_all_N5000_r5_p10_J15_bcmk = rbind(cor_all_N5000_r5_p10_J15_bcmk, temp)
    
    
    
  }
}


###########################################################################
##################### RMSD method #########################################
###########################################################################
results_path <- "simulation-output/rmsd_S4/"
setwd(results_path)

N=20000
rho=2
p=20
J=30
cor_all_N5000_r5_p10_J15_rmsd_0.05 = NULL
discovery_N5000_r5_p10_J15_rmsd_0.05 = NULL


for (seed in 80:89){
  RMSD_threshold = 0.05
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    cor_all_N5000_r5_p10_J15_rmsd_0.05 = rbind(cor_all_N5000_r5_p10_J15_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_discovery_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    discovery_N5000_r5_p10_J15_rmsd_0.05 = rbind(discovery_N5000_r5_p10_J15_rmsd_0.05, temp)
  }
}

apply(discovery_N5000_r5_p10_J15_rmsd_0.05, 2, mean)
mean(discovery_N5000_r5_p10_J15_rmsd_0.05[,2]/discovery_N5000_r5_p10_J15_rmsd_0.05[,1])

N=20000
rho=2
p=20
J=30
cor_all_N5000_r5_p10_J15_rmsd_0.1 = NULL
discovery_N5000_r5_p10_J15_rmsd_0.1 = NULL
for (seed in 80:89){
  RMSD_threshold = 0.1
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    cor_all_N5000_r5_p10_J15_rmsd_0.1 = rbind(cor_all_N5000_r5_p10_J15_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_discovery_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    discovery_N5000_r5_p10_J15_rmsd_0.1 = rbind(discovery_N5000_r5_p10_J15_rmsd_0.1, temp)
  }
}

apply(discovery_N5000_r5_p10_J15_rmsd_0.1, 2, mean)
mean(discovery_N5000_r5_p10_J15_rmsd_0.1[,2]/discovery_N5000_r5_p10_J15_rmsd_0.1[,1])


N=20000
rho=2
p=20
J=30
cor_all_N5000_r5_p10_J15_rmsd_0.15 = NULL
discovery_N5000_r5_p10_J15_rmsd_0.15 = NULL
for (seed in 80:89){
  RMSD_threshold = 0.15
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    cor_all_N5000_r5_p10_J15_rmsd_0.15 = rbind(cor_all_N5000_r5_p10_J15_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_discovery_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    discovery_N5000_r5_p10_J15_rmsd_0.15 = rbind(discovery_N5000_r5_p10_J15_rmsd_0.15, temp)
  }
}

apply(discovery_N5000_r5_p10_J15_rmsd_0.15, 2, mean)
mean(discovery_N5000_r5_p10_J15_rmsd_0.15[,2]/discovery_N5000_r5_p10_J15_rmsd_0.15[,1])


g=matrix(0, J, p)

set.seed(4)

for(j in 1:J){
  inx = sample(1:p, p/2-1, replace=FALSE)
  g[j, inx] = runif(p/2-1, rho-1, rho+1)
}
nonzero_list = apply(g!=0, 2, sum)
exceeding_cols = which(nonzero_list - ceiling(J/2) + 1 >= 0)
rep_times = nonzero_list[exceeding_cols] - ceiling(J/2) + 1
additional_cols = rep(exceeding_cols, rep_times)

for(c in additional_cols){
  nonzero_rows = which(g[, c]!=0)
  g[sample(nonzero_rows, 1), c] =0
}


true_0 = sum(g==0)
false_0 = sum(g!=0)

## RMSD 0.05
pred_0 = mean(discovery_N5000_r5_p10_J15_rmsd_0.05[,1])
pred_1 = J*p - pred_0

FN = mean(discovery_N5000_r5_p10_J15_rmsd_0.05[,2])
FNR = FN/false_0
TN = pred_0 - FN
FP = true_0 - TN
FPR = FP/true_0
round(c(pred_0, FPR, FNR ), 2)

## RMSD 0.10
pred_0 = mean(discovery_N5000_r5_p10_J15_rmsd_0.1[,1])
pred_1 = J*p - pred_0

FN = mean(discovery_N5000_r5_p10_J15_rmsd_0.1[,2])
FNR = FN/false_0
TN = pred_0 - FN
FP = true_0 - TN
FPR = FP/true_0
round(c(pred_0, FPR, FNR ), 2)

## RMSD 0.15
pred_0 = mean(discovery_N5000_r5_p10_J15_rmsd_0.15[,1])
pred_1 = J*p - pred_0

FN = mean(discovery_N5000_r5_p10_J15_rmsd_0.15[,2])
FNR = FN/false_0
TN = pred_0 - FN
FP = true_0 - TN
FPR = FP/true_0
round(c(pred_0, FPR, FNR ), 2)


data_list <- list(
  "RMSD 0.05" = cor_all_N5000_r5_p10_J15_rmsd_0.05[,1],
  "RMSD 0.10" = cor_all_N5000_r5_p10_J15_rmsd_0.1[,1],
  "RMSD 0.15" = cor_all_N5000_r5_p10_J15_rmsd_0.15[,1],
  "Baseline" = cor_all_N5000_r5_p10_J15_bcmk[,1],
  "Proposed" = cor_all_N5000_r5_p10_J15[,1]
)

# Convert to data frame, handling unequal lengths
plot_data <- data.frame(
  Method = rep(names(data_list), sapply(data_list, length)),
  Value = unlist(data_list)
)

# Create ordered factor for proper x-axis ordering
plot_data$Method <- factor(plot_data$Method,
                           levels = c("RMSD 0.05", "RMSD 0.10", "RMSD 0.15", "Baseline", "Proposed"))
# Create the plot
ggplot(plot_data, aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(width = 0.7, alpha = 0.8, outlier.size = 1.5, outlier.alpha = 0.6) +
  scale_fill_manual(values = c("RMSD 0.05" = "#6baed6",  # Light blue
                               "RMSD 0.10" = "#4292c6",    # Medium blue
                               "RMSD 0.15" = "#2171b5",   # Dark blue
                               "Baseline" = "#59A14F",     # Neutral gray
                               "Proposed" = "#d95f02"     # Vibrant orange (stands out)
                               )) +
  labs(title = "Setting 4 - N = 20000, J = 30, p = 20",
       y = "Kendall Rank Correlation") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 15)),
    axis.title.x =  element_blank(),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(1,1,1,1), "cm"),
    panel.border = element_rect(colour = "gray90", fill = NA, linewidth = 0.5)
  ) +
  ylim(-0.5, 1) 
