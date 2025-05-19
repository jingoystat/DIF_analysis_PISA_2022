library(ggplot2)
library(tidyr)

results_path <- "simulation-output/proposed_S1/"
setwd(results_path)

tran_g_dat=NULL
a_dat=NULL
b_dat=NULL
mu_dat=NULL
sigma_dat=NULL
cor_all_N10000_r5_p10_J15=NULL

N=10000
p=10
J=15
rho=2
for (seed in 80:89){
  if(file.exists(paste0('seed',seed,'_tran_gamma_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))){
    temp_p10_rho0.3_N10000 = read.csv(paste0('seed',seed,'_tran_gamma_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    tran_g_dat = rbind(tran_g_dat, temp_p10_rho0.3_N10000)
    
    temp_p10_rho0.3_N10000 = read.csv(paste0('seed',seed,'_alpha_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    a_dat = rbind(a_dat, temp_p10_rho0.3_N10000)
    
    temp_p10_rho0.3_N10000 = read.csv(paste0('seed',seed,'_tran_d_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    b_dat = rbind(b_dat, temp_p10_rho0.3_N10000)
    
    temp_p10_rho0.3_N10000 = read.csv(paste0('seed',seed,'_tran_mu_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    mu_dat = rbind(mu_dat, temp_p10_rho0.3_N10000)
    
    temp_p10_rho0.3_N10000 = read.csv(paste0('seed',seed,'_sigma_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    sigma_dat = rbind(sigma_dat, temp_p10_rho0.3_N10000)
    
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    cor_all_N10000_r5_p10_J15 = rbind(cor_all_N10000_r5_p10_J15, temp)
    
    
   
  }
}



###########################################################################
##################### Benchmark method ####################################
###########################################################################
results_path <- "simulation-output/baseline_S1/"
setwd(results_path)
a_dat_bcmk=NULL
b_dat_bcmk=NULL
mu_dat_bcmk=NULL
sigma_dat_bcmk=NULL
cor_all_N10000_r5_p10_J15_bcmk=NULL

N=10000
p=10
J=15
rho=2
for (seed in 80:89){
  if(file.exists(paste0('seed',seed,'_alpha_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    cor_all_N10000_r5_p10_J15_bcmk = rbind(cor_all_N10000_r5_p10_J15_bcmk, temp)
    
  }
}


###########################################################################
##################### RMSD method #########################################
###########################################################################
results_path <- "simulation-output/rmsd_S1/"
setwd(results_path)

N=10000
rho=2
p=10
J=15
cor_all_N10000_r5_p10_J15_rmsd_0.05 = NULL
discovery_N10000_r5_p10_J15_rmsd_0.05 = NULL


for (seed in 80:89){
  RMSD_threshold = 0.05
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    cor_all_N10000_r5_p10_J15_rmsd_0.05 = rbind(cor_all_N10000_r5_p10_J15_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_discovery_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    discovery_N10000_r5_p10_J15_rmsd_0.05 = rbind(discovery_N10000_r5_p10_J15_rmsd_0.05, temp)
  }
}

N=10000
rho=2
p=10
J=15
cor_all_N10000_r5_p10_J15_rmsd_0.1 = NULL
discovery_N10000_r5_p10_J15_rmsd_0.1 = NULL
for (seed in 80:89){
  RMSD_threshold = 0.1
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    cor_all_N10000_r5_p10_J15_rmsd_0.1 = rbind(cor_all_N10000_r5_p10_J15_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_discovery_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    discovery_N10000_r5_p10_J15_rmsd_0.1 = rbind(discovery_N10000_r5_p10_J15_rmsd_0.1, temp)
  }
}

apply(discovery_N10000_r5_p10_J15_rmsd_0.1, 2, mean)
mean(discovery_N10000_r5_p10_J15_rmsd_0.1[,2]/discovery_N10000_r5_p10_J15_rmsd_0.1[,1])


N=10000
rho=2
p=10
J=15
cor_all_N10000_r5_p10_J15_rmsd_0.15 = NULL
discovery_N10000_r5_p10_J15_rmsd_0.15 = NULL
for (seed in 80:89){
  RMSD_threshold = 0.15
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    cor_all_N10000_r5_p10_J15_rmsd_0.15 = rbind(cor_all_N10000_r5_p10_J15_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_discovery_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    discovery_N10000_r5_p10_J15_rmsd_0.15 = rbind(discovery_N10000_r5_p10_J15_rmsd_0.15, temp)
  }
}

apply(discovery_N10000_r5_p10_J15_rmsd_0.15, 2, mean)
mean(discovery_N10000_r5_p10_J15_rmsd_0.15[,2]/discovery_N10000_r5_p10_J15_rmsd_0.15[,1])

#### True Gamma #####

J=15
p=10
g=matrix(0, J, p)

set.seed(3)

row_indices_1 <- 1:(J/5)
row_indices_2 <- (J/5+2):(floor(J/2)+J/5) #rep(1:p, times = floor(J/2))  # Each column appears 4 times
row_indices_3 <- c(J/5+1, (floor(J/2)+J/5+1):(floor(J/2)+J/5+1)) #rep(1:p, times = floor(J/2))  # Each column appears 4 times
row_indices_4 <- (floor(J/2)+J/5+2):J

col_indices_1 <- (p/2+1):p
col_indices_2 <- 2:(p/2) 
col_indices_3 <- c(1, (p/2+1):(p/2+2))

l_r1 <- length(row_indices_1)
l_r2 <- length(row_indices_2)
l_r3 <- length(row_indices_3)
l_r4 <- length(row_indices_4)

l_c1 <- length(col_indices_1)
l_c2 <- length(col_indices_2)
l_c3 <- length(col_indices_3)

g[row_indices_1, col_indices_1] <- matrix(runif(l_r1*l_c1, rho-1, rho+1))  # Assign a random value
g[row_indices_2, col_indices_2] <- matrix(runif(l_r2*l_c2, rho-1, rho+1))  # Assign a random value
g[row_indices_3, col_indices_3] <- matrix(runif(l_r3*l_c3, rho-1, rho+1))  # Assign a random value
g[row_indices_4, col_indices_1] <- matrix(runif(l_r4*l_c1, rho-1, rho+1))  # Assign a random value


zero_col=c(6,9,7,7,9,6,9)
for(j in c(row_indices_1, row_indices_4)){
  g[j, zero_col[1]]=0
  zero_col = zero_col[-1]
}


true_0 = sum(g==0)
false_0 = sum(g!=0)

## RMSD 0.05
pred_0 = mean(discovery_N10000_r5_p10_J15_rmsd_0.05[,1])
pred_1 = J*p - pred_0

FN = mean(discovery_N10000_r5_p10_J15_rmsd_0.05[,2])
FNR = FN/false_0
TN = pred_0 - FN
FP = true_0 - TN
FPR = FP/true_0
round(c(pred_0, FPR, FNR ), 2)

## RMSD 0.10
pred_0 = mean(discovery_N10000_r5_p10_J15_rmsd_0.1[,1])
pred_1 = J*p - pred_0

FN = mean(discovery_N10000_r5_p10_J15_rmsd_0.1[,2])
FNR = FN/false_0
TN = pred_0 - FN
FP = true_0 - TN
FPR = FP/true_0
round(c(pred_0, FPR, FNR ), 2)

## RMSD 0.15
pred_0 = mean(discovery_N10000_r5_p10_J15_rmsd_0.15[,1])
pred_1 = J*p - pred_0

FN = mean(discovery_N10000_r5_p10_J15_rmsd_0.15[,2])
FNR = FN/false_0
TN = pred_0 - FN
FP = true_0 - TN
FPR = FP/true_0
round(c(pred_0, FPR, FNR ), 2)


plot_data <- data.frame(
  Value = c(cor_all_N10000_r5_p10_J15_rmsd_0.05[,1],
            cor_all_N10000_r5_p10_J15_rmsd_0.1[,1],
            cor_all_N10000_r5_p10_J15_rmsd_0.15[,1],
            cor_all_N10000_r5_p10_J15_bcmk[,1],
            cor_all_N10000_r5_p10_J15[,1]),
  Method = rep(c("RMSD 0.05", "RMSD 0.10", "RMSD 0.15", "Baseline", "Proposed"),
              each = length(cor_all_N10000_r5_p10_J15_rmsd_0.05[,1]))
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
  labs(title = "Setting 1 - N = 10000, J = 15, p = 10",
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
  ylim(0, 1) 

