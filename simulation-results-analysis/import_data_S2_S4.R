results_path <- "simulation-output/proposed_S2/"
setwd(results_path)

##### n = 10000, J = 15, p = 10 #######

#1 p10 rho0.3 N10000
tran_g_dat=NULL
a_dat=NULL
b_dat=NULL
mu_dat=NULL
sigma_dat=NULL

g_MSE_N10000_r2_p10_J15=NULL
a_MSE_N10000_r2_p10_J15=NULL
b_MSE_N10000_r2_p10_J15=NULL
mu_MSE_N10000_r2_p10_J15=NULL
sigma_MSE_N10000_r2_p10_J15=NULL
pvalgam_N10000_r2_p10_J15=NULL
cor_all_N10000_r2_p10_J15=NULL

N=10000
p=10
J=15
rho=2
for (seed in 80:89){
  if(file.exists(paste0('seed',seed,'_tran_gamma_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    cor_all_N10000_r2_p10_J15 = rbind(cor_all_N10000_r2_p10_J15, temp)
    
    temp = read.csv(paste0('seed',seed,'_g_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    g_MSE_N10000_r2_p10_J15 = rbind(g_MSE_N10000_r2_p10_J15, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    a_MSE_N10000_r2_p10_J15 = rbind(a_MSE_N10000_r2_p10_J15, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    b_MSE_N10000_r2_p10_J15 = rbind(b_MSE_N10000_r2_p10_J15, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    mu_MSE_N10000_r2_p10_J15 = rbind(mu_MSE_N10000_r2_p10_J15, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    sigma_MSE_N10000_r2_p10_J15 = rbind(sigma_MSE_N10000_r2_p10_J15, temp)
    
    # temp = read.csv(paste0('seed',seed,'_pvaluegam_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'), header = F)
    # pvalgam_N10000_r2_p10_J15 = rbind(pvalgam_N10000_r2_p10_J15, temp)
    
  }
}




results_path <- "simulation-output/baseline_S2/"
setwd(results_path)

#1 p10 rho0.3 N10000

a_MSE_N10000_r2_p10_J15_bcmk=NULL
b_MSE_N10000_r2_p10_J15_bcmk=NULL
mu_MSE_N10000_r2_p10_J15_bcmk=NULL
sigma_MSE_N10000_r2_p10_J15_bcmk=NULL

cor_all_N10000_r2_p10_J15_bcmk=NULL

N=10000
p=10
J=15
rho=2
for (seed in 80:89){
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))){
    
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    cor_all_N10000_r2_p10_J15_bcmk = rbind(cor_all_N10000_r2_p10_J15_bcmk, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    a_MSE_N10000_r2_p10_J15_bcmk = rbind(a_MSE_N10000_r2_p10_J15_bcmk, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    b_MSE_N10000_r2_p10_J15_bcmk = rbind(b_MSE_N10000_r2_p10_J15_bcmk, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    mu_MSE_N10000_r2_p10_J15_bcmk = rbind(mu_MSE_N10000_r2_p10_J15_bcmk, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    sigma_MSE_N10000_r2_p10_J15_bcmk = rbind(sigma_MSE_N10000_r2_p10_J15_bcmk, temp)
    
    
  }
}

###########################################################################
##################### RMSD method #########################################
###########################################################################
results_path <- "simulation-output/rmsd_S2"
setwd(results_path)

N=10000
rho=2
p=10
J=15
cor_all_N10000_r2_p10_J15_rmsd_0.05 = NULL
g_MSE_N10000_r2_p10_J15_rmsd_0.05 = NULL
a_MSE_N10000_r2_p10_J15_rmsd_0.05 = NULL
b_MSE_N10000_r2_p10_J15_rmsd_0.05 = NULL
mu_MSE_N10000_r2_p10_J15_rmsd_0.05 = NULL
sigma_MSE_N10000_r2_p10_J15_rmsd_0.05 = NULL

for (seed in 80:89){
  RMSD_threshold = 0.05
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    cor_all_N10000_r2_p10_J15_rmsd_0.05 = rbind(cor_all_N10000_r2_p10_J15_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_g_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    g_MSE_N10000_r2_p10_J15_rmsd_0.05 = rbind(g_MSE_N10000_r2_p10_J15_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    a_MSE_N10000_r2_p10_J15_rmsd_0.05 = rbind(a_MSE_N10000_r2_p10_J15_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    b_MSE_N10000_r2_p10_J15_rmsd_0.05 = rbind(b_MSE_N10000_r2_p10_J15_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    mu_MSE_N10000_r2_p10_J15_rmsd_0.05 = rbind(mu_MSE_N10000_r2_p10_J15_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    sigma_MSE_N10000_r2_p10_J15_rmsd_0.05 = rbind(sigma_MSE_N10000_r2_p10_J15_rmsd_0.05, temp)
    
  }
}

N=10000
rho=2
p=10
J=15
cor_all_N10000_r2_p10_J15_rmsd_0.1 = NULL
g_MSE_N10000_r2_p10_J15_rmsd_0.1 = NULL
a_MSE_N10000_r2_p10_J15_rmsd_0.1 = NULL
b_MSE_N10000_r2_p10_J15_rmsd_0.1 = NULL
mu_MSE_N10000_r2_p10_J15_rmsd_0.1 = NULL
sigma_MSE_N10000_r2_p10_J15_rmsd_0.1 = NULL

for (seed in 80:89){
  RMSD_threshold = 0.1
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    cor_all_N10000_r2_p10_J15_rmsd_0.1 = rbind(cor_all_N10000_r2_p10_J15_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_g_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    g_MSE_N10000_r2_p10_J15_rmsd_0.1 = rbind(g_MSE_N10000_r2_p10_J15_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    a_MSE_N10000_r2_p10_J15_rmsd_0.1 = rbind(a_MSE_N10000_r2_p10_J15_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    b_MSE_N10000_r2_p10_J15_rmsd_0.1 = rbind(b_MSE_N10000_r2_p10_J15_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    mu_MSE_N10000_r2_p10_J15_rmsd_0.1 = rbind(mu_MSE_N10000_r2_p10_J15_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    sigma_MSE_N10000_r2_p10_J15_rmsd_0.1 = rbind(sigma_MSE_N10000_r2_p10_J15_rmsd_0.1, temp)
    
  }
}

N=10000
rho=2
p=10
J=15
cor_all_N10000_r2_p10_J15_rmsd_0.15 = NULL
g_MSE_N10000_r2_p10_J15_rmsd_0.15 = NULL
a_MSE_N10000_r2_p10_J15_rmsd_0.15 = NULL
b_MSE_N10000_r2_p10_J15_rmsd_0.15 = NULL
mu_MSE_N10000_r2_p10_J15_rmsd_0.15 = NULL
sigma_MSE_N10000_r2_p10_J15_rmsd_0.15 = NULL

for (seed in 80:89){
  RMSD_threshold = 0.15
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    cor_all_N10000_r2_p10_J15_rmsd_0.15 = rbind(cor_all_N10000_r2_p10_J15_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_g_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    g_MSE_N10000_r2_p10_J15_rmsd_0.15 = rbind(g_MSE_N10000_r2_p10_J15_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    a_MSE_N10000_r2_p10_J15_rmsd_0.15 = rbind(a_MSE_N10000_r2_p10_J15_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    b_MSE_N10000_r2_p10_J15_rmsd_0.15 = rbind(b_MSE_N10000_r2_p10_J15_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    mu_MSE_N10000_r2_p10_J15_rmsd_0.15 = rbind(mu_MSE_N10000_r2_p10_J15_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    sigma_MSE_N10000_r2_p10_J15_rmsd_0.15 = rbind(sigma_MSE_N10000_r2_p10_J15_rmsd_0.15, temp)
    
  }
}



##### n = 20000, J = 15, p = 10 #######

results_path <- "simulation-output/proposed_S2/"
setwd(results_path)

#1 p10 rho0.3 N20000
tran_g_dat=NULL
a_dat=NULL
b_dat=NULL
mu_dat=NULL
sigma_dat=NULL

g_MSE_N20000_r2_p10_J15=NULL
a_MSE_N20000_r2_p10_J15=NULL
b_MSE_N20000_r2_p10_J15=NULL
mu_MSE_N20000_r2_p10_J15=NULL
sigma_MSE_N20000_r2_p10_J15=NULL
pvalgam_N20000_r2_p10_J15=NULL

cor_all_N20000_r2_p10_J15=NULL

N=20000
p=10
J=15
rho=2
for (seed in 80:89){
  if(file.exists(paste0('seed',seed,'_tran_gamma_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))){
    # temp_p10_rho0.3_N20000 = read.csv(paste0('seed',seed,'_tran_gamma_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    # tran_g_dat = rbind(tran_g_dat, temp_p10_rho0.3_N20000)
    # 
    # temp_p10_rho0.3_N20000 = read.csv(paste0('seed',seed,'_alpha_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    # a_dat = rbind(a_dat, temp_p10_rho0.3_N20000)
    # 
    # temp_p10_rho0.3_N20000 = read.csv(paste0('seed',seed,'_tran_d_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    # b_dat = rbind(b_dat, temp_p10_rho0.3_N20000)
    # 
    # temp_p10_rho0.3_N20000 = read.csv(paste0('seed',seed,'_tran_mu_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    # mu_dat = rbind(mu_dat, temp_p10_rho0.3_N20000)
    # 
    # temp_p10_rho0.3_N20000 = read.csv(paste0('seed',seed,'_sigma_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    # sigma_dat = rbind(sigma_dat, temp_p10_rho0.3_N20000)
    
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    cor_all_N20000_r2_p10_J15 = rbind(cor_all_N20000_r2_p10_J15, temp)
    
    temp = read.csv(paste0('seed',seed,'_g_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    g_MSE_N20000_r2_p10_J15 = rbind(g_MSE_N20000_r2_p10_J15, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    a_MSE_N20000_r2_p10_J15 = rbind(a_MSE_N20000_r2_p10_J15, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    b_MSE_N20000_r2_p10_J15 = rbind(b_MSE_N20000_r2_p10_J15, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    mu_MSE_N20000_r2_p10_J15 = rbind(mu_MSE_N20000_r2_p10_J15, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    sigma_MSE_N20000_r2_p10_J15 = rbind(sigma_MSE_N20000_r2_p10_J15, temp)
    
    # temp = read.csv(paste0('seed',seed,'_pvaluegam_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'), header = F)
    # pvalgam_N20000_r2_p10_J15 = rbind(pvalgam_N20000_r2_p10_J15, temp)
    
  }
}

#eval_typeI_gam(pvalgam_N20000_r2_p10_J15, g=g, sig)


results_path <- "simulation-output/baseline_S2"
setwd(results_path)

#1 p10 rho0.3 N20000

a_MSE_N20000_r2_p10_J15_bcmk=NULL
b_MSE_N20000_r2_p10_J15_bcmk=NULL
mu_MSE_N20000_r2_p10_J15_bcmk=NULL
sigma_MSE_N20000_r2_p10_J15_bcmk=NULL

cor_all_N20000_r2_p10_J15_bcmk=NULL

N=20000
p=10
J=15
rho=2
for (seed in 80:89){
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))){
    
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    cor_all_N20000_r2_p10_J15_bcmk = rbind(cor_all_N20000_r2_p10_J15_bcmk, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    a_MSE_N20000_r2_p10_J15_bcmk = rbind(a_MSE_N20000_r2_p10_J15_bcmk, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    b_MSE_N20000_r2_p10_J15_bcmk = rbind(b_MSE_N20000_r2_p10_J15_bcmk, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    mu_MSE_N20000_r2_p10_J15_bcmk = rbind(mu_MSE_N20000_r2_p10_J15_bcmk, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    sigma_MSE_N20000_r2_p10_J15_bcmk = rbind(sigma_MSE_N20000_r2_p10_J15_bcmk, temp)
    
    
  }
}

###########################################################################
##################### RMSD method #########################################
###########################################################################
results_path <- "simulation-output/rmsd_S2"
setwd(results_path)

N=20000
rho=2
p=10
J=15
cor_all_N20000_r2_p10_J15_rmsd_0.05 = NULL
g_MSE_N20000_r2_p10_J15_rmsd_0.05 = NULL
a_MSE_N20000_r2_p10_J15_rmsd_0.05 = NULL
b_MSE_N20000_r2_p10_J15_rmsd_0.05 = NULL
mu_MSE_N20000_r2_p10_J15_rmsd_0.05 = NULL
sigma_MSE_N20000_r2_p10_J15_rmsd_0.05 = NULL

for (seed in 80:89){
  RMSD_threshold = 0.05
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    cor_all_N20000_r2_p10_J15_rmsd_0.05 = rbind(cor_all_N20000_r2_p10_J15_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_g_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    g_MSE_N20000_r2_p10_J15_rmsd_0.05 = rbind(g_MSE_N20000_r2_p10_J15_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    a_MSE_N20000_r2_p10_J15_rmsd_0.05 = rbind(a_MSE_N20000_r2_p10_J15_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    b_MSE_N20000_r2_p10_J15_rmsd_0.05 = rbind(b_MSE_N20000_r2_p10_J15_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    mu_MSE_N20000_r2_p10_J15_rmsd_0.05 = rbind(mu_MSE_N20000_r2_p10_J15_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    sigma_MSE_N20000_r2_p10_J15_rmsd_0.05 = rbind(sigma_MSE_N20000_r2_p10_J15_rmsd_0.05, temp)
    
  }
}

N=20000
rho=2
p=10
J=15
cor_all_N20000_r2_p10_J15_rmsd_0.1 = NULL
g_MSE_N20000_r2_p10_J15_rmsd_0.1 = NULL
a_MSE_N20000_r2_p10_J15_rmsd_0.1 = NULL
b_MSE_N20000_r2_p10_J15_rmsd_0.1 = NULL
mu_MSE_N20000_r2_p10_J15_rmsd_0.1 = NULL
sigma_MSE_N20000_r2_p10_J15_rmsd_0.1 = NULL

for (seed in 80:89){
  RMSD_threshold = 0.1
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    cor_all_N20000_r2_p10_J15_rmsd_0.1 = rbind(cor_all_N20000_r2_p10_J15_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_g_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    g_MSE_N20000_r2_p10_J15_rmsd_0.1 = rbind(g_MSE_N20000_r2_p10_J15_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    a_MSE_N20000_r2_p10_J15_rmsd_0.1 = rbind(a_MSE_N20000_r2_p10_J15_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    b_MSE_N20000_r2_p10_J15_rmsd_0.1 = rbind(b_MSE_N20000_r2_p10_J15_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    mu_MSE_N20000_r2_p10_J15_rmsd_0.1 = rbind(mu_MSE_N20000_r2_p10_J15_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    sigma_MSE_N20000_r2_p10_J15_rmsd_0.1 = rbind(sigma_MSE_N20000_r2_p10_J15_rmsd_0.1, temp)
    
  }
}

N=20000
rho=2
p=10
J=15
cor_all_N20000_r2_p10_J15_rmsd_0.15 = NULL
g_MSE_N20000_r2_p10_J15_rmsd_0.15 = NULL
a_MSE_N20000_r2_p10_J15_rmsd_0.15 = NULL
b_MSE_N20000_r2_p10_J15_rmsd_0.15 = NULL
mu_MSE_N20000_r2_p10_J15_rmsd_0.15 = NULL
sigma_MSE_N20000_r2_p10_J15_rmsd_0.15 = NULL

for (seed in 80:89){
  RMSD_threshold = 0.15
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    cor_all_N20000_r2_p10_J15_rmsd_0.15 = rbind(cor_all_N20000_r2_p10_J15_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_g_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    g_MSE_N20000_r2_p10_J15_rmsd_0.15 = rbind(g_MSE_N20000_r2_p10_J15_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    a_MSE_N20000_r2_p10_J15_rmsd_0.15 = rbind(a_MSE_N20000_r2_p10_J15_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    b_MSE_N20000_r2_p10_J15_rmsd_0.15 = rbind(b_MSE_N20000_r2_p10_J15_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    mu_MSE_N20000_r2_p10_J15_rmsd_0.15 = rbind(mu_MSE_N20000_r2_p10_J15_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    sigma_MSE_N20000_r2_p10_J15_rmsd_0.15 = rbind(sigma_MSE_N20000_r2_p10_J15_rmsd_0.15, temp)
    
  }
}


##### n = 10000, J = 30, p = 20 #######

results_path <- "simulation-output/proposed_S4"
setwd(results_path)

#1 p10 rho0.3 N10000
tran_g_dat=NULL
a_dat=NULL
b_dat=NULL
mu_dat=NULL
sigma_dat=NULL

g_MSE_N10000_r2_p20_J30=NULL
a_MSE_N10000_r2_p20_J30=NULL
b_MSE_N10000_r2_p20_J30=NULL
mu_MSE_N10000_r2_p20_J30=NULL
sigma_MSE_N10000_r2_p20_J30=NULL
pvalgam_N10000_r2_p20_J30=NULL

cor_all_N10000_r2_p20_J30=NULL

N=10000
p=20
J=30
rho=2
for (seed in 80:89){
  if(file.exists(paste0('seed',seed,'_tran_gamma_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    cor_all_N10000_r2_p20_J30 = rbind(cor_all_N10000_r2_p20_J30, temp)
    
    temp = read.csv(paste0('seed',seed,'_g_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    g_MSE_N10000_r2_p20_J30 = rbind(g_MSE_N10000_r2_p20_J30, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    a_MSE_N10000_r2_p20_J30 = rbind(a_MSE_N10000_r2_p20_J30, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    b_MSE_N10000_r2_p20_J30 = rbind(b_MSE_N10000_r2_p20_J30, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    mu_MSE_N10000_r2_p20_J30 = rbind(mu_MSE_N10000_r2_p20_J30, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    sigma_MSE_N10000_r2_p20_J30 = rbind(sigma_MSE_N10000_r2_p20_J30, temp)
    
    # temp = read.csv(paste0('seed',seed,'_pvaluegam_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'), header = F)
    # pvalgam_N10000_r2_p20_J30 = rbind(pvalgam_N10000_r2_p20_J30, temp)
    
  }
}

#eval_typeI_gam(pvalgam_N10000_r2_p20_J30, g=g, sig)


results_path <- "simulation-output/baseline_S4"
setwd(results_path)

#1 p10 rho0.3 N10000

a_MSE_N10000_r2_p20_J30_bcmk=NULL
b_MSE_N10000_r2_p20_J30_bcmk=NULL
mu_MSE_N10000_r2_p20_J30_bcmk=NULL
sigma_MSE_N10000_r2_p20_J30_bcmk=NULL

cor_all_N10000_r2_p20_J30_bcmk=NULL

N=10000
p=20
J=30
rho=2
for (seed in 80:89){
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))){
    
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    cor_all_N10000_r2_p20_J30_bcmk = rbind(cor_all_N10000_r2_p20_J30_bcmk, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    a_MSE_N10000_r2_p20_J30_bcmk = rbind(a_MSE_N10000_r2_p20_J30_bcmk, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    b_MSE_N10000_r2_p20_J30_bcmk = rbind(b_MSE_N10000_r2_p20_J30_bcmk, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    mu_MSE_N10000_r2_p20_J30_bcmk = rbind(mu_MSE_N10000_r2_p20_J30_bcmk, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    sigma_MSE_N10000_r2_p20_J30_bcmk = rbind(sigma_MSE_N10000_r2_p20_J30_bcmk, temp)
    
    
  }
}

###########################################################################
##################### RMSD method #########################################
###########################################################################
results_path <- "simulation-output/rmsd_S4"
setwd(results_path)

N=10000
rho=2
p=20
J=30
cor_all_N10000_r2_p20_J30_rmsd_0.05 = NULL
g_MSE_N10000_r2_p20_J30_rmsd_0.05 = NULL
a_MSE_N10000_r2_p20_J30_rmsd_0.05 = NULL
b_MSE_N10000_r2_p20_J30_rmsd_0.05 = NULL
mu_MSE_N10000_r2_p20_J30_rmsd_0.05 = NULL
sigma_MSE_N10000_r2_p20_J30_rmsd_0.05 = NULL

for (seed in 80:89){
  RMSD_threshold = 0.05
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    cor_all_N10000_r2_p20_J30_rmsd_0.05 = rbind(cor_all_N10000_r2_p20_J30_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_g_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    g_MSE_N10000_r2_p20_J30_rmsd_0.05 = rbind(g_MSE_N10000_r2_p20_J30_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    a_MSE_N10000_r2_p20_J30_rmsd_0.05 = rbind(a_MSE_N10000_r2_p20_J30_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    b_MSE_N10000_r2_p20_J30_rmsd_0.05 = rbind(b_MSE_N10000_r2_p20_J30_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    mu_MSE_N10000_r2_p20_J30_rmsd_0.05 = rbind(mu_MSE_N10000_r2_p20_J30_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    sigma_MSE_N10000_r2_p20_J30_rmsd_0.05 = rbind(sigma_MSE_N10000_r2_p20_J30_rmsd_0.05, temp)
    
  }
}

N=10000
rho=2
p=20
J=30
cor_all_N10000_r2_p20_J30_rmsd_0.1 = NULL
g_MSE_N10000_r2_p20_J30_rmsd_0.1 = NULL
a_MSE_N10000_r2_p20_J30_rmsd_0.1 = NULL
b_MSE_N10000_r2_p20_J30_rmsd_0.1 = NULL
mu_MSE_N10000_r2_p20_J30_rmsd_0.1 = NULL
sigma_MSE_N10000_r2_p20_J30_rmsd_0.1 = NULL

for (seed in 80:89){
  RMSD_threshold = 0.1
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    cor_all_N10000_r2_p20_J30_rmsd_0.1 = rbind(cor_all_N10000_r2_p20_J30_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_g_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    g_MSE_N10000_r2_p20_J30_rmsd_0.1 = rbind(g_MSE_N10000_r2_p20_J30_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    a_MSE_N10000_r2_p20_J30_rmsd_0.1 = rbind(a_MSE_N10000_r2_p20_J30_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    b_MSE_N10000_r2_p20_J30_rmsd_0.1 = rbind(b_MSE_N10000_r2_p20_J30_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    mu_MSE_N10000_r2_p20_J30_rmsd_0.1 = rbind(mu_MSE_N10000_r2_p20_J30_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    sigma_MSE_N10000_r2_p20_J30_rmsd_0.1 = rbind(sigma_MSE_N10000_r2_p20_J30_rmsd_0.1, temp)
    
  }
}

N=10000
rho=2
p=20
J=30
cor_all_N10000_r2_p20_J30_rmsd_0.15 = NULL
g_MSE_N10000_r2_p20_J30_rmsd_0.15 = NULL
a_MSE_N10000_r2_p20_J30_rmsd_0.15 = NULL
b_MSE_N10000_r2_p20_J30_rmsd_0.15 = NULL
mu_MSE_N10000_r2_p20_J30_rmsd_0.15 = NULL
sigma_MSE_N10000_r2_p20_J30_rmsd_0.15 = NULL

for (seed in 80:89){
  RMSD_threshold = 0.15
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    cor_all_N10000_r2_p20_J30_rmsd_0.15 = rbind(cor_all_N10000_r2_p20_J30_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_g_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    g_MSE_N10000_r2_p20_J30_rmsd_0.15 = rbind(g_MSE_N10000_r2_p20_J30_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    a_MSE_N10000_r2_p20_J30_rmsd_0.15 = rbind(a_MSE_N10000_r2_p20_J30_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    b_MSE_N10000_r2_p20_J30_rmsd_0.15 = rbind(b_MSE_N10000_r2_p20_J30_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    mu_MSE_N10000_r2_p20_J30_rmsd_0.15 = rbind(mu_MSE_N10000_r2_p20_J30_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    sigma_MSE_N10000_r2_p20_J30_rmsd_0.15 = rbind(sigma_MSE_N10000_r2_p20_J30_rmsd_0.15, temp)
    
  }
}

##### n = 20000, J = 30, p = 20 #######


results_path <- "simulation-output/proposed_S4"
setwd(results_path)

#1 p10 rho0.3 N20000
tran_g_dat=NULL
a_dat=NULL
b_dat=NULL
mu_dat=NULL
sigma_dat=NULL

g_MSE_N20000_r2_p20_J30=NULL
a_MSE_N20000_r2_p20_J30=NULL
b_MSE_N20000_r2_p20_J30=NULL
mu_MSE_N20000_r2_p20_J30=NULL
sigma_MSE_N20000_r2_p20_J30=NULL
pvalgam_N20000_r2_p20_J30=NULL

cor_all_N20000_r2_p20_J30=NULL

N=20000
p=20
J=30
rho=2
for (seed in 80:89){
  if(file.exists(paste0('seed',seed,'_tran_gamma_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))){
    # temp_p10_rho0.3_N20000 = read.csv(paste0('seed',seed,'_tran_gamma_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    # tran_g_dat = rbind(tran_g_dat, temp_p10_rho0.3_N20000)
    # 
    # temp_p10_rho0.3_N20000 = read.csv(paste0('seed',seed,'_alpha_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    # a_dat = rbind(a_dat, temp_p10_rho0.3_N20000)
    # 
    # temp_p10_rho0.3_N20000 = read.csv(paste0('seed',seed,'_tran_d_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    # b_dat = rbind(b_dat, temp_p10_rho0.3_N20000)
    # 
    # temp_p10_rho0.3_N20000 = read.csv(paste0('seed',seed,'_tran_mu_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    # mu_dat = rbind(mu_dat, temp_p10_rho0.3_N20000)
    # 
    # temp_p10_rho0.3_N20000 = read.csv(paste0('seed',seed,'_sigma_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    # sigma_dat = rbind(sigma_dat, temp_p10_rho0.3_N20000)
    
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    cor_all_N20000_r2_p20_J30 = rbind(cor_all_N20000_r2_p20_J30, temp)
    
    temp = read.csv(paste0('seed',seed,'_g_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    g_MSE_N20000_r2_p20_J30 = rbind(g_MSE_N20000_r2_p20_J30, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    a_MSE_N20000_r2_p20_J30 = rbind(a_MSE_N20000_r2_p20_J30, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    b_MSE_N20000_r2_p20_J30 = rbind(b_MSE_N20000_r2_p20_J30, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    mu_MSE_N20000_r2_p20_J30 = rbind(mu_MSE_N20000_r2_p20_J30, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    sigma_MSE_N20000_r2_p20_J30 = rbind(sigma_MSE_N20000_r2_p20_J30, temp)
    
    # temp = read.csv(paste0('seed',seed,'_pvaluegam_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'), header = F)
    # pvalgam_N20000_r2_p20_J30 = rbind(pvalgam_N20000_r2_p20_J30, temp)
    
  }
}

#eval_typeI_gam(pvalgam_N20000_r2_p20_J30, g=g, sig)


results_path <- "simulation-output/baseline_S4"
setwd(results_path)

#1 p10 rho0.3 N20000

a_MSE_N20000_r2_p20_J30_bcmk=NULL
b_MSE_N20000_r2_p20_J30_bcmk=NULL
mu_MSE_N20000_r2_p20_J30_bcmk=NULL
sigma_MSE_N20000_r2_p20_J30_bcmk=NULL

cor_all_N20000_r2_p20_J30_bcmk=NULL

N=20000
p=20
J=30
rho=2
for (seed in 80:89){
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))){
    
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    cor_all_N20000_r2_p20_J30_bcmk = rbind(cor_all_N20000_r2_p20_J30_bcmk, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    a_MSE_N20000_r2_p20_J30_bcmk = rbind(a_MSE_N20000_r2_p20_J30_bcmk, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    b_MSE_N20000_r2_p20_J30_bcmk = rbind(b_MSE_N20000_r2_p20_J30_bcmk, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    mu_MSE_N20000_r2_p20_J30_bcmk = rbind(mu_MSE_N20000_r2_p20_J30_bcmk, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'),header = F)
    sigma_MSE_N20000_r2_p20_J30_bcmk = rbind(sigma_MSE_N20000_r2_p20_J30_bcmk, temp)
    
    
  }
}

###########################################################################
##################### RMSD method #########################################
###########################################################################
results_path <- "simulation-output/rmsd_S4"
setwd(results_path)

N=20000
rho=2
p=20
J=30
cor_all_N20000_r2_p20_J30_rmsd_0.05 = NULL
g_MSE_N20000_r2_p20_J30_rmsd_0.05 = NULL
a_MSE_N20000_r2_p20_J30_rmsd_0.05 = NULL
b_MSE_N20000_r2_p20_J30_rmsd_0.05 = NULL
mu_MSE_N20000_r2_p20_J30_rmsd_0.05 = NULL
sigma_MSE_N20000_r2_p20_J30_rmsd_0.05 = NULL

for (seed in 80:89){
  RMSD_threshold = 0.05
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    cor_all_N20000_r2_p20_J30_rmsd_0.05 = rbind(cor_all_N20000_r2_p20_J30_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_g_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    g_MSE_N20000_r2_p20_J30_rmsd_0.05 = rbind(g_MSE_N20000_r2_p20_J30_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    a_MSE_N20000_r2_p20_J30_rmsd_0.05 = rbind(a_MSE_N20000_r2_p20_J30_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    b_MSE_N20000_r2_p20_J30_rmsd_0.05 = rbind(b_MSE_N20000_r2_p20_J30_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    mu_MSE_N20000_r2_p20_J30_rmsd_0.05 = rbind(mu_MSE_N20000_r2_p20_J30_rmsd_0.05, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    sigma_MSE_N20000_r2_p20_J30_rmsd_0.05 = rbind(sigma_MSE_N20000_r2_p20_J30_rmsd_0.05, temp)
    
  }
}

N=20000
rho=2
p=20
J=30
cor_all_N20000_r2_p20_J30_rmsd_0.1 = NULL
g_MSE_N20000_r2_p20_J30_rmsd_0.1 = NULL
a_MSE_N20000_r2_p20_J30_rmsd_0.1 = NULL
b_MSE_N20000_r2_p20_J30_rmsd_0.1 = NULL
mu_MSE_N20000_r2_p20_J30_rmsd_0.1 = NULL
sigma_MSE_N20000_r2_p20_J30_rmsd_0.1 = NULL

for (seed in 80:89){
  RMSD_threshold = 0.1
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    cor_all_N20000_r2_p20_J30_rmsd_0.1 = rbind(cor_all_N20000_r2_p20_J30_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_g_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    g_MSE_N20000_r2_p20_J30_rmsd_0.1 = rbind(g_MSE_N20000_r2_p20_J30_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    a_MSE_N20000_r2_p20_J30_rmsd_0.1 = rbind(a_MSE_N20000_r2_p20_J30_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    b_MSE_N20000_r2_p20_J30_rmsd_0.1 = rbind(b_MSE_N20000_r2_p20_J30_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    mu_MSE_N20000_r2_p20_J30_rmsd_0.1 = rbind(mu_MSE_N20000_r2_p20_J30_rmsd_0.1, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    sigma_MSE_N20000_r2_p20_J30_rmsd_0.1 = rbind(sigma_MSE_N20000_r2_p20_J30_rmsd_0.1, temp)
    
  }
}

N=20000
rho=2
p=20
J=30
cor_all_N20000_r2_p20_J30_rmsd_0.15 = NULL
g_MSE_N20000_r2_p20_J30_rmsd_0.15 = NULL
a_MSE_N20000_r2_p20_J30_rmsd_0.15 = NULL
b_MSE_N20000_r2_p20_J30_rmsd_0.15 = NULL
mu_MSE_N20000_r2_p20_J30_rmsd_0.15 = NULL
sigma_MSE_N20000_r2_p20_J30_rmsd_0.15 = NULL

for (seed in 80:89){
  RMSD_threshold = 0.15
  if(file.exists(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))){
    temp = read.csv(paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    cor_all_N20000_r2_p20_J30_rmsd_0.15 = rbind(cor_all_N20000_r2_p20_J30_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_g_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    g_MSE_N20000_r2_p20_J30_rmsd_0.15 = rbind(g_MSE_N20000_r2_p20_J30_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    a_MSE_N20000_r2_p20_J30_rmsd_0.15 = rbind(a_MSE_N20000_r2_p20_J30_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    b_MSE_N20000_r2_p20_J30_rmsd_0.15 = rbind(b_MSE_N20000_r2_p20_J30_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    mu_MSE_N20000_r2_p20_J30_rmsd_0.15 = rbind(mu_MSE_N20000_r2_p20_J30_rmsd_0.15, temp)
    
    temp = read.csv(paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv') ,header = F)
    sigma_MSE_N20000_r2_p20_J30_rmsd_0.15 = rbind(sigma_MSE_N20000_r2_p20_J30_rmsd_0.15, temp)
    
  }
}
