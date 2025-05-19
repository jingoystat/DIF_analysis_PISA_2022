source("functions_rmsd.R")

#args <- commandArgs(TRUE)
seed=80 # 80 82 84 86 88
rho=2 # 2
N=10000 # 10000 20000
J=15 # 15
p=10 # 10
RMSD_threshold=0.05 # 0.05 0.1 0.15

a=rep(c(1, 1.2, 1.4,1.6,1.8), J/5)
b=rep(c(0.8, 0.2, -0.4, -1, 1), J/5)
g=matrix(0, J, p)

set.seed(3)
row_indices_1 <- 1:(J/5)
row_indices_2 <- (J/5+2):(floor(J/2)+J/5)  
row_indices_3 <- c(J/5+1, (floor(J/2)+J/5+1):(floor(J/2)+J/5+1))  
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

x=matrix(0, nrow = N, ncol=p)
size = N/p
for(j in 1:p){
  x[((j-1)*size+1):(j*size),j] = 1
}


mu.vec=rep(c( 0.06, -0.06, -0.08, -0.1, -0.12, -0.2, -0.22, -0.24, -0.26, -0.28), p/10)
mu.vec = (mu.vec - mean(mu.vec))
sigma.vec=rep(c(0.8, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 1.1), p/10)

r=coverage_proposed_inference(N, sig=0.05, beta.vec=b, alpha.vec=a, gamma=g, mu=mu.vec, sigma=sigma.vec, x=x, ite=1, RMSD_threshold, sig_list=seq(0, 1, 0.02), seed)


cor_all = c()
for(i in 1:dim(r$mu)[1]){
  cor_all = c(cor_all, cor(as.numeric(r$mu[i,]), mu.vec, method = "kendall")
  )
}

write.table(g, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_true_g_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))
write.table(cor_all, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))

write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_alpha_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_beta_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_gamma_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_mu_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))
write.table(r$sigma, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_sigma_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))
write.table(r$discovery, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_discovery_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))

write.table(r$g_MSE, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_g_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))
write.table(r$a_MSE, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))
write.table(r$b_MSE, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))
write.table(r$mu_MSE, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))
write.table(r$sigma_MSE, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'_threshold',RMSD_threshold,'.csv'))

