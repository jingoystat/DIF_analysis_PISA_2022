source("functions.R")

#args <- commandArgs(TRUE)
seed=80 # 80 81 82 83 84 85 86 87 88 89
rho=2 # 2
N=10000 # 10000 20000
J=30 # 30
p=20 # 20

a=rep(c(1, 1.2, 1.4,1.6,1.8), J/5)
b=rep(c(0.8, 0.2, -0.4, -1, 1), J/5)
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


x=matrix(0, nrow = N, ncol=p)
size = N/p
for(j in 1:p){
  x[((j-1)*size+1):(j*size),j] = 1
}


mu.vec=c( 0.07, -0.06, -0.08, -0.09, -0.1, -0.12, -0.14, -0.15, -0.16, -0.18, 
          -0.2, -0.21, -0.22, -0.24, -0.26, -0.28, -0.29, -0.3, -0.31, -0.32)
mu.vec = (mu.vec - mean(mu.vec))
sigma.vec=rep(c(0.8, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 1.1), p/10)

r=coverage_proposed_inference(seed=seed, N, sig=0.05, beta.vec=b, alpha.vec=a, gamma=g, mu=mu.vec, sigma=sigma.vec, x=x, ite=1, sig_list=seq(0, 1, 0.02))


cor_all = c()
for(i in 1:dim(r$tran_mu)[1]){
  cor_all = c(cor_all, cor(as.numeric(r$tran_mu[i,]), mu.vec, method = "kendall")
  )
}

write.table(cor_all, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_cor_all_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))
write.table(r$alpha, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_alpha_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))
write.table(r$beta, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_beta_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))
write.table(r$gamma, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_gamma_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))
write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_mu_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))
write.table(r$sigma, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_sigma_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))
write.table(r$tran_g, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_tran_gamma_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))
write.table(r$tran_mu, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_tran_mu_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))
write.table(r$tran_d, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_tran_d_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))

write.table(r$g_MSE, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_g_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))
write.table(r$a_MSE, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_a_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))
write.table(r$b_MSE, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_b_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))
write.table(r$mu_MSE, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_mu_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))
write.table(r$sigma_MSE, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('seed',seed,'_sigma_MSE_N', N,'_rho',rho,'_J',J,'_p',p,'.csv'))


