library(vegan)
library(gplots)
library(MASS)
library(fastGHQuad)
library(quantreg)
library(tidyverse)
library(mirt)

#Use computer scored items only
item_math <- c(
  'CM033Q01S',
  'CM474Q01S',
  'CM155Q01S',
  'CM155Q04S',
  'CM411Q01S',
  'CM411Q02S',
  'CM803Q01S',
  'CM442Q02S',
  'CM034Q01S',
  'CM305Q01S',
  'CM496Q01S',
  'CM496Q02S',
  'CM423Q01S',
  'CM192Q01S',
  'CM603Q01S',
  'CM571Q01S',
  'CM564Q01S',
  'CM564Q02S',
  'CM447Q01S',
  'CM273Q01S',
  'CM408Q01S',
  'CM420Q01S',
  'CM446Q01S',
  'DM446Q02C',
  'CM559Q01S',
  'DM828Q02C',
  'CM828Q03S',
  'CM464Q01S',
  'CM800Q01S',
  'CM982Q01S',
  'CM982Q02S',
  'CM982Q03S',
  'CM982Q04S',
  'CM992Q01S',
  'CM992Q02S',
  'DM992Q03C',
  'CM915Q01S',
  'CM915Q02S',
  'CM906Q01S',
  'DM00KQ02C',
  'CM909Q01S',
  'CM909Q02S',
  'CM909Q03S',
  'CM949Q01S',
  'CM949Q02S',
  'CM00GQ01S',
  'DM955Q01C',
  'DM955Q02C',
  'DM998Q02C',
  'CM905Q01S',
  'DM905Q02C',
  'CM919Q01S',
  'CM919Q02S',
  'CM954Q01S',
  'DM954Q02C',
  'CM954Q04S',
  'CM943Q01S',
  'CM943Q02S',
  'DM953Q02C',
  'CM953Q03S',
  'CM948Q01S',
  'CM948Q02S',
  'CM948Q03S',
  'CM936Q01S',
  'DM936Q02C',
  'CM939Q02S',
  'CM967Q01S',
  'CM967Q03S',
  'CMA133Q01S',
  'CMA102Q01S',
  'DMA102Q02C',
  'CMA102Q03S',
  'DMA103Q01C',
  'CMA103Q02S',
  'CMA103Q04S',
  'CMA105Q01S',
  'CMA105Q02S',
  'CMA105Q03S',
  'CMA105Q05S',
  'CMA107Q01S',
  'CMA107Q02S',
  'CMA149Q02S',
  'DMA149Q03C',
  'DMA153Q01C',
  'CMA153Q02S',
  'CMA114Q01S',
  'CMA121Q03S',
  'CMA136Q03S',
  'CMA119Q01S',
  'CMA119Q02S',
  'CMA120Q01S',
  'CMA120Q03S',  
  'CMA127Q03S',
  'CMA123Q02S',
  'CMA139Q01S',
  'CMA139Q02S',
  'CMA157Q02S',
  'CMA143Q01S',
  'CMA143Q02S',
  'CMA143Q03S',
  'CMA143Q04S',
  'CMA144Q01S',
  'CMA101Q02S',
  'CMA135Q01S',
  'CMA135Q02S',
  'CMA135Q04S',
  'CMA137Q01S',
  'CMA137Q03S',
  'CMA161Q02S',
  'CMA108Q01S',
  'CMA108Q02S',
  'CMA108Q03S',
  'CMA146Q01S',
  'CMA146Q02S',
  'CMA110Q01S',
  'CMA110Q02S',
  'CMA116Q01S',
  'CMA117Q01S',
  'CMA117Q02S',
  'CMA151Q02S',
  'CMA125Q03S',
  'CMA145Q01S',
  'CMA145Q02S',
  'CMA129Q01S',
  'CMA129Q03S',
  'CMA131Q01S',
  'CMA131Q02S',
  'CMA131Q04S',
  'CMA150Q01S',
  'CMA150Q02S',
  'CMA160Q03S',
  'CMA140Q01S',
  'CMA140Q02S',
  'CMA140Q03S',
  'CMA132Q02S',
  'CMA147Q01S',
  'CMA147Q02S',
  'CMA147Q03S',
  'CMA147Q04S',
  'CMA134Q03S',
  'CMA148Q01S',
  'CMA148Q02S',
  'CMA154Q01S',
  'CMA154Q02S',
  'CMA158Q01S',
  'CMA109Q01S',
  'CMA109Q02S',
  'CMA109Q03S',
  'CMA111Q01S',
  'CMA111Q02S',
  'DMA111Q03C',
  'CMA113Q01S',
  'CMA113Q03S',
  'CMA115Q01S',
  'CMA115Q02S',
  'CMA112Q01S',
  'CMA126Q03S',
  'CMA126Q04S',
  'CMA138Q01S',
  'CMA138Q02S',
  'CMA130Q02S',
  'CMA130Q03S',
  'CMA162Q01S',
  'CMA141Q01S',
  'CMA141Q02S',
  'CMA141Q03S',
  'CMA141Q04S',
  'CMA142Q01S',
  'CMA124Q01S')



solve_optimization_quantile_c1 <- function(gamma, a, b) {
  J <- nrow(gamma) # Number of rows (j indices)
  p <- ncol(gamma) # Number of columns (k indices)
  
  # Create the response vector
  y <- as.vector(gamma) + rep(b, p)
  
  # Create the design matrix for the predictors
  X1 <- matrix(0, nrow = J*p, ncol = p)
  X2 <- matrix(0, nrow = J*p, ncol = J)
  
  # Fill in the design matrix
  for(k in 1:p){
    for(j in 1:J){
      X1[(k-1)*J+j, k] = a[j]
      X2[(k-1)*J+j, j] = 1
    }
  }
  X1 = X1[, -1]
  #print(X1)
  #print(X2)
  # Run quantile regression at the median (tau = 0.5)
  fit <- rq(y ~ -1 + X1 + X2, tau=0.5) # Use -1 to remove intercept since we handle it manually
  
  # Extract the coefficients for d_j and c_k
  coefficients <- coef(fit)
  c_optimal <- coefficients[1:(p-1)]
  d_optimal <- coefficients[p:(p+J-1)]
  
  return(list(c = c_optimal, d = d_optimal))
}


length(item_math)#169

#Read in data
dat=read.csv(file = 'PISA_2022_math_sas.csv',  sep = ",", header = T)
dat=dat[,-1]
colnames(dat)=NULL
dat=data.matrix(dat)

x=read.csv(file = 'countries_2022_math_sas.csv', sep = ",", header = T)
x=x[,-1]
colnames(x)=NULL
x=data.matrix(x)
covariates=x
weights=read.csv('PISA_2022_math_samplingweights_sas.csv', header = T)[,1]


N=dim(dat)[1]
J=dim(dat)[2]
p=dim(x)[2]

####################################
group_size=N/p
x_group=c()
for(i in 1:N){
  x_group = c(x_group, which(x[i,]!=0))
}

dat1=as.data.frame(cbind(x_group, dat))
mod <- TAM::tam.mml.2pl(dat1[,-1], group=dat1[,1] )
rmod <- CDM::IRT.RMSD(mod)

#RMSD_mtx = rmsd(a=alpha.vec,  b=beta.vec, g=gamma, x=x, dat=dat, z, mu=mu, sigma=sigma, w)
RMSD_mtx = rmod$RMSD[,2:(p+1)]

source("functions_math_rmsd.R")

args <- commandArgs(TRUE)
RMSD_threshold=as.numeric(args[[1]])

##################################

a_ini=read.csv('a_ini_math_2022.csv', header=F)[,1]
b_ini=rnorm(J) 
g_ini=matrix(runif(J*p), J, p)
mu_ini=runif(p, -1, 1)
sigma_ini=runif(p, 0.2, 1)

a=a_ini
b=b_ini
g=as.matrix(g_ini)
mu=mu_ini
sigma=sigma_ini


gauss=gaussHermiteData(n=35)#We use 35 Gauss-Hermite quadrature for integral approximation in E-step
w=gauss$w
z=gauss$x

r=EM_2PL_inference_lrt(a=a_ini,  b=b_ini, g=g, x=x, dat=dat, z, mu=mu, sigma=sigma, w, ite=ite, RMSD_mtx, RMSD_threshold=RMSD_threshold, tol = 0.1)
r$mu
r$sigma
r$alpha.vec

write.table(r$mu, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('mu_math_RMSD_2022_',RMSD_threshold,'.csv'))
write.table(r$sigma, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('sigma_math_RMSD_2022_',RMSD_threshold,'.csv'))
write.table(r$alpha.vec, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('a_math_RMSD_2022_',RMSD_threshold,'.csv'))
write.table(r$beta.vec, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('b_math_RMSD_2022_',RMSD_threshold,'.csv'))
write.table(r$gamma.vec, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('gamma_math_RMSD_2022_',RMSD_threshold,'.csv'))
write.table(r$post, sep=",",  col.names=FALSE, row.names=FALSE, file = paste0('post_math_RMSD_2022_',RMSD_threshold,'.csv'))


