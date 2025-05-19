library(fastGHQuad)
library(quantreg)
library(MASS)
library(foreach)
library(doParallel)

gauss=gaussHermiteData(n=35)#We use 35 Gaussian quadrature
w=gauss$w
z=gauss$x

if_zero=function(x.vec){
  if (all(x.vec==0)){
    return(1)
  }
  else{return(0)}
}

### Estimation ###
#Transform Gauss Hermite Quadrature z into th, th_k=sqrt(2)*sigma*z_k+mu for k=1,...,n for all i=1,...,N.
#Now x is a N by p binary matrix
#mu and sigma are a p-dim vectors
#Return matrix Theta of dim N by n.
tran=function(z, x, mu, sigma){
  n=length(z)
  N=dim(x)[1]
  p=dim(x)[2]
  mean.vec=apply((rep(1, N)%*%t(mu))*x, 1, sum) #N dim vector
  var.vec=apply((rep(1, N)%*%t(sigma^2))*x, 1, sum) #N dim vector
  #zero.vec=apply(x, 1, if_zero)
  #var.vec=var.vec+zero.vec
  sd.vec=var.vec^0.5
  Theta=sqrt(2)*(sd.vec%*%t(rep(1, n)))*(rep(1, N)%*%t(z))+mean.vec%*%t(rep(1, n)) #N by n 
  return(Theta)
}

evaluate_prod=function(v){
  L=length(v)
  return(exp(sum(v, na.rm = TRUE)))
}



#Evaluate the posterior quadrature values, return a matrix N by n.
#a, b are alpha and beta of length J. dat is of dim N by J. g is gamma of dim J by p.
#Now x is a N by p binary matrix
#mu and sigma are a p-dim vectors
#z, w are the Gauss Hermite quadratures and weights of length n, 
quad=function(a, b, g, x, dat, z, mu, sigma, w){ #z is quadrature of dim n
  n=length(z)
  N=dim(dat)[1]
  J=dim(dat)[2]
  p=dim(x)[2]
  Theta=tran(z, x, mu, sigma) #Transform Gauss Hermite Quadrature z into th, th_i=sqrt(2)*z_i+mu for i=1,...,n
  th=Theta[, 1]
  ta=rep(1, N)%*%t(a)
  tb=rep(1, N)%*%t(b)
  tg=x%*%t(g)
  
  temp=(th%*%t(rep(1, J)))*ta+ tb + tg #N by J
  prob=exp(temp)/(1+exp(temp)) #dim N by J
  q=apply(log((prob^dat)*((1-prob)^(1-dat))), 1, evaluate_prod) #length N vector
  for (k in 2:n){ #matrix evaluate
    th=Theta[, k]
    temp=(th%*%t(rep(1, J)))*ta+tb+tg
    prob=exp(temp)/(1+exp(temp))
    q_c=apply(log((prob^dat)*((1-prob)^(1-dat))), 1, evaluate_prod)
    q=cbind(q, q_c)
  }
  Z=rowSums(q*(rep(1, N)%*%t(w)))
  return((q*(rep(1, N)%*%t(w)))/(Z%*%t(rep(1, n))))
}






#Evaluate target function at j=1 (i.e. anchor item), phi_j=(a_j, b_j, g_j1,...,g_jp), i.e. negative log-likelihood at phi_j
#Y_j=Y[,j]
#weights are sampling weights
target_function1=function(phi_j, x, Y_j, weights, z,  mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  #a_j=phi_j[1]
  b_j=phi_j
  g_j=rep(0, p)
  a_j=1
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu_pres, sigma_pres)
  #mean.v=apply(x*(rep(1, N)%*%t(mu)), 1, sum)
  #var.v=apply(x*(rep(1, N)%*%t(sigma^2)), 1, sum)
  #sd.v=var.v^0.5
  #t1=-0.5*((Theta-mean.v%*%t(rep(1,n)))/(sd.v%*%t(rep(1,n))))^2 #normal term in log-likelihood  #remove prior term
  temp=Theta*a_j+b_j*(rep(1,N)%*%t(rep(1, n)))+(x%*%g_j)%*%t(rep(1, n)) #N by n
  t2=Y_j%*%t(rep(1, n))*temp-log(1+exp(temp)) #The remaining term in log-likelihood
  #log_likelihood=t1+t2 #N by n
  log_likelihood=(weights%*%t(rep(1, n)))*t2 #N by n
  res=-sum(q*log_likelihood, na.rm = TRUE)
  return(res/N)
}

#Evaluate target function at phi_j=(a_j, b_j) given (g_j1, ... g_jp) at previous step values.
#Y_j=Y[,j]
#weights are sampling weights
target_function_ab=function(ab, g, x, Y_j, weights, z,  mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  a_j=ab
  b_j=0
  g_j=g
  #  g_j[1]=0
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu_pres, sigma_pres)
  #mean.v=apply(x*(rep(1, N)%*%t(mu)), 1, sum)
  #var.v=apply(x*(rep(1, N)%*%t(sigma^2)), 1, sum)
  #sd.v=var.v^0.5
  #t1=-0.5*((Theta-mean.v%*%t(rep(1,n)))/(sd.v%*%t(rep(1,n))))^2 #normal term in log-likelihood  #remove prior term
  temp=Theta*a_j+b_j*(rep(1,N)%*%t(rep(1, n)))+(x%*%g_j)%*%t(rep(1, n)) #N by n
  t2=Y_j%*%t(rep(1, n))*temp-log(1+exp(temp)) #The remaining term in log-likelihood
  #log_likelihood=t1+t2 #N by n
  log_likelihood=(weights%*%t(rep(1, n)))*t2 #N by n
  res=-sum(q*log_likelihood, na.rm = TRUE)
  return(res/N)
}

#Evaluate target function at g_jk, given (a_j, b_j, g_j1,...g_j(k-1),g_j(k+1),...g_jp) at previous step values,
#Y_k=Y[[i in group k],j]
#ab=(a_j, b_j)
#gp=(g_j1,...g_j(k-1))
#ga=(g_j(k+1),...g_jp)
#gk=g[k]
#weights are sampling weights
#ind_k is the indices of individuals belonging to group k
target_function_g1=function(ab, k, gk, gp, x_k, Y_k, weights, z, ind_k, mu_pres, sigma_pres, q){
  N_k=dim(x_k)[1]
  p=dim(x_k)[2]
  n=length(z)
  a_j=1
  b_j=0
  gp[k]=gk
  g_j=gp
  #  g_j[1] = -sum(g_j[-1])
  
  q_k=q[ind_k, ]
  weight_k=weights[ind_k]
  
  #Evaluate log-likelihood
  Theta=tran(z, x_k, mu_pres, sigma_pres)
  #mean.v=apply(x*(rep(1, N)%*%t(mu)), 1, sum)
  #var.v=apply(x*(rep(1, N)%*%t(sigma^2)), 1, sum)
  #sd.v=var.v^0.5
  #t1=-0.5*((Theta-mean.v%*%t(rep(1,n)))/(sd.v%*%t(rep(1,n))))^2 #normal term in log-likelihood  #remove prior term
  temp=Theta*a_j+b_j*(rep(1,N_k)%*%t(rep(1, n)))+(x_k%*%g_j)%*%t(rep(1, n)) #N_k by n
  t2=Y_k%*%t(rep(1, n))*temp-log(1+exp(temp)) #The remaining term in log-likelihood
  log_likelihood=(weight_k%*%t(rep(1, n)))*t2 #N_k by n
  res=-sum(q_k*log_likelihood, na.rm = TRUE)
  return(res/N_k)
}



target_function_g=function(ab, k, gk, gp, x_k, Y_k, weights, z, ind_k, mu_pres, sigma_pres, q){
  N_k=dim(x_k)[1]
  p=dim(x_k)[2]
  n=length(z)
  a_j=ab[1]
  b_j=0
  gp[k]=gk
  g_j=gp
  #  g_j[1] = -sum(g_j[-1])
  
  q_k=q[ind_k, ]
  weight_k=weights[ind_k]
  
  #Evaluate log-likelihood
  Theta=tran(z, x_k, mu_pres, sigma_pres)
  #mean.v=apply(x*(rep(1, N)%*%t(mu)), 1, sum)
  #var.v=apply(x*(rep(1, N)%*%t(sigma^2)), 1, sum)
  #sd.v=var.v^0.5
  #t1=-0.5*((Theta-mean.v%*%t(rep(1,n)))/(sd.v%*%t(rep(1,n))))^2 #normal term in log-likelihood  #remove prior term
  temp=Theta*a_j+b_j*(rep(1,N_k)%*%t(rep(1, n)))+(x_k%*%g_j)%*%t(rep(1, n)) #N_k by n
  t2=Y_k%*%t(rep(1, n))*temp-log(1+exp(temp)) #The remaining term in log-likelihood
  log_likelihood=(weight_k%*%t(rep(1, n)))*t2 #N_k by n
  res=-sum(q_k*log_likelihood, na.rm = TRUE)
  return(res/N_k)
}



target_mean=function(para_k, k, x, weights, z, mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  
  #  mu_pres[1]=0
  me=mu_pres
  me[k]=para_k
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu_pres, sigma_pres)
  mean.v=apply(x*(rep(1, N)%*%t(me)), 1, sum, na.rm = TRUE)
  var.v=apply(x*(rep(1, N)%*%t(sigma_pres^2)), 1, sum, na.rm = TRUE)
  #zero.vec=apply(x, 1, if_zero)
  #var.v=var.v+zero.vec
  sd.v=var.v^0.5
  t1=-0.5*((Theta-mean.v%*%t(rep(1,n)))/(sd.v%*%t(rep(1,n))))^2 #normal term in log-likelihood  #remove prior term
  log_likelihood=(weights%*%t(rep(1,n)))*t1 #N by n
  res=-sum(q*log_likelihood, na.rm = TRUE)
  return(res/N)
}

target_var=function(para_k, k, x, weights, z, mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  
  sig=sigma_pres
  sig[k]=para_k
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu_pres, sigma_pres)
  mean.v=apply(x*(rep(1, N)%*%t(mu_pres)), 1, sum, na.rm = TRUE)
  var.v=apply(x*(rep(1, N)%*%t(sig^2)), 1, sum, na.rm = TRUE)
  #zero.vec=apply(x, 1, if_zero)
  #var.v=var.v+zero.vec
  sd.v=var.v^0.5
  t1=-0.5*((Theta-mean.v%*%t(rep(1,n)))/(sd.v%*%t(rep(1,n))))^2 #N by n
  t2=-0.5*log(2*pi*(var.v)) #N dim vector
  res=sum(t2*weights, na.rm = TRUE)+ sum(q*t1*(weights%*%t(rep(1,n))), na.rm = TRUE)
  return(-res/N)
}





#Evaluate gradient with respect to phi
grad_phi1=function(phi_j, x, Y_j, weights, z, mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  #a_j=phi_j[1]
  b_j=phi_j
  g_j=rep(0, p)
  a_j=1
  
  #Evaluate gradient
  Theta=tran(z, x, mu_pres, sigma_pres)
  temp=Theta*a_j+b_j*(rep(1,N)%*%t(rep(1, n)))+(x%*%g_j)%*%t(rep(1, n)) #N by n
  prob=exp(temp)/(1+exp(temp))
  
  da=sum(Theta*(Y_j%*%t(rep(1, n))-prob)*q*(weights%*%t(rep(1,n))), na.rm = TRUE)
  db=sum((Y_j%*%t(rep(1, n))-prob)*q*(weights%*%t(rep(1,n))), na.rm = TRUE)
  return(-db/N)
}


#Evaluate gradient with respect to a_j, b_j for j not anchor items.
grad_phi_ab=function(ab, g, x, Y_j, weights, z, mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  a_j=ab
  b_j=0
  g_j=g
  # g_j[1]=0
  
  #Evaluate gradient
  Theta=tran(z, x, mu_pres, sigma_pres)
  #mean.v=apply(x*(rep(1, N)%*%t(mu)), 1, sum)
  #var.v=apply(x*(rep(1, N)%*%t(sigma^2)), 1, sum)
  #sd.v=var.v^0.5
  temp=Theta*a_j+b_j*(rep(1,N)%*%t(rep(1, n)))+(x%*%g_j)%*%t(rep(1, n)) #N by n
  prob=exp(temp)/(1+exp(temp))
  
  da=sum(Theta*(Y_j%*%t(rep(1, n))-prob)*q*(weights%*%t(rep(1,n))), na.rm = TRUE)
  #db=sum((Y_j%*%t(rep(1, n))-prob)*q*(weights%*%t(rep(1,n))), na.rm = TRUE)
  return(-da/N)
}


#Evaluate gradient with respect to g_k
grad_phi_g1=function(ab, k, gk, gp, x_k, Y_k, weights, z, ind_k,  mu_pres, sigma_pres, q){
  N_k=dim(x_k)[1]
  p=dim(x_k)[2]
  n=length(z)
  a_j=1
  b_j=0
  gp[k]=gk
  g_j=gp
  #  g_j[1] = -sum(g_j[-1])
  
  q_k=q[ind_k, ]
  weight_k=weights[ind_k]
  
  #Evaluate gradient
  Theta=tran(z, x_k, mu_pres, sigma_pres)
  #mean.v=apply(x*(rep(1, N)%*%t(mu)), 1, sum)
  #var.v=apply(x*(rep(1, N)%*%t(sigma^2)), 1, sum)
  #sd.v=var.v^0.5
  temp=Theta*a_j+b_j*(rep(1,N_k)%*%t(rep(1, n)))+(x_k%*%g_j)%*%t(rep(1, n)) #N by n
  prob=exp(temp)/(1+exp(temp))
  dg=sum((Y_k%*%t(rep(1, n))-prob)*q_k*(weight_k%*%t(rep(1, n))), na.rm = TRUE)
  return(-dg/N_k)
}


grad_phi_g=function(ab, k, gk, gp, x_k, Y_k, weights, z, ind_k,  mu_pres, sigma_pres, q){
  N_k=dim(x_k)[1]
  p=dim(x_k)[2]
  n=length(z)
  a_j=ab[1]
  b_j=0
  gp[k]=gk
  g_j=gp
  #  g_j[1] = -sum(g_j[-1])
  
  q_k=q[ind_k, ]
  weight_k=weights[ind_k]
  
  #Evaluate gradient
  Theta=tran(z, x_k, mu_pres, sigma_pres)
  #mean.v=apply(x*(rep(1, N)%*%t(mu)), 1, sum)
  #var.v=apply(x*(rep(1, N)%*%t(sigma^2)), 1, sum)
  #sd.v=var.v^0.5
  temp=Theta*a_j+b_j*(rep(1,N_k)%*%t(rep(1, n)))+(x_k%*%g_j)%*%t(rep(1, n)) #N by n
  prob=exp(temp)/(1+exp(temp))
  dg=sum((Y_k%*%t(rep(1, n))-prob)*q_k*(weight_k%*%t(rep(1, n))), na.rm = TRUE)
  return(-dg/N_k)
}




grad_mean=function(para_k, k, x, weights, z, mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  
  #  mu_pres[1]=-sum(mu_pres[-1])
  me=mu_pres
  me[k]=para_k
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu_pres, sigma_pres)
  mean.v=apply(x*(rep(1, N)%*%t(me)), 1, sum, na.rm = TRUE)
  var.v=apply(x*(rep(1, N)%*%t(sigma_pres^2)), 1, sum, na.rm = TRUE)
  # zero.vec=apply(x, 1, if_zero)
  # var.v=var.v+zero.vec
  dmu_k=sum(((x[,k]/var.v)%*%t(rep(1,n)))*(Theta-mean.v%*%t(rep(1, n)))*q*(weights%*%t(rep(1,n))), na.rm = TRUE)
  return(-dmu_k/N)
}


grad_var=function(para_k, k, x, weights, z, mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  
  sig=sigma_pres
  sig[k]=para_k
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu_pres, sigma_pres)
  mean.v=apply(x*(rep(1, N)%*%t(mu_pres)), 1, sum, na.rm = TRUE)
  var.v=apply(x*(rep(1, N)%*%t(sig^2)), 1, sum, na.rm = TRUE)
  # zero.vec=apply(x, 1, if_zero)
  # var.v=var.v+zero.vec
  t1=sum(-weights*x[,k]/2/var.v, na.rm = TRUE)
  t2=sum(((x[,k]%*%t(rep(1,n)))*q*(weights%*%t(rep(1,n)))*(Theta-mean.v%*%t(rep(1,n)))^2)/2/((var.v%*%t(rep(1,n)))^2), na.rm = TRUE)
  dsig_k=t1+t2
  return(-dsig_k/N)
}


#Evaluate marginal log-likelihood
mml=function(a, b, g, x, dat, weights, z, mu, sigma, w){ #q is quadrature
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  J=dim(dat)[2]
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu, sigma) #N by n
  th=Theta[, 1]
  temp=(th%*%t(a))+rep(1, N)%*%t(b)+x%*%t(g) #N by J
  prob=exp(temp)/(1+exp(temp))
  res=apply((prob^dat)*((1-prob)^(1-dat)), 1, prod, na.rm=TRUE)*w[1] ## summation of log instead
  for (k in 2:n){
    th=Theta[, k]
    temp=(th%*%t(a))+rep(1, N)%*%t(b)+x%*%t(g) #N by J
    prob=exp(temp)/(1+exp(temp))
    res_r=apply((prob^dat)*((1-prob)^(1-dat)), 1, prod, na.rm=TRUE)*w[k]
    res=res+res_r
  }
  res=-sum((weights%*%t(rep(1, J)))*log(res*(pi)^(-0.5)))
  return(res)
}



#weights are sampling weights
EM_2PL_inference <- function(a, b, g, x, dat, weights, z, mu, sigma, w, ite, tol = 0.001){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  J=dim(dat)[2]
  phi0=matrix(0, J, 2+p)
  phi1=cbind(a, b, g)
  mu0=rep(0, p)
  mu1=mu
  sigma0=rep(0, p)
  sigma1=sigma
  ind_group=list()
  dat_group=list()
  x_group=list()
  for (k in 1:p){
    ind_c=which(x[,k]==1)
    ind_group[[k]]=ind_c
    dat_group[[k]]=dat[ind_c,]
    x_group[[k]]=x[ind_c,]
  }
  #MLL0=0
  #MLL = mml(a, b, g, x, dat, z, mu, sigma, w)
  while(max(c(abs(phi0-phi1)), abs(mu0-mu1), abs(sigma0-sigma1))>tol){#max
    #MLL0=MLL
    phi0=phi1
    mu0=mu1
    sigma0=sigma1
    
    #E-step
    quadrature=quad(a=phi0[,1], b=phi0[,2], g=phi0[ ,c(3:(2+p))], x, dat, z, mu0, sigma0, w)
    
    #M-step
    ###Add the sampling weights to optim Y_j
    phi1[1, 1]=0.5170212
    mu1 = rep(0, p)
    phi1[, 2] = rep(0, J)
    
    for(j in 1:J){#update phi_j one by one
      if(j==1){
        for (k in 1:p){
          par_updates_g = optim(par=phi1[j, 2+k], fn=target_function_g1, gr=grad_phi_g1, method = "L-BFGS-B", ab=phi1[j, c(1,2)], k=k, gp=phi1[j, 3:(2+p)], x_k=x_group[[k]], Y_k=dat_group[[k]][,j],ind_k=ind_group[[k]], weights=weights, z=z, mu_pres=mu1,  sigma_pres=sigma1, q=quadrature, lower=-5, upper=5)
          phi1[j, 2+k]=par_updates_g$par
        }
      }
      else{
        for (l in 1:ite){
          par_updates_ab = optim(par=phi1[j, 1], fn=target_function_ab, gr=grad_phi_ab, method = "L-BFGS-B", g=phi1[j, 3:(2+p)], x=x, Y_j=dat[,j], weights=weights, z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=-5, upper=5)
          phi1[j, 1]=par_updates_ab$par
          for (k in 1:p){
            par_updates_g = optim(par=phi1[j, 2+k], fn=target_function_g, gr=grad_phi_g, method = "L-BFGS-B", ab=phi1[j, c(1,2)], k=k, gp=phi1[j, 3:(2+p)], x_k=x_group[[k]], Y_k=dat_group[[k]][,j], weights=weights, ind_k=ind_group[[k]], z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=-5, upper=5)
            phi1[j, 2+k]=par_updates_g$par
          }
        }
      }
    }
    #ab, k, gk, gp, x_k, Y_k, z, ind_k,  mu_pres, sigma_pres, q
    #To update mu one by one
    # for (k in 2:p){
    #   par_updates = optim(par=mu0[k], fn=target_mean, gr=grad_mean, method = "L-BFGS-B", k=k,x=x, z=z, weights=weights, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=-5, upper=5)
    #   mu1[k]=par_updates$par
    # }
    
    #To update sigma one by one
    for (k in 1:p){
      par_updates = optim(par=sigma0[k], fn=target_var, gr=grad_var, method = "L-BFGS-B", k=k,x=x, z=z, weights=weights, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=0.1, upper=5)
      sigma1[k]=par_updates$par
    }
    
    
    MLL = mml(a=phi1[,1], b=phi1[,2], g=phi1[, c(3:(2+p))], x, dat, weights, z, mu=mu1, sigma=sigma1, w)
    print(MLL) 
  }
  list(mu=mu1, sigma=sigma1, alpha.vec=phi1[, 1], beta.vec=phi1[, 2], gamma.vec = phi1[, c(3:(p+2))], post=quadrature);
}


set.seed(3)


#Read in data
dat=read.csv(file = 'PISA_2022_science_sas.csv',  sep = ",", header = T)
dat=dat[,-1]
colnames(dat)=NULL
dat=data.matrix(dat)

x=read.csv(file = 'countries_2022_science_sas.csv', sep = ",", header = T)
x=x[,-1]
colnames(x)=NULL
x=data.matrix(x)
covariates=x
weights=read.csv('PISA_2022_science_samplingweights_sas.csv', header = T)[,1]


N=dim(dat)[1]
J=dim(dat)[2]
p=dim(x)[2]


a_ini=read.csv('a_ini_science_2022.csv', header=F)[,1]
b_ini=rep(0, J)
g_ini=matrix(0, J, p)
mu_ini=rep(0, p)
sigma_ini=runif(p, 0.2, 1)


a=a_ini
b=b_ini
g=g_ini
mu=mu_ini
sigma=sigma_ini
r=EM_2PL_inference(a=a, b=b, g=g, x=x, dat=dat, weights=weights, z, mu=mu, sigma=sigma, w, ite=1, tol = 0.05)

g_hat=r$gamma.vec
mu_hat=r$mu
sigma_hat=r$sigma
a_hat=r$alpha.vec
b_hat=r$beta.vec

write.table(g_hat, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamma_science_2022_d0.csv')
write.table(a_hat, sep=",",  col.names=FALSE, row.names=FALSE, file = 'alpha_science_2022_d0.csv')
write.table(b_hat, sep=",",  col.names=FALSE, row.names=FALSE, file = 'beta_science_2022_d0.csv')
write.table(mu_hat, sep=",",  col.names=FALSE, row.names=FALSE, file = 'mu_science_2022_d0.csv')
write.table(sigma_hat, sep=",",  col.names=FALSE, row.names=FALSE, file = 'sigma_science_2022_d0.csv')
