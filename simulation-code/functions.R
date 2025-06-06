#install.packages("Rcpp")
library(Rcpp)
#install.packages("fastGHQuad")
library(fastGHQuad)
library(quantreg)
library(MASS)
library(mirt)

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
  Z=rowSums(q*(rep(1, N)%*%t(w)), na.rm = TRUE)
  return((q*(rep(1, N)%*%t(w)))/(Z%*%t(rep(1, n))))
}






#Evaluate target function at j=1 (i.e. anchor item), phi_j=(a_j, b_j, g_j1,...,g_jp), i.e. negative log-likelihood at phi_j
#Y_j=Y[,j]
target_function1=function(phi_j, x, Y_j, z,  mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  #a_j=phi_j[1]
  b_j=phi_j
  g_j=rep(0, p)
  a_j=1
  
  
  #Evaluate parameter dependent log-likelihood
  Theta=tran(z, x, mu_pres, sigma_pres)
  temp=Theta*a_j+b_j*(rep(1,N)%*%t(rep(1, n)))+(x%*%g_j)%*%t(rep(1, n)) #N by n
  t2=Y_j%*%t(rep(1, n))*temp-log(1+exp(temp)) #The parameter dependent term in log-likelihood
  res=-sum(q*t2, na.rm = TRUE)
  return(res/N)
}

#Evaluate target function at phi_j=(a_j, b_j) given (g_j1, ... g_jp) at previous step values.
#Y_j=Y[,j]
target_function_ab=function(ab, g, x, Y_j, z,  mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  a_j=ab
  b_j=0
  g_j=g
 
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu_pres, sigma_pres)
  temp=Theta*a_j+b_j*(rep(1,N)%*%t(rep(1, n)))+(x%*%g_j)%*%t(rep(1, n)) #N by n
  t2=Y_j%*%t(rep(1, n))*temp-log(1+exp(temp)) #The remaining term in log-likelihood
  res=-sum(q*t2, na.rm = TRUE)
  return(res/N)
}

#Evaluate target function at g_jk, given (a_j, b_j, g_j1,...g_j(k-1),g_j(k+1),...g_jp) at previous step values,
#Y_k=Y[[i in group k],j]
#ab=(a_j, b_j)
#gp=(g_j1,...g_j(k-1))
#ga=(g_j(k+1),...g_jp)
#gk=g[k]
#ind_k is the indices of individuals belonging to group k
target_function_g1=function(ab, k, gk, gp, x_k, Y_k, z, ind_k, mu_pres, sigma_pres, q){
  N_k=dim(x_k)[1]
  p=dim(x_k)[2]
  n=length(z)
  a_j=1
  b_j=0
  gp[k]=gk
  g_j=gp
  q_k=q[ind_k, ]
  
  #Evaluate parameter dependent log-likelihood
  Theta=tran(z, x_k, mu_pres, sigma_pres)
  temp=Theta*a_j+b_j*(rep(1,N_k)%*%t(rep(1, n)))+(x_k%*%g_j)%*%t(rep(1, n)) #N_k by n
  t2=Y_k%*%t(rep(1, n))*temp-log(1+exp(temp)) #The parameter dependent term in log-likelihood
  res=-sum(q_k*t2, na.rm = TRUE)
  return(res/N_k)
}


target_function_g=function(ab, k, gk, gp, x_k, Y_k, z, ind_k, mu_pres, sigma_pres, q){
  N_k=dim(x_k)[1]
  p=dim(x_k)[2]
  n=length(z)
  a_j=ab[1]
  b_j=0
  gp[k]=gk
  g_j=gp
 
  
  q_k=q[ind_k, ]
  
  #Evaluate parameter dependent log-likelihood
  Theta=tran(z, x_k, mu_pres, sigma_pres)
  temp=Theta*a_j+b_j*(rep(1,N_k)%*%t(rep(1, n)))+(x_k%*%g_j)%*%t(rep(1, n)) #N_k by n
  t2=Y_k%*%t(rep(1, n))*temp-log(1+exp(temp)) #The parameter dependent term in log-likelihood
  res=-sum(q_k*t2, na.rm = TRUE)
  return(res/N_k)
}



target_mean=function(para_k, k, x, z, mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  
  mu_pres[1]=0
  me=mu_pres
  me[k]=para_k
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu_pres, sigma_pres)
  mean.v=apply(x*(rep(1, N)%*%t(me)), 1, sum)
  var.v=apply(x*(rep(1, N)%*%t(sigma_pres^2)), 1, sum)
  sd.v=var.v^0.5
  t1=-0.5*((Theta-mean.v%*%t(rep(1,n)))/(sd.v%*%t(rep(1,n))))^2 #normal term in log-likelihood  #remove prior term
  log_likelihood=t1 #N by n
  res=-sum(q*log_likelihood, na.rm = TRUE)
  return(res/N)
}

target_var=function(para_k, k, x, z, mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  
  #sigma_pres[1]=1
  sig=sigma_pres
  sig[k]=para_k
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu_pres, sigma_pres)
  mean.v=apply(x*(rep(1, N)%*%t(mu_pres)), 1, sum)
  var.v=apply(x*(rep(1, N)%*%t(sig^2)), 1, sum)
  sd.v=var.v^0.5
  t1=-0.5*((Theta-mean.v%*%t(rep(1,n)))/(sd.v%*%t(rep(1,n))))^2 #normal term in log-likelihood  #remove prior term
  t2=sum(-0.5*log(2*pi*(var.v)), na.rm = TRUE)
  log_likelihood=t1 #N by n
  res=-sum(q*t1)-t2
  return(res/N)
}





#Evaluate gradient with respect to phi
grad_phi1=function(phi_j, x, Y_j, z, mu_pres, sigma_pres, q){
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
  
  da=sum(Theta*(Y_j%*%t(rep(1, n))-prob)*q, na.rm = TRUE)
  db=sum((Y_j%*%t(rep(1, n))-prob)*q, na.rm = TRUE)
  return(-db/N)
}


#Evaluate gradient with respect to a_j, b_j for j not anchor items.
grad_phi_ab=function(ab, g, x, Y_j, z, mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  a_j=ab
  b_j=0
  g_j=g

  
  #Evaluate gradient
  Theta=tran(z, x, mu_pres, sigma_pres)
  temp=Theta*a_j+b_j*(rep(1,N)%*%t(rep(1, n)))+(x%*%g_j)%*%t(rep(1, n)) #N by n
  prob=exp(temp)/(1+exp(temp))
  
  da=sum(Theta*(Y_j%*%t(rep(1, n))-prob)*q, na.rm = TRUE)
#  db=sum((Y_j%*%t(rep(1, n))-prob)*q, na.rm = TRUE)
  return(-da/N)
}


#Evaluate gradient with respect to g_k
grad_phi_g1=function(ab, k, gk, gp, x_k, Y_k, z, ind_k,  mu_pres, sigma_pres, q){
  N_k=dim(x_k)[1]
  p=dim(x_k)[2]
  n=length(z)
  a_j=1
  b_j=0
  gp[k]=gk
  g_j=gp

  
  q_k=q[ind_k, ]
  
  #Evaluate gradient
  Theta=tran(z, x_k, mu_pres, sigma_pres)
  temp=Theta*a_j+b_j*(rep(1,N_k)%*%t(rep(1, n)))+(x_k%*%g_j)%*%t(rep(1, n)) #N by n
  prob=exp(temp)/(1+exp(temp))
  dg=sum((Y_k%*%t(rep(1, n))-prob)*q_k, na.rm = TRUE)
  return(-dg/N_k)
}

grad_phi_g=function(ab, k, gk, gp, x_k, Y_k, z, ind_k,  mu_pres, sigma_pres, q){
  N_k=dim(x_k)[1]
  p=dim(x_k)[2]
  n=length(z)
  a_j=ab[1]
  b_j=0
  gp[k]=gk
  g_j=gp

  
  q_k=q[ind_k, ]
  
  #Evaluate gradient
  Theta=tran(z, x_k, mu_pres, sigma_pres)
  temp=Theta*a_j+b_j*(rep(1,N_k)%*%t(rep(1, n)))+(x_k%*%g_j)%*%t(rep(1, n)) #N by n
  prob=exp(temp)/(1+exp(temp))
  dg=sum((Y_k%*%t(rep(1, n))-prob)*q_k, na.rm = TRUE)
  return(-dg/N_k)
}




grad_mean=function(para_k, k, x, z, mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  
  mu_pres[1]=0
  me=mu_pres
  me[k]=para_k
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu_pres, sigma_pres)
  mean.v=apply(x*(rep(1, N)%*%t(me)), 1, sum)
  var.v=apply(x*(rep(1, N)%*%t(sigma_pres^2)), 1, sum)
  #zero.vec=apply(x, 1, if_zero)
  #var.v=var.v+zero.vec
  dmu_k=sum(((x[,k]/var.v)%*%t(rep(1,n)))*(Theta-mean.v%*%t(rep(1, n)))*q, na.rm = TRUE)
  return(-dmu_k/N)
}


grad_var=function(para_k, k, x, z, mu_pres, sigma_pres, q){
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  
  #sigma_pres[1]=1
  sig=sigma_pres
  sig[k]=para_k
  
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu_pres, sigma_pres)
  mean.v=apply(x*(rep(1, N)%*%t(mu_pres)), 1, sum)
  var.v=apply(x*(rep(1, N)%*%t(sig^2)), 1, sum)
  #zero.vec=apply(x, 1, if_zero)
  #var.v=var.v+zero.vec
  t1=sum(-x[,k]/2/var.v, na.rm = TRUE)
  t2=sum(((x[,k]%*%t(rep(1,n)))*q*(Theta-mean.v%*%t(rep(1,n)))^2)/2/((var.v%*%t(rep(1,n)))^2), na.rm = TRUE)
  dsig_k=t1+t2
  return(-dsig_k/N)
}


#Evaluate marginal log-likelihood
mml=function(a, b, g, x, dat, z, mu, sigma, w){ #q is quadrature
  N=dim(x)[1]
  p=dim(x)[2]
  n=length(z)
  J=dim(dat)[2]
  
  #Evaluate log-likelihood
  Theta=tran(z, x, mu, sigma) #N by n
  th=Theta[, 1]
  temp=(th%*%t(a))+rep(1, N)%*%t(b)+x%*%t(g) #N by J
  prob=exp(temp)/(1+exp(temp))
  res=apply((prob^dat)*((1-prob)^(1-dat)), 1, prod, na.rm=TRUE)*w[1]
  for (k in 2:n){
    th=Theta[, k]
    temp=(th%*%t(a))+rep(1, N)%*%t(b)+x%*%t(g) #N by J
    prob=exp(temp)/(1+exp(temp))
    res_r=apply((prob^dat)*((1-prob)^(1-dat)), 1, prod, na.rm=TRUE)*w[k]
    res=res+res_r
  }
  res=-sum(log(res*(pi)^(-0.5)))
  return(res)
}


EM_2PL_inference <- function(a, b, g, x, dat, z, mu, sigma, w, ite, tol = 0.001){
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
  while(max(c(abs(phi0-phi1)), abs(mu0-mu1), abs(sigma0-sigma1))>tol){#max
    phi0=phi1
    mu0=mu1
    sigma0=sigma1
    
    #E-step
    quadrature=quad(a=phi0[,1], b=phi0[,2], g=phi0[,c(3:(2+p))], x, dat, z, mu0, sigma0, w)
    
    #M-step
    phi1[1, 1]=1
    mu1 = rep(0, p)
    phi1[, 2] = rep(0, J)
   
    
    for(j in 1:J){#update phi_j one by one
      if(j==1){
        #par_updates = optim(par=phi1[1, 2], fn=target_function1, gr=grad_phi1, method = "L-BFGS-B", x=x, Y_j=dat[,1], z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=-5, upper=5)
        for (k in 1:p){
          par_updates_g = optim(par=phi1[j, 2+k], fn=target_function_g1, gr=grad_phi_g1, method = "L-BFGS-B", ab=phi1[j, c(1,2)], k=k, gp=phi1[j, 3:(2+p)], x_k=x_group[[k]], Y_k=dat_group[[k]][,j],ind_k=ind_group[[k]], z=z, mu_pres=mu1,  sigma_pres=sigma1, q=quadrature, lower=-5, upper=5)
          phi1[j, 2+k]=par_updates_g$par
        }
        
      }
      else{
        for (l in 1:ite){
          par_updates_ab = optim(par=phi1[j, 1], fn=target_function_ab, gr=grad_phi_ab, method = "L-BFGS-B", g=phi1[j, 3:(2+p)], x=x, Y_j=dat[,j], z=z, mu_pres=mu1,  sigma_pres=sigma1, q=quadrature, lower=-5, upper=5)
          phi1[j, 1]=par_updates_ab$par
          for (k in 1:p){
            par_updates_g = optim(par=phi1[j, 2+k], fn=target_function_g, gr=grad_phi_g, method = "L-BFGS-B", ab=phi1[j, c(1,2)], k=k, gp=phi1[j, 3:(2+p)], x_k=x_group[[k]], Y_k=dat_group[[k]][,j],ind_k=ind_group[[k]], z=z, mu_pres=mu1,  sigma_pres=sigma1, q=quadrature, lower=-5, upper=5)
            phi1[j, 2+k]=par_updates_g$par
          }
        }
      }
    }
    
    
    #ab, k, gk, gp, x_k, Y_k, z, ind_k,  mu_pres, sigma_pres, q
    
    
    #To update mu one by one
    #for (k in 2:p){
    #  par_updates = optim(par=mu0[k], fn=target_mean, gr=grad_mean, method = "L-BFGS-B", k=k,x=x, z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=-5, upper=5)
    #  mu1[k]=par_updates$par
    #}
    
    #To update sigma one by one
    for (k in 1:p){
      par_updates = optim(par=sigma0[k], fn=target_var, gr=grad_var, method = "L-BFGS-B", k=k,x=x, z=z, mu_pres=mu0,  sigma_pres=sigma0, q=quadrature, lower=0.1, upper=5)
      sigma1[k]=par_updates$par
    }
    
    MLL = mml(a=phi1[,1], b=phi1[,2], g=phi1[, c(3:(2+p))], x, dat, z, mu=mu1, sigma=sigma1, w)
    print(MLL) 
  }
  list(mu=mu1, sigma=sigma1, alpha.vec=phi1[, 1], beta.vec=phi1[, 2], gamma.vec = phi1[, c(3:(p+2))], post=quadrature);
}



# Define the loss function
loss_function <- function(params, gamma, a, b, J, p) {
  c <- params[1:p]          # Extract c_k parameters
  d <- params[(p+1):(p+J)]  # Extract d_j parameters
  
  total_loss <- 0
  for (k in 1:p) {
    for (j in 1:J) {
      total_loss <- total_loss + abs(b[j] + gamma[j, k] + c[k] * a[j] - d[j])
    }
  }
  
  return(total_loss)
}



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


eva_MSE=function(est, tru){#tru is a vector of the true parameter, est is a matrix of estimated parameter
  N=dim(est)[1]
  squared_err=NULL
  for (i in 1:N){
    vec=(as.numeric(est[i,])-tru)^2
    squared_err=rbind(squared_err, vec)
  }
  #true_mat=rep(1, N)%*%t(tru)
  #est=as.matrix(est)
  #squared_error_mat=(est-tru)^2
  mse=apply(squared_err, 2, mean)
  return(as.numeric(mse))
  #return(squared_err)
}

eva_MSE_gamma=function(est, tru, J=15){#tru is a vector of the true parameter, est is a matrix of estimated parameter
  p=dim(est)[2]
  N=dim(est)[1]/J
  squared_err=NULL
  
  for (i in 1:N){
    g_mtx = est[((i-1)*J+1):(i*J),]
    err = mean((c(t(g_mtx)) - tru)^2)
    
    vec=(as.numeric(c(t(g_mtx)))-tru)^2
    squared_err=rbind(squared_err, err)
  }
  mse=mean(squared_err)
  return(as.numeric(mse))
  #return(squared_err)
}




coverage_proposed_inference=function(seed, N, sig, beta.vec,  alpha.vec, gamma, mu, sigma, x, ite,  sig_list){
  set.seed(seed)
  J=length(alpha.vec)
  p=dim(x)[2]
  n_sig=length(sig_list)
  gamm=c(t(gamma))
  ind1_gam=which(gamm!=0)
  ind0_gam=which(gamm==0)
  n_p_gam=length(ind1_gam)
  n_n_gam=length(ind0_gam)
  
  cover_gam=NULL
  p_value_gam=NULL
  alpha=NULL
  beta=NULL
  rgamma=NULL
  TPR_gam=NULL
  FPR_gam=NULL
  r_mu=NULL
  r_sigma=NULL
  p_val_hotelling=NULL
  
  tran_mu = NULL
  tran_g = NULL
  tran_d = NULL
  for (l in 1:10){
    print(l)
    set.seed(seed*10 + l)
  
    #Generate data
    mu.v=apply((x*rep(1, N))%*%(mu), 1, sum)
    var.v=apply((x*(rep(1, N)))%*%(sigma^2), 1, sum)
    sigma.v=var.v^0.5
    th=rnorm(n=N,mean=mu.v, sd=sigma.v)
    temp=(th%*%t(rep(1, J)))*(rep(1, N)%*%t(alpha.vec))+rep(1, N)%*%t(beta.vec)+x%*%t(gamma) #N by J
    prob=exp(temp)/(1+exp(temp)) #dim N by J
    dat = matrix(0, N, J);
    dat[] = rbinom(N*J, 1, prob);
   
    r=EM_2PL_inference(a=alpha.vec,  b=beta.vec, g=gamma, x=x, dat=dat, z, mu=mu, sigma=sigma, w,ite,  tol = 0.01)

    g_hat=r$gamma.vec
    mu_hat=r$mu
    sigma_hat=r$sigma
    a_hat=r$alpha.vec
    b_hat=r$beta.vec
    
    r_mu=rbind(r_mu, mu_hat)
    r_sigma=rbind(r_sigma, sigma_hat)
    alpha=rbind(alpha, a_hat)
    beta=rbind(beta, b_hat)
    rgamma=rbind(rgamma, g_hat)

    gam=r$gamma.vec
    opt_quantile_c1 = solve_optimization_quantile_c1(gam, r$alpha.vec, b_hat)
    opt_d = opt_quantile_c1$d
    opt_c = c(0, opt_quantile_c1$c)
    mean_c = mean(opt_c)
    opt_c = opt_c - mean_c
    opt_d = opt_d + mean_c*r$alpha.vec
      
    transform_p = gam - outer(r$alpha.vec, opt_c) + b_hat - opt_d
    transform_mu = mu_hat + opt_c
    transform_d = opt_d

    tran_mu = rbind(tran_mu, transform_mu)
    tran_g = rbind(tran_g, transform_p)
    tran_d = rbind(tran_d, transform_d)
    
    
    transform_p=c(t(transform_p)) #J*p dim vector
    
  }
  
   g_MSE=mean(eva_MSE_gamma(est=tran_g, tru=c(t(gamma))))
   a_MSE=mean(eva_MSE(est=alpha, tru=alpha.vec))
   b_MSE=mean(eva_MSE(est=tran_d, tru=beta.vec))
   mu_MSE=mean(eva_MSE(est=tran_mu, tru=mu))
   sigma_MSE=mean(eva_MSE(est=r_sigma, tru=sigma))
  
  return(list( alpha=alpha, beta=beta, gamma=rgamma, mu=r_mu, sigma=r_sigma,
              tran_mu=tran_mu, tran_g=tran_g, tran_d=tran_d, g_MSE=g_MSE, a_MSE=a_MSE, 
              b_MSE=b_MSE, mu_MSE=mu_MSE, sigma_MSE=sigma_MSE))
}



