library(vegan)
library(gplots)
library(MASS)
library(fastGHQuad)
library(quantreg)
library(tidyverse)
library(sf)
library(rvest)
library(stringr)
library(scales)
library(xtable)
library(haven)
library(ggplot2)


##To read in the sampling weights
df_weights=read_sas('STU_QQQ_SAS/CY08MSP_STU_QQQ.SAS7BDAT')
sampling_weights=df_weights$W_FSTUWT

##To read in the cognitive data file
df=read_sas('STU_COG_SAS/CY08MSP_STU_COG.SAS7BDAT')

#Double-check the individual indices matches each other
student_id_questionnaire=df_weights$CNTSTUID
student_id_cog=df$CNTSTUID
all(student_id_questionnaire==student_id_cog) #true

#Re-arrange the sampling weights according to the student id in cog file
ind_weights=match(student_id_cog, student_id_questionnaire)
all(student_id_questionnaire[ind_weights]==student_id_cog) #true
sampling_weights=sampling_weights[ind_weights]

#Join the sampling weights to the response data file
df=cbind.data.frame(df, sampling_weights)

temp_science <- read.csv("science_names.csv", sep=',')
temp_science = temp_science[,1]
item_science<- c()
for(item in temp_science){
  print(unique(df[, item]))
  if(setequal(unique(df[, item]), c(NA, 0, 1))){
    item_science = c(item_science, item)
  }
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



length(item_science)#103


g_hat <- as.matrix(read.csv("gamma_science_2022_d0.csv", header = F))
a_hat <- read.csv("alpha_science_2022_d0.csv", header = F)[,1]
b_hat <- read.csv("beta_science_2022_d0.csv", header = F)[,1]
mu_hat <- read.csv("mu_science_2022_d0.csv", header = F)[,1]
sigma_hat <- read.csv("sigma_science_2022_d0.csv", header = F)[,1]

g_hat = as.matrix(g_hat)
transform_p=NULL
mu_tran=c()

opt_quantile_c1 = solve_optimization_quantile_c1(g_hat, a_hat, b_hat)
opt_d = opt_quantile_c1$d
opt_c = c(0, opt_quantile_c1$c)
mean_c = mean(opt_c)
opt_c = opt_c - mean_c
opt_d = opt_d + mean_c*a_hat

transform_p = g_hat - outer(a_hat, opt_c) + b_hat - opt_d
transform_mu = mu_hat + opt_c
transform_d = opt_d

xtable(transform_p[,1:18], type = "latex", file = "compare_latent_skill.tex", digits=2)
xtable(transform_p[,19:37], type = "latex", file = "compare_latent_skill.tex", digits=2)


trans_a_matrix = matrix(a_hat, ncol=2, byrow=F)
output_a_matrix = cbind(1:52, trans_a_matrix[,1], 53:104,trans_a_matrix[,2])
xtable(as.matrix(output_a_matrix), type = "latex", file = "compare_latent_skill.tex", digits=2)


trans_d_matrix = matrix(transform_d, ncol=2, byrow=F)
output_d_matrix = cbind(1:52, trans_d_matrix[,1], 53:104,trans_d_matrix[,2])
xtable(as.matrix(output_d_matrix), type = "latex", file = "compare_latent_skill.tex", digits=2)

transform_p=c(t(transform_p)) #J*p dim vector

#write.table(transform_p, sep=",",  col.names=FALSE, row.names=FALSE, file = 'gamma_transformed_science_2022.csv')
#write.table(transform_mu, sep=",",  col.names=FALSE, row.names=FALSE, file = 'mu_transformed_science_2022.csv')

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

gam_val=read.csv('gamma_transformed_science_2022_d0.csv', header=F)
gam_val=gam_val[,1]


#Similarity matrix

l2_norm=function(vector){
  return((sum(vector^2))^0.5)
}
country=c('Australia', 'Austria', 'Belgium', 'Canada', 'Chile', 'Colombia', 'Costa Rica', 'Czech Republic', 'Denmark',        
          'Estonia',  'Finland', 'France', 'Germany', 'Greece', 'Hungary', 'Iceland', 'Ireland', 'Israel', 'Italy', 'Japan', 'Korea',          
          'Latvia', 'Lithuania', 'Mexico', 'Netherlands', 'New Zealand','Norway', 'Poland', 'Portugal', 'Slovak Republic', 'Slovenia', 'Spain',          
          'Sweden', 'Switzerland', 'Turkey', 'United Kingdom', 'United States')

g_est=matrix(gam_val, byrow = T, nrow=J, ncol=p)


#Plot similarity matrix
#Evaluate a distance matrix of dimension 37 by 37
#dist_australia=apply(g_est, 2, l2_norm)
distance_mat=matrix(0, 37, 37)
for (k1 in 1:p){
  for (k2 in 1:p){
    dis=l2_norm(g_est[,k1]-g_est[,k2])
    distance_mat[k1, k2]=dis
  }
}


#Proposed method
mu_hat_full=as.matrix(transform_mu)
sigma_hat_full=sigma_hat

theta_ord_full=order(mu_hat_full, decreasing = T)
country_ord_full=country[theta_ord_full]
country_ord_full

round(mu_hat_full[theta_ord_full],3)
round(sigma_hat_full[theta_ord_full], 3)

# [1] "Japan"           "Korea"           "Estonia"         "Czech Republic"  "Australia"       "Ireland"        
# [7] "Switzerland"     "Poland"          "New Zealand"     "Finland"         "Austria"         "Germany"        
# [13] "Canada"          "Belgium"         "Latvia"          "Sweden"          "United Kingdom"  "Spain"          
# [19] "United States"   "Portugal"        "Hungary"         "France"          "Lithuania"       "Italy"          
# [25] "Slovenia"        "Norway"          "Denmark"         "Netherlands"     "Turkey"          "Slovak Republic"
# [31] "Israel"          "Chile"           "Iceland"         "Greece"          "Colombia"        "Mexico"         
# [37] "Costa Rica"   


### RMSD 0.05

mu_hat_0.05=read.csv( 'mu_science_RMSD_2022_0.05.csv', header=F)[,1]
theta_ord_0.05=order(mu_hat_0.05, decreasing = T)
country_ord_0.05=country[theta_ord_0.05]
country_ord_0.05
sigma_hat_0.05=read.csv( 'sigma_science_RMSD_2022_0.05.csv', header=F)[,1]

round(mu_hat_0.05[theta_ord_0.05],3)
round(sigma_hat_0.05[theta_ord_0.05], 3)


xtable(cbind(country_ord_full, round(mu_hat_full[theta_ord_full],3), round(sigma_hat_full[theta_ord_full], 3), 1:37,
             country_ord_0.05, round(mu_hat_0.05[theta_ord_0.05],3), round(sigma_hat_0.05[theta_ord_0.05], 3)), 
       type = "latex", file = "compare_latent_skill.tex", digits=3)

# [1] "Japan"           "Korea"           "Estonia"         "Czech Republic"  "Australia"       "Poland"         
# [7] "Switzerland"     "Ireland"         "New Zealand"     "Austria"         "Finland"         "United Kingdom" 
# [13] "Canada"          "Germany"         "Belgium"         "Latvia"          "Sweden"          "United States"  
# [19] "Spain"           "Portugal"        "Hungary"         "Slovenia"        "France"          "Italy"          
# [25] "Lithuania"       "Netherlands"     "Norway"          "Denmark"         "Turkey"          "Israel"         
# [31] "Slovak Republic" "Chile"           "Iceland"         "Greece"          "Colombia"        "Mexico"         
# [37] "Costa Rica"     

#To evaluate the Kendall's rank correlations
cor(mu_hat_full, mu_hat_0.05, method = "kendall")
#0.9459459


#### RMSD 0.1 
mu_hat_0.1=read.csv( 'mu_science_RMSD_2022_0.1.csv', header=F)[,1]
sigma_hat_0.1=read.csv( 'sigma_science_RMSD_2022_0.1.csv', header=F)[,1]

theta_ord_0.1=order(mu_hat_0.1, decreasing = T)
country_ord_0.1=country[theta_ord_0.1]
country_ord_0.1


round(mu_hat_0.1[theta_ord_0.1],3)
round(sigma_hat_0.1[theta_ord_0.1], 3)


# [1] "Japan"           "Korea"           "Estonia"         "Czech Republic"  "Australia"       "Poland"         
# [7] "Ireland"         "Switzerland"     "New Zealand"     "Austria"         "Finland"         "Germany"        
# [13] "United Kingdom"  "Spain"           "Latvia"          "Sweden"          "Canada"          "Belgium"        
# [19] "United States"   "Portugal"        "Hungary"         "France"          "Italy"           "Norway"         
# [25] "Lithuania"       "Netherlands"     "Slovenia"        "Denmark"         "Israel"          "Turkey"         
# [31] "Chile"           "Slovak Republic" "Iceland"         "Greece"          "Colombia"        "Costa Rica"     
# [37] "Mexico"         

#### RMSD 0.15 ######
mu_hat_0.15=read.csv( 'mu_science_RMSD_2022_0.15.csv', header=F)[,1]
sigma_hat_0.15=read.csv( 'sigma_science_RMSD_2022_0.15.csv', header=F)[,1]

#mu_hat_2PL=c(mu_hat_2PL)
theta_ord_0.15=order(mu_hat_0.15, decreasing = T)
country_ord_0.15=country[theta_ord_0.15]
country_ord_0.15

round(mu_hat_0.15[theta_ord_0.15],3)
round(sigma_hat_0.15[theta_ord_0.15], 3)

# [1] "Japan"           "Korea"           "Estonia"         "Czech Republic"  "Australia"       "Poland"         
# [7] "Ireland"         "Switzerland"     "New Zealand"     "Austria"         "Finland"         "Sweden"         
# [13] "Germany"         "United Kingdom"  "Belgium"         "Canada"          "Latvia"          "Spain"          
# [19] "United States"   "Portugal"        "Hungary"         "France"          "Italy"           "Lithuania"      
# [25] "Slovenia"        "Norway"          "Netherlands"     "Denmark"         "Turkey"          "Israel"         
# [31] "Chile"           "Slovak Republic" "Iceland"         "Greece"          "Colombia"        "Mexico"         
# [37] "Costa Rica"     

## baseline 

mu_hat_baseline=read.csv( 'mu_science_2022_baseline.csv', header=F)[,1]
sigma_hat_baseline=read.csv( 'sigma_science_2022_baseline.csv', header=F)[,1]

theta_ord_baseline=order(mu_hat_baseline, decreasing = T)
country_ord_baseline=country[theta_ord_baseline]
country_ord_baseline

round(mu_hat_baseline[theta_ord_baseline],3)
round(sigma_hat_baseline[theta_ord_baseline],3)

# [1] "Japan"           "Korea"           "Estonia"         "Australia"       "Czech Republic"  "Poland"         
# [7] "Ireland"         "Switzerland"     "Austria"         "New Zealand"     "Finland"         "Germany"        
# [13] "United Kingdom"  "Sweden"          "Belgium"         "Canada"          "Spain"           "Latvia"         
# [19] "United States"   "Portugal"        "Hungary"         "France"          "Italy"           "Lithuania"      
# [25] "Norway"          "Netherlands"     "Slovenia"        "Denmark"         "Turkey"          "Israel"         
# [31] "Chile"           "Slovak Republic" "Iceland"         "Greece"          "Colombia"        "Mexico"         
# [37] "Costa Rica"       


xtable(cbind(country_ord_0.1, round(mu_hat_0.1[theta_ord_0.1],3), round(sigma_hat_0.1[theta_ord_0.1], 3), 1:37,
             country_ord_0.15, round(mu_hat_0.15[theta_ord_0.15],3), round(sigma_hat_0.15[theta_ord_0.15], 3), 1:37,
             country_ord_baseline, round(mu_hat_baseline[theta_ord_baseline],3), round(sigma_hat_baseline[theta_ord_baseline], 3)), 
       type = "latex", file = "compare_latent_skill.tex", digits=3)


library(mize)
cost_fun <- function(R, D) {
  diff2 <- (R - D) ^ 2
  sum(diff2) * 0.5
}

g_est=matrix(gam_val, byrow = T, nrow=J, ncol=p)

l2_norm=function(vector){
  return((sum(vector^2))^0.5)
}

country = c("AUS", "AUT", "BEL", "CAN", "CHL", "COL", "CRI", "CZE", "DNK", "EST", 
            "FIN", "FRA", "DEU", "GRC", "HUN", "ISL", "IRL", "ISR", "ITA", "JPN", 
            "KOR", "LVA", "LTU", "MEX", "NLD", "NZL", "NOR", "POL", "PRT", "SVK", 
            "SVN", "ESP", "SWE", "CHE", "TUR", "GBR", "USA")

distance_mat=matrix(0, 37, 37)
for (k1 in 1:p){
  for (k2 in 1:p){
    dis=l2_norm(g_est[,k1]-g_est[,k2])
    distance_mat[k1, k2]=dis
  }
}
dist_mat <- as.matrix(distance_mat)
colnames(dist_mat) <- country
rownames(dist_mat) <- country

cost_grad <- function(R, D, y) {
  K <- (R - D) / (D + 1.e-10)
  
  G <- matrix(nrow = nrow(y), ncol = ncol(y))
  
  for (i in 1:nrow(y)) {
    dyij <- sweep(-y, 2, -y[i, ])
    G[i, ] <- apply(dyij * K[, i], 2, sum)
  }
  
  as.vector(t(G)) * -2
}

mmds_fn <- function(par) {
  R <- as.matrix(dist_mat)
  y <- matrix(par, ncol = 2, byrow = TRUE)
  D <- as.matrix(stats::dist(y))
  
  cost_fun(R, D)
}

mmds_gr <- function(par) {
  R <- as.matrix(dist_mat)
  y <- matrix(par, ncol = 2, byrow = TRUE)
  D <- as.matrix(stats::dist(y))
  
  cost_grad(R, D, y)
}

set.seed(39)
ed0 <- runif(dim(dist_mat)[2] * 2)

res_euro <- mize(ed0, list(fn = mmds_fn, gr = mmds_gr), 
                 method = "L-BFGS", verbose = TRUE, 
                 grad_tol = 1e-8, check_conv_every = 10)

plot_mmds <- function(coords, dist, ...) {
  if (methods::is(coords, "numeric")) {
    coords <- matrix(coords, ncol = 2, byrow = T)
  }
  graphics::plot(coords, type = 'n')
  graphics::text(coords[, 1], coords[, 2], labels = country, ...)
}


plot_mmds(res_euro$par, dist_mat, cex = 0.5)

plot_mtx <- matrix(res_euro$par, ncol = 2, byrow = TRUE)


color=c('green', 'purple', 'purple', 'tan1', 'gray33','gray33', 'gray33', 'red', 'blue','red',
        'blue', 'purple', 'purple','purple', 'red', 'blue', 'purple', 'black', 'red',
        'brown', 'brown', 'red','red', 'gray33', 'purple', 'green', 'blue',
        'red', 'purple', 'red','red', 'red', 'blue', 'purple', 'black','purple', 'tan1')

index=1:37
index=index[-c(37, 32, 30, 9, 34,  13, 21, 25, 26, 27, 14, 2, 3, 4, 15, 16, 11, 10, 24, 29, 17, 6)]
par(cex=0.95, mar=c(0.1, 0.1, 0.1, 0.1))
plot(plot_mtx[,1], plot_mtx[,2], col=color, axes = TRUE, pch = 19)
text(plot_mtx[index, 1], plot_mtx[index, 2], country[index], pos=1, cex=0.6, offset=0.4)
index_left=c(24,  13, 2, 10, 4)
index_right=c( 34,  15, 30, 3, 27)
index_up=c(29,21, 37, 26, 32, 16,14, 11,17,6,9, 25)

text(x=plot_mtx[index_left,1], y=plot_mtx[index_left,2], labels=c(country[index_left]), pos=2, cex=0.6, offset=0.4)
text(x=plot_mtx[index_right,1], y=plot_mtx[index_right,2], labels=c(country[index_right]), pos=4, cex=0.6, offset=0.4)
text(x=plot_mtx[index_up,1], y=plot_mtx[index_up,2], labels=c(country[index_up]), pos=3, cex=0.6, offset=0.4)


legend("bottomright", legend = c('Nordic', 'West Europe', 'Central&Eastern&Southern Europe','Middle East', 'North America',
                                 'Oceanic','Asia', 'South America'), cex=0.9, text.font=4.5, bty="n",
       lwd = 3, col = c('blue', 'purple', 'red', 'black', 'tan1','green','brown','gray33'))



# Function to rotate coordinates by a specified angle (in degrees)
rotate_coords <- function(coords, angle_deg) {
  angle_rad <- angle_deg * pi / 180
  rotation_matrix <- matrix(c(cos(angle_rad), -sin(angle_rad), 
                              sin(angle_rad), cos(angle_rad)), 
                            nrow = 2)
  rotated <- coords %*% rotation_matrix
  return(rotated)
}

# Apply rotation to your coordinates (30 degrees in this example)
rotation_angle <- 50  # Change this to your desired rotation angle
rotated_mtx <- rotate_coords(plot_mtx, rotation_angle)

# Plot with rotated coordinates
par(cex = 0.9, mar = c(0.1, 0.1, 0.1, 0.1))
plot(rotated_mtx[,1], rotated_mtx[,2], col = color, pch = 19,
     xlab = "", ylab = "", bty = "n", axes = FALSE)

# Add text labels with adjusted positions
text(rotated_mtx[index, 1], rotated_mtx[index, 2], country[index], 
     pos = 1, cex = 0.6, offset = 0.4)
text(rotated_mtx[index_left,1], rotated_mtx[index_left,2], 
     country[index_left], pos = 2, cex = 0.6, offset = 0.4)
text(rotated_mtx[index_right,1], rotated_mtx[index_right,2], 
     country[index_right], pos = 4, cex = 0.6, offset = 0.4)
text(rotated_mtx[index_up,1], rotated_mtx[index_up,2], 
     country[index_up], pos = 3, cex = 0.6, offset = 0.4)

# Add legend (position might need adjustment after rotation)
legend("topleft", 
       legend = c('Nordic', 'West Europe', 'Central&Eastern&Southern Europe',
                  'Middle East', 'North America', 'Oceanic', 'Asia', 'South America'), 
       cex = 0.7, text.font = 4.5, bty = "n",
       pch = 19, col = c('blue', 'purple', 'red', 'black', 'tan1', 'green', 'brown', 'gray33'))





