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
library(mize)
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

country=c('Australia', 'Austria', 'Belgium', 'Canada', 'Chile', 'Colombia', 'Costa Rica', 'Czech Republic', 'Denmark',        
          'Estonia',  'Finland', 'France', 'Germany', 'Greece', 'Hungary', 'Iceland', 'Ireland', 'Israel', 'Italy', 'Japan', 'Korea',          
          'Latvia', 'Lithuania', 'Mexico', 'Netherlands', 'New Zealand','Norway', 'Poland', 'Portugal', 'Slovak Republic', 'Slovenia', 'Spain',          
          'Sweden', 'Switzerland', 'Turkey', 'United Kingdom', 'United States')

length(item_math)#169


g_hat <- as.matrix(read.csv("gamma_math_2022_d0.csv", header = F))
a_hat <- read.csv("alpha_math_2022_d0.csv", header = F)[,1]
b_hat <- read.csv("beta_math_2022_d0.csv", header = F)[,1]
mu_hat <- read.csv("mu_math_2022_d0.csv", header = F)[,1]
sigma_hat <- read.csv("sigma_math_2022_d0.csv", header = F)[,1]

#g_hat = as.matrix(g_hat)
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


trans_d_matrix = matrix(transform_d, ncol=3, byrow=F)
output_d_matrix = cbind(1:57, trans_d_matrix[,1], 58:114,trans_d_matrix[,2], 115:169, trans_d_matrix[,3] )
xtable(as.matrix(output_d_matrix), type = "latex", file = "compare_latent_skill.tex", digits=2)

trans_a_matrix = matrix(a_hat, ncol=3, byrow=F)
output_a_matrix = cbind(1:57, trans_a_matrix[,1], 58:114,trans_a_matrix[,2], 115:169, trans_a_matrix[,3] )
xtable(as.matrix(output_a_matrix), type = "latex", file = "compare_latent_skill.tex", digits=2)

transform_p=c(t(transform_p)) #J*p dim vector

gam_val=read.csv('gamma_transformed_math_2022_d0.csv', header=F)
gam_val=gam_val[,1]

#### histogram ####

df <- data.frame(gam_val = gam_val)

ggplot(df, aes(x = gam_val)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "white") +
  labs(
    #title = "Histogram of γ Values",
    x = expression(hat(gamma)[jk]),
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )



### Ranking 0.05


mu_hat_0.05=read.csv( 'mu_math_RMSD_2022_0.05.csv', header=F)[,1]
#mu_hat_2PL=c(mu_hat_2PL)
theta_ord_0.05=order(mu_hat_0.05, decreasing = T)
country_ord_0.05=country[theta_ord_0.05]
country_ord_0.05
sigma_hat_0.05=read.csv( 'sigma_math_RMSD_2022_0.05.csv', header=F)[,1]

round(mu_hat_0.05[theta_ord_0.05],3)
round(sigma_hat_0.05[theta_ord_0.05], 3)

# [1] "Japan"           "Korea"           "Estonia"         "Switzerland"     "Czech Republic"  "Poland"         
# [7] "Belgium"         "Austria"         "Netherlands"     "Australia"       "Ireland"         "Sweden"         
# [13] "Latvia"          "Germany"         "Canada"          "Spain"           "Hungary"         "Denmark"        
# [19] "New Zealand"     "United Kingdom"  "Finland"         "Italy"           "Portugal"        "Norway"         
# [25] "Slovenia"        "Slovak Republic" "Lithuania"       "France"          "Iceland"         "Israel"         
# [31] "United States"   "Turkey"          "Chile"           "Greece"          "Costa Rica"      "Mexico"         
# [37] "Colombia"     

#To evaluate the Kendall's rank correlations
cor(mu_hat, mu_hat_0.05, method = "kendall")
#0.963964


#Compare difference in country ranking for 2PL model with and without DIF
#2PL without DIF
mu_hat_0.1=read.csv( 'mu_math_RMSD_2022_0.1.csv', header=F)[,1]
sigma_hat_0.1=read.csv( 'sigma_math_RMSD_2022_0.1.csv', header=F)[,1]

#mu_hat_2PL=c(mu_hat_2PL)
theta_ord_0.1=order(mu_hat_0.1, decreasing = T)
country_ord_0.1=country[theta_ord_0.1]
country_ord_0.1


round(mu_hat_0.1[theta_ord_0.1],3)
round(sigma_hat_0.1[theta_ord_0.1], 3)


# [1] "Japan"           "Korea"           "Estonia"         "Switzerland"     "Czech Republic"  "Poland"         
# [7] "Belgium"         "Austria"         "Netherlands"     "Australia"       "Ireland"         "Sweden"         
# [13] "Latvia"          "Germany"         "Canada"          "New Zealand"     "United Kingdom"  "Spain"          
# [19] "Hungary"         "Denmark"         "Finland"         "Italy"           "Norway"          "Portugal"       
# [25] "Slovak Republic" "Slovenia"        "Lithuania"       "France"          "Iceland"         "Israel"         
# [31] "United States"   "Turkey"          "Greece"          "Chile"           "Mexico"          "Colombia"       
# [37] "Costa Rica"          


#Compare difference in country ranking for 2PL model with and without DIF
#2PL without DIF
mu_hat_0.15=read.csv( 'mu_math_RMSD_2022_0.15.csv', header=F)[,1]
sigma_hat_0.15=read.csv( 'sigma_math_RMSD_2022_0.15.csv', header=F)[,1]

#mu_hat_2PL=c(mu_hat_2PL)
theta_ord_0.15=order(mu_hat_0.15, decreasing = T)
country_ord_0.15=country[theta_ord_0.15]
country_ord_0.15


round(mu_hat_0.15[theta_ord_0.15],3)
round(sigma_hat_0.15[theta_ord_0.15], 3)

# [1] "Japan"           "Korea"           "Estonia"         "Switzerland"     "Czech Republic"  "Poland"         
# [7] "Belgium"         "Austria"         "Netherlands"     "Australia"       "Ireland"         "Sweden"         
# [13] "Germany"         "Canada"          "Latvia"          "United Kingdom"  "New Zealand"     "Hungary"        
# [19] "Spain"           "Italy"           "Denmark"         "Finland"         "Portugal"        "Norway"         
# [25] "Slovenia"        "France"          "Slovak Republic" "Lithuania"       "Iceland"         "Israel"         
# [31] "United States"   "Turkey"          "Greece"          "Chile"           "Mexico"          "Colombia"       
# [37] "Costa Rica"    
## baseline 


mu_hat_baseline=read.csv( 'mu_math_2022_baseline.csv',header=F)[,1]
sigma_hat_baseline=read.csv( 'sigma_math_2022_baseline.csv', header=F)[,1]

theta_ord_baseline=order(mu_hat_baseline, decreasing = T)
country_ord_baseline=country[theta_ord_baseline]
country_ord_baseline

round(mu_hat_baseline[theta_ord_baseline],3)
round(sigma_hat_baseline[theta_ord_baseline],3)


# [1] "Japan"           "Korea"           "Estonia"         "Switzerland"     "Czech Republic"  "Poland"         
# [7] "Belgium"         "Austria"         "Netherlands"     "Australia"       "Ireland"         "Sweden"         
# [13] "Germany"         "Latvia"          "Canada"          "United Kingdom"  "Spain"           "New Zealand"    
# [19] "Hungary"         "Italy"           "Denmark"         "Finland"         "Portugal"        "Norway"         
# [25] "Slovenia"        "Slovak Republic" "France"          "Lithuania"       "Israel"          "Iceland"        
# [31] "United States"   "Turkey"          "Greece"          "Chile"           "Mexico"          "Colombia"       
# [37] "Costa Rica"     


xtable(cbind(country_ord_0.1, round(mu_hat_0.1[theta_ord_0.1],3), round(sigma_hat_0.1[theta_ord_0.1], 3), 1:37,
             country_ord_0.15, round(mu_hat_0.15[theta_ord_0.15],3), round(sigma_hat_0.15[theta_ord_0.15], 3), 1:37,
             country_ord_baseline, round(mu_hat_baseline[theta_ord_baseline],3), round(sigma_hat_baseline[theta_ord_baseline], 3)), 
       type = "latex", file = "compare_latent_skill.tex", digits=3)



#[1] "Japan"           "Korea"           "Estonia"         "Switzerland"     "Czech Republic"  "Poland"         
# [7] "Belgium"         "Austria"         "Netherlands"     "Australia"       "Ireland"         "Sweden"         
# [13] "Germany"         "Latvia"          "Canada"          "United Kingdom"  "Spain"           "New Zealand"    
# [19] "Hungary"         "Italy"           "Denmark"         "Finland"         "Portugal"        "Norway"         
# [25] "Slovenia"        "Slovak Republic" "France"          "Lithuania"       "Israel"          "Iceland"        
# [31] "United States"   "Turkey"          "Greece"          "Chile"           "Mexico"          "Colombia"       
# [37] "Costa Rica"     

#Proposed method
mu_hat_full=transform_mu
#mu_hat_full=c(0, mu_hat_full)
theta_ord_full=order(mu_hat_full, decreasing = T)
country_ord_full=country[theta_ord_full]
country_ord_full

# "Japan"           "Korea"           "Estonia"         "Switzerland"     "Czech Republic"  "Belgium"        
# "Poland"          "Netherlands"     "Austria"         "Ireland"         "Australia"       "Sweden"         
# "Germany"         "Latvia"          "Canada"          "Hungary"         "United Kingdom"  "Spain"          
# "New Zealand"     "Denmark"         "Finland"         "Italy"           "Portugal"        "Norway"         
# "Slovenia"        "Slovak Republic" "Lithuania"       "France"          "Iceland"         "Israel"         
# "United States"   "Turkey"          "Greece"          "Chile"           "Mexico"          "Colombia"       
# "Costa Rica"     



###### MDS ###########

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

set.seed(40)
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
rotation_angle <- 60  # Change this to your desired rotation angle
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
legend("topright", 
       legend = c('Nordic', 'West Europe', 'Central&Eastern&Southern Europe',
                  'Middle East', 'North America', 'Oceanic', 'Asia', 'South America'), 
       cex = 0.7, text.font = 4.5, bty = "n",
       pch = 19, col = c('blue', 'purple', 'red', 'black', 'tan1', 'green', 'brown', 'gray33'))


