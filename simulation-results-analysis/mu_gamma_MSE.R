library(ggplot2)

### p = 10 ####
p10_s1_n10000=c(mean(mu_MSE_N10000_r2_p10_J15[,1]),
         mean(mu_MSE_N10000_r2_p10_J15_rmsd_0.05[,1]),
         mean(mu_MSE_N10000_r2_p10_J15_rmsd_0.1[,1]),
         mean(mu_MSE_N10000_r2_p10_J15_rmsd_0.15[,1]),
         mean(mu_MSE_N10000_r2_p10_J15_bcmk[,1]))
round(p10_s1_n10000, 3)

p10_s1_n20000=c(mean(mu_MSE_N20000_r2_p10_J15[,1]),
                mean(mu_MSE_N20000_r2_p10_J15_rmsd_0.05[,1]),
                mean(mu_MSE_N20000_r2_p10_J15_rmsd_0.1[,1]),
                mean(mu_MSE_N20000_r2_p10_J15_rmsd_0.15[,1]),
                mean(mu_MSE_N20000_r2_p10_J15_bcmk[,1]))
round(p10_s1_n20000, 3)

### p = 20 ####

p20_s1_n10000=c(mean(mu_MSE_N10000_r2_p20_J30[,1]),
                mean(mu_MSE_N10000_r2_p20_J30_rmsd_0.05[,1]),
                mean(mu_MSE_N10000_r2_p20_J30_rmsd_0.1[,1]),
                mean(mu_MSE_N10000_r2_p20_J30_rmsd_0.15[,1]),
                mean(mu_MSE_N10000_r2_p20_J30_bcmk[,1]))
round(p20_s1_n10000, 3)

p20_s1_n20000=c(mean(mu_MSE_N20000_r2_p20_J30[,1]),
                mean(mu_MSE_N20000_r2_p20_J30_rmsd_0.05[,1]),
                mean(mu_MSE_N20000_r2_p20_J30_rmsd_0.1[,1]),
                mean(mu_MSE_N20000_r2_p20_J30_rmsd_0.15[,1]),
                mean(mu_MSE_N20000_r2_p20_J30_bcmk[,1]))
round(p20_s1_n20000, 3)

p10_s1_n10000=c(mean(sigma_MSE_N10000_r2_p10_J15[,1]),
                mean(sigma_MSE_N10000_r2_p10_J15_rmsd_0.05[,1]),
                mean(sigma_MSE_N10000_r2_p10_J15_rmsd_0.1[,1]),
                mean(sigma_MSE_N10000_r2_p10_J15_rmsd_0.15[,1]),
                mean(sigma_MSE_N10000_r2_p10_J15_bcmk[,1]))
round(p10_s1_n10000, 3)

p10_s1_n20000=c(mean(sigma_MSE_N20000_r2_p10_J15[,1]),
                mean(sigma_MSE_N20000_r2_p10_J15_rmsd_0.05[,1]),
                mean(sigma_MSE_N20000_r2_p10_J15_rmsd_0.1[,1]),
                mean(sigma_MSE_N20000_r2_p10_J15_rmsd_0.15[,1]),
                mean(sigma_MSE_N20000_r2_p10_J15_bcmk[,1]))
round(p10_s1_n20000, 3)

### p = 20 ####

p20_s1_n10000=c(mean(sigma_MSE_N10000_r2_p20_J30[,1]),
                mean(sigma_MSE_N10000_r2_p20_J30_rmsd_0.05[,1]),
                mean(sigma_MSE_N10000_r2_p20_J30_rmsd_0.1[,1]),
                mean(sigma_MSE_N10000_r2_p20_J30_rmsd_0.15[,1]),
                mean(sigma_MSE_N10000_r2_p20_J30_bcmk[,1]))
round(p20_s1_n10000, 3)

p20_s1_n20000=c(mean(sigma_MSE_N20000_r2_p20_J30[,1]),
                mean(sigma_MSE_N20000_r2_p20_J30_rmsd_0.05[,1]),
                mean(sigma_MSE_N20000_r2_p20_J30_rmsd_0.1[,1]),
                mean(sigma_MSE_N20000_r2_p20_J30_rmsd_0.15[,1]),
                mean(sigma_MSE_N20000_r2_p20_J30_bcmk[,1]))
round(p20_s1_n20000, 3)

p10_s1_n10000=c(mean(g_MSE_N10000_r2_p10_J15[,1]),
                mean(g_MSE_N10000_r2_p10_J15_rmsd_0.05[,1]),
                mean(g_MSE_N10000_r2_p10_J15_rmsd_0.1[,1]),
                mean(g_MSE_N10000_r2_p10_J15_rmsd_0.15[,1]))
round(p10_s1_n10000, 3)

p10_s1_n20000=c(mean(g_MSE_N20000_r2_p10_J15[,1]),
                mean(g_MSE_N20000_r2_p10_J15_rmsd_0.05[,1]),
                mean(g_MSE_N20000_r2_p10_J15_rmsd_0.1[,1]),
                mean(g_MSE_N20000_r2_p10_J15_rmsd_0.15[,1]))
#mean(g_MSE_N20000_r2_p10_J15_bcmk[,1]))
round(p10_s1_n20000, 3)

### p = 20 ####

p20_s1_n10000=c(mean(g_MSE_N10000_r2_p20_J30[,1]),
                mean(g_MSE_N10000_r2_p20_J30_rmsd_0.05[,1]),
                mean(g_MSE_N10000_r2_p20_J30_rmsd_0.1[,1]),
                mean(g_MSE_N10000_r2_p20_J30_rmsd_0.15[,1]))
round(p20_s1_n10000, 3)

p20_s1_n20000=c(mean(g_MSE_N20000_r2_p20_J30[,1]),
                mean(g_MSE_N20000_r2_p20_J30_rmsd_0.05[,1]),
                mean(g_MSE_N20000_r2_p20_J30_rmsd_0.1[,1]),
                mean(g_MSE_N20000_r2_p20_J30_rmsd_0.15[,1]))
round(p20_s1_n20000, 3)


