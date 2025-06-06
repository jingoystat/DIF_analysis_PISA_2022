library(haven)

##To read in the sampling weights
df_weights=read_sas('~/STU_QQQ_SAS/CY08MSP_STU_QQQ.SAS7BDAT')
sampling_weights=df_weights$W_FSTUWT

##To read in the cognitive data file
df=read_sas('~/STU_COG_SAS/CY08MSP_STU_COG.SAS7BDAT')

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

temp_science <- read.csv("science_names.csv", sep=',', header=F)
temp_science = temp_science[,1]
item_science<- c()
for(item in temp_science){
  print(unique(df[, item]))
  if(setequal(unique(df[, item]), c(NA, 0, 1))){
    item_science = c(item_science, item)
  }
}

data_science=df[which(df$OECD==1), c(item_science, 'CNT', 'sampling_weights')] 
dim(data_science)#295157    106
unique(data_science$CNT) #37 countries




#Remove items with all NAs
ind_na_col=c()
for (i in 1:106){
  if (sum(is.na(data_science[[i]]))==295157){
    ind_na_col=c(ind_na_col, i)
  }
}
ind_na_col
#There is no item with all NA's

#Remove examinees with all NAs
#First convert into a dataframe
dat=as.data.frame(data_science)
ind_na_row=which(rowSums(is.na(dat[,1:106])) == 104) #na
ind_na_row
dat=dat[-ind_na_row, ]
dim(dat) #124015    106



#To ensure reliable analysis, we only use items with more than 10000 observations and we only consider binary scored items

ind_1000_col=c()
obs = c()
for (i in 1:104){
  if (sum(!is.na(dat[[i]]))>1000){
    ind_1000_col=c(ind_1000_col, i)
    obs = c(obs, sum(!is.na(dat[[i]])))
  }
}
ind_1000_col
obs


#We also removed all the items that have at least one of the 37 OECD countries that does not have 
#any response to; those are items with indices 44-51

dat=dat[ , names(dat) %in% c(item_science, 'CNT', 'sampling_weights')]
write.csv(dat$sampling_weights, file = 'PISA_2022_science_samplingweights_sas.csv')


response=dat[,1:104]
for (i in 1:104){
  response[,i]=as.numeric(response[,i])
}
write.csv(response, file = 'PISA_2022_science_sas.csv')



#Create x feature matrix
country=c('Australia', 'Austria', 'Belgium', 'Canada', 'Chile', 'Colombia', 'Costa Rica', 'Czech Republic', 'Denmark',        
          'Estonia',  'Finland', 'France', 'Germany', 'Greece', 'Hungary', 'Iceland', 'Ireland', 'Israel', 'Italy', 'Japan', 'Korea',          
          'Latvia', 'Lithuania', 'Mexico', 'Netherlands', 'New Zealand','Norway', 'Poland', 'Portugal', 'Slovak Republic', 'Slovenia', 'Spain',          
          'Sweden', 'Switzerland', 'Turkey', 'United Kingdom', 'United States')

country=dat$CNT
oecd_country=c("AUS", "AUT", "BEL", "CAN", "CHL", "COL", "CRI", "CZE", "DNK", "EST", "FIN", "FRA",
               "DEU", "GRC", "HUN", "ISL", "IRL", "ISR", "ITA", "JPN", "KOR", "LVA", "LTU", "MEX",
               "NLD", "NZL", "NOR", "POL", "PRT", "SVK", "SVN", "ESP", "SWE", "CHE", "TUR", "GBR", "USA")



ind_con=c()
for (i in 1:length(country)){
  ind=which(oecd_country==country[i])
  ind_con=c(ind_con, ind)
}

dat_country=matrix(0, nrow = length(country), ncol = 37)
for (i in 1: length(country)){
  dat_country[i, ind_con[i]]=1
}


#dat_country=dat_country[,-1]#Remove the baseline

write.csv(dat_country, file = 'countries_2022_science_sas.csv')
