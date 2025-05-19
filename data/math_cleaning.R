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


data_math=df[which(df$OECD==1), c(item_math, 'CNT', 'sampling_weights')] 
dim(data_math)#295157    171
unique(data_math$CNT) #37 countries




#Remove items with all NAs
ind_na_col=c()
for (i in 1:169){
  if (sum(is.na(data_math[[i]]))==295157){
    ind_na_col=c(ind_na_col, i)
  }
}
ind_na_col
#There is no item with all NA's

#Remove examinees with all NAs
#First convert into a dataframe
dat=as.data.frame(data_math)
ind_na_row=which(rowSums(is.na(dat[,1:169])) == 167) #18019
dat=dat[-ind_na_row, ]
dim(dat) #277138    171



#To ensure reliable analysis, we only use items with more than 10000 observations and we only consider binary scored items

ind_1000_col=c()
obs = c()
for (i in 1:169){
  if (sum(!is.na(data_math[[i]]))>1000){
    ind_1000_col=c(ind_1000_col, i)
    obs = c(obs, sum(!is.na(data_math[[i]])))
  }
}
ind_1000_col
obs


#We also removed all the items that have at least one of the 37 OECD countries that does not have 
#any response to; those are items with indices 44-51

dat=dat[ , names(dat) %in% c(item_math, 'CNT', 'sampling_weights')]
write.csv(dat$sampling_weights, file = 'PISA_2022_math_samplingweights_sas.csv')


response=dat[,1:169]
for (i in 1:169){
  response[,i]=as.numeric(response[,i])
}
write.csv(response, file = 'PISA_2022_math_sas.csv')



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

write.csv(dat_country, file = 'countries_2022_math_sas.csv')







































