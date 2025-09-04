

library(TwoSampleMR)
library(MRPRESSO) 
library(ggplot2)

#####Read harmonized data of diet and aortic dissection#####
dat_diet_AD <-read.csv(file="diet_AD_xth.csv")

res_diet_AD<- mr(dat_diet_AD)
res_diet_AD_OR <-generate_odds_ratios(res_diet_AD) 



#Positive results were selected for sensitivity analysis and 
#visualization of results, respectively 
dat_fish_AD <-read.csv(file="fish_AD_xth.csv")
res_fish_AD<- mr(dat_fish_AD)


dat_ofish_AD <-read.csv(file="ofish_AD_xth.csv")
res_ofish_AD<- mr(dat_ofish_AD)

####Fish consumption  and   Aorta Dissection####

#Get effect values for each SNP individually 
res_single_fAD <- mr_singlesnp(dat_fish_AD)
res_single_fAD

#leave-one-out analysis 
res_loo_fAD <- mr_leaveoneout(dat_fish_AD)
res_loo_fAD 


#Scatter plot 
p1 <- mr_scatter_plot(res_fish_AD, dat_fish_AD)
p1[[1]]

ggsave(p1[[1]], file="p1.tif", width=7, height=7)


#Forest plot
p2 <- mr_forest_plot(res_single_fAD)
p2[[1]]

ggsave(p2[[1]], file="p2.tif", width=7, height=20)


#leave-one-out analysis plot
p3 <- mr_leaveoneout_plot(res_loo_fAD)
p3[[1]]

ggsave(p3[[1]], file="p3.tif", width=7, height=7)


#Funnel plot 
p4 <- mr_funnel_plot(res_single_fAD)
p4[[1]]

ggsave(p4[[1]], file="p4.tif", width=7, height=7)


#heterogeneity test
mr_heterogeneity(dat_fish_AD) 


#MR PRESSO
run_mr_presso(dat_fish_AD,NbDistribution = 1000)

#horizontal pleiotropy test
mr_pleiotropy_test(dat_fish_AD)  


####oily Fish consumption  and   Aorta Dissection####

#Get effect values for each SNP individually 
res_single_ofAD <- mr_singlesnp(dat_ofish_AD)
res_single_ofAD

#leave-one-out analysis 
res_loo_ofAD <- mr_leaveoneout(dat_ofish_AD)
res_loo_ofAD 


#Scatter plot 
p1 <- mr_scatter_plot(res_ofish_AD, dat_ofish_AD)
p1[[1]]

ggsave(p1[[1]], file="po1.tif", width=7, height=7)


#Forest plot
p2 <- mr_forest_plot(res_single_ofAD)
p2[[1]]

ggsave(p2[[1]], file="p02.tif", width=7, height=20)


#leave-one-out analysis plot
p3 <- mr_leaveoneout_plot(res_loo_ofAD)
p3[[1]]

ggsave(p3[[1]], file="p03.tif", width=7, height=7)


#Funnel plot 
p4 <- mr_funnel_plot(res_single_ofAD)
p4[[1]]

ggsave(p4[[1]], file="p04.tif", width=7, height=7)


#heterogeneity test
mr_heterogeneity(dat_ofish_AD) 


#MR_PRESSO
run_mr_presso(dat_ofish_AD,NbDistribution = 1000)

#horizontal pleiotropy test
mr_pleiotropy_test(dat_ofish_AD) 



#####Read harmonized data of diet and Aortic Aneurysm#####
dat_diet_AA <-read.csv(file="diet_AA_xth.csv")

res_diet_AA<- mr(dat_diet_AA)
res_diet_AA_OR <-generate_odds_ratios(res_diet_AA) 


#Positive results were selected for sensitivity analysis and 
#visualization of results, respectively 

dat_Dfruit_AA <-read.csv(file="Dfruit_AA_xth.csv")
res_Dfruit_AA<- mr(dat_Dfruit_AA)


dat_ofish_AA <-read.csv(file="ofish_AA_xth.csv")
res_ofish_AA<- mr(dat_ofish_AA)


dat_sala_AA <-read.csv(file="sala_AA_xth.csv")
res_sala_AA<- mr(dat_sala_AA)

####Dried fruit consumption and Aortic Aneurysm ####

#Get effect values for each SNP individually 
res_single_DfruitAA <- mr_singlesnp(dat_Dfruit_AA)
res_single_DfruitAA 

#leave-one-out analysis 
res_loo_DfruitAA <- mr_leaveoneout(dat_Dfruit_AA)
res_loo_DfruitAA  


#Scatter plot  
p1 <- mr_scatter_plot(res_Dfruit_AA, dat_Dfruit_AA)
p1[[1]]

ggsave(p1[[1]], file="p11.tif", width=7, height=7)


#Forest plot
p2 <- mr_forest_plot(res_single_DfruitAA)
p2[[1]]

ggsave(p2[[1]], file="p12.tif", width=7, height=20)


#leave-one-out analysis plot
p3 <- mr_leaveoneout_plot(res_loo_DfruitAA)
p3[[1]]

ggsave(p3[[1]], file="p13.tif", width=7, height=7)


#Funnel plot
p4 <- mr_funnel_plot(res_single_DfruitAA)
p4[[1]]

ggsave(p4[[1]], file="p14.tif", width=7, height=7)


#heterogeneity test
mr_heterogeneity(dat_Dfruit_AA) 


#MR PRESSO
run_mr_presso(dat_Dfruit_AA,NbDistribution = 1000)

#horizontal pleiotropy test
mr_pleiotropy_test(dat_Dfruit_AA)  


#### Oily fish consumption and Aortic Aneurysm ####

#Get effect values for each SNP individually 
res_single_ofishAA <- mr_singlesnp(dat_ofish_AA)
res_single_ofishAA 

#leave-one-out analysis 
res_loo_ofishAA <- mr_leaveoneout(dat_ofish_AA)
res_loo_ofishAA  


#Scatter plot  
p1 <- mr_scatter_plot(res_ofish_AA, dat_ofish_AA)
p1[[1]]

ggsave(p1[[1]], file="p011.tif", width=7, height=7)


#Forest plot
p2 <- mr_forest_plot(res_single_ofishAA)
p2[[1]]

ggsave(p2[[1]], file="p012.tif", width=7, height=20)


#leave-one-out analysis plot
p3 <- mr_leaveoneout_plot(res_loo_ofishAA)
p3[[1]]

ggsave(p3[[1]], file="p013.tif", width=7, height=7)


#Funnel plot
p4 <- mr_funnel_plot(res_single_ofishAA)
p4[[1]]

ggsave(p4[[1]], file="p014.tif", width=7, height=7)


#heterogeneity test
mr_heterogeneity(dat_ofish_AA) 


#MR_PRESSO
run_mr_presso(dat_ofish_AA,NbDistribution = 1000)

#horizontal pleiotropy test
mr_pleiotropy_test(dat_ofish_AA)  



#### Salad consumption and Aortic Aneurysm ####

#Get effect values for each SNP individually 
res_single_salaAA <- mr_singlesnp(dat_sala_AA)
res_single_salaAA 

#leave-one-out analysis 
res_loo_salaAA <- mr_leaveoneout(dat_sala_AA)
res_loo_salaAA  


#Scatter plot  
p1 <- mr_scatter_plot(res_sala_AA, dat_sala_AA)
p1[[1]]

ggsave(p1[[1]], file="p021.tif", width=7, height=7)


#Forest plot
p2 <- mr_forest_plot(res_single_salaAA)
p2[[1]]

ggsave(p2[[1]], file="p022.tif", width=7, height=20)


#leave-one-out analysis plot 
p3 <- mr_leaveoneout_plot(res_loo_salaAA)
p3[[1]]

ggsave(p3[[1]], file="p023.tif", width=7, height=7)


#Funnel plot
p4 <- mr_funnel_plot(res_single_salaAA)
p4[[1]]

ggsave(p4[[1]], file="p024.tif", width=7, height=7)


#heterogeneity test
mr_heterogeneity(dat_sala_AA) 


#MR_PRESSO
run_mr_presso(dat_sala_AA,NbDistribution = 1000)

#horizontal pleiotropy test
mr_pleiotropy_test(dat_sala_AA)  



