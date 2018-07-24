#JE model with Montagu template: getting the cohort column
#loading library:
library(xlsx)

#loading data:
input_data_path = 'D:/OUCRU/Hannah/JEV/JEV_model/data_JE_clean/Montagu_data/test2/'
naive_pop = read.xlsx2(paste(input_data_path,'naive_pop_2000_2044.xlsx', sep = ""),1)
naive_pop_less_15_2000_2044 = naive_pop[as.numeric(as.character((naive_pop$age_from))) < 15,]
load(paste(input_data_path,'cam_pop_less_15_2000_2044.RData', sep = ""))
load(paste(input_data_path,'rou_pop_less_15_2000_2044.RData', sep = ""))

Det_Mon_template = read.csv(paste(input_data_path,"deterministic_burden_template_JE-OUCRU-Clapham.csv", sep = ""))
Sto_Mon_template = read.csv(paste(input_data_path,"stochastic_burden_template_JE-OUCRU-Clapham.csv", sep = ""))

countr_ls_Det_temp = as.character(unique(Det_Mon_template$country))
countr_ls_Sto_temp = as.character(unique(Sto_Mon_template$country))
countr_ls_demo = as.character(unique(naive_pop_less_15_2000_2044$country))
 
countr_ls_demo == countr_ls_Det_temp  #Checking if the order of countries are wrong.
countr_ls_demo == countr_ls_Sto_temp
l_1_coun_temp = 15*45
Det_order_sim = match(countr_ls_demo,countr_ls_Det_temp)
Det_order = as.vector(sapply(Det_order_sim, FUN = f <- function(x){seq((x-1)*l_1_coun_temp + 1, x*l_1_coun_temp)}))
Sto_order_sim = match(countr_ls_demo,countr_ls_Sto_temp)
Sto_order = as.vector(sapply(Sto_order_sim, FUN = f <- function(x){seq((x-1)*l_1_coun_temp + 1, x*l_1_coun_temp)}))
###Fill-in information of "cohort_size" colummn
##Determistic:
#No vaccine scenario:
Det_No_vac_pop = as.numeric(as.vector(t(naive_pop_less_15_2000_2044[,-c(1,2,3)])))[Det_order]
Det_No_vac_template = Det_Mon_template
Det_No_vac_template$cohort_size = Det_No_vac_pop

#Campaign vaccine scenario:
Det_cam_vac_pop = unlist(lapply(cam_pop_less_15, FUN = f<-function(x)as.numeric(as.vector(t(x)))))[Det_order]
Det_cam_vac_template = Det_Mon_template
Det_cam_vac_template$cohort_size = Det_cam_vac_pop

#Routine vaccine scenario:
Det_rou_vac_pop = unlist(lapply(rou_pop_less_15, FUN = f<-function(x)as.numeric(as.vector(t(x)))))[Det_order]
Det_rou_vac_template = Det_Mon_template
Det_rou_vac_template$cohort_size = Det_rou_vac_pop

##Stochastic:
runs = length(unique(Sto_Mon_template$run_id))

#No vaccine scenario:
Sto_No_vac_template = Sto_Mon_template
Sto_No_vac_template$cohort_size = rep(Det_No_vac_pop, each = runs)

#Campaign vaccine scenario:
Sto_cam_vac_template = Sto_Mon_template
Sto_cam_vac_template$cohort_size = rep(Det_cam_vac_pop, each = runs)

#Routine vaccine scenario:
Sto_rou_vac_template = Sto_Mon_template
Sto_rou_vac_template$cohort_size = rep(Det_rou_vac_pop, each = runs)

#save template:
write.csv(Det_No_vac_template, file = paste(input_data_path, "Det_No_vac_template.csv"), row.names = F)
write.csv(Det_cam_vac_template, file = paste(input_data_path, "Det_Cam_vac_template.csv"), row.names = F)
write.csv(Det_rou_vac_template, file = paste(input_data_path, "Det_Rou_vac_template.csv"), row.names = F)

write.csv(Sto_No_vac_template, file = paste(input_data_path, "Sto_No_vac_template.csv"), row.names = F)
write.csv(Sto_cam_vac_template, file = paste(input_data_path, "Sto_Cam_vac_template.csv"), row.names = F)
write.csv(Sto_rou_vac_template, file = paste(input_data_path, "Sto_Rou_vac_template.csv"), row.names = F)

