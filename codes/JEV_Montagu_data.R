rm(list = ls())
#Montagu data - test version 3 - 16 countries => generate naive population, campaign population
#Loading library:
library(xlsx)
#Loading data:
folder_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_raw/Montagu-data" 
demo_data = read.csv(paste(folder_path, "/201710gavi-2_dds-201710_int_pop_both.csv", sep = ""))
cov_cam_gavi = read.csv(paste(folder_path, "/coverage_201710gavi-2_je-campaign-gavi.csv", sep = ""))
cov_rou_gavi = read.csv(paste(folder_path, "/coverage_201710gavi-2_je-routine-gavi.csv", sep = ""))

#List of country:
country_ls = as.character(unique(cov_cam_gavi$country_code))

#Selected year:
year_sel = 1950:2100
age_sel = 0:99

##stratify data by countries:

#Country demographic data:
col_sel_demo = c("country", "age_from", "age_to", paste("X", year_sel, sep = ""))
coun_demo_all = data.frame(matrix(vector(), nrow = 0, ncol = length(col_sel_demo), dimnames = list(c(),col_sel_demo)), stringsAsFactors = F)

for(i in 1:length(country_ls)){
  coun_name = country_ls[i]

  coun_demo_r = demo_data[demo_data$country_code == coun_name,]
  coun_demo = coun_demo_r[order(coun_demo_r$year),]
  
  max_table_length = length(age_sel)
  year_age_pop = c()
  for(year in year_sel){
    year_age_pop_r = coun_demo[coun_demo$year == year & coun_demo$age_from %in% age_sel, "value"]
    length(year_age_pop_r) = max_table_length
    year_age_pop = cbind(year_age_pop, year_age_pop_r)
  }
  coun_demo_final = data.frame(coun_demo[coun_demo$year == year & coun_demo$age_from %in% age_sel, c("country_code", "age_from", "age_to")],
                               year_age_pop)
  colnames(coun_demo_final) = col_sel_demo

  coun_demo_all = rbind(coun_demo_all, coun_demo_final)
}

#Country coverage data campaign:
coun_cam_all = cov_cam_gavi[cov_cam_gavi$year %in% year_sel,
                            c("country","activity_type", "year", "age_first", "age_last", "target", "coverage")]

#Country coverage data routine:
coun_rou_all = cov_rou_gavi[cov_rou_gavi$year %in% year_sel
                            ,c("country", "year", "age_first", "age_last", "target", "coverage")]

##different population by scenario:
#No vaccine scenario:
naive_pop = coun_demo_all
naive_pop_age_sel = naive_pop[naive_pop$age_from == age_sel,]


naive_pop_age_sel_sum = c()
for(ii in 1:length(country_ls)){
  cou_age_sel = colSums(naive_pop_age_sel[naive_pop_age_sel$country_code %in% country_ls[ii],4:length(naive_pop_age_sel)])
  naive_pop_age_sel_sum = rbind(naive_pop_age_sel_sum, cou_age_sel)
}
naive_pop_age_sel_sum = data.frame(country_ls, naive_pop_age_sel_sum)
  
#Routine scenario:
rou_pop_age_sel_sum = c() #sum data
rou_pop_age_sel = list() #No sum data:
for(rou in 1:length(country_ls)){
  #the coverage data:
  coun_rou_cov = coun_rou_all$coverage[coun_rou_all$country == country_ls[rou]]
  coun_rou_cov[is.na(coun_rou_cov)] = 0
  
  #population in selected age group of the country:
  cou_age_sel_age = naive_pop_age_sel[naive_pop_age_sel$country %in% country_ls[rou],4:length(naive_pop_age_sel)]
  
  #generate the routine data - vaccinate the first age group, let the vaccinated ppl age:
  cou_vac_age_sel = data.frame(matrix(0, nrow = nrow(cou_age_sel_age), ncol = ncol(cou_age_sel_age),
                                      dimnames = list(age_sel, colnames(cou_age_sel_age))))
  for(vac in 1:ncol(cou_age_sel_age)){
    vac_ppl = cou_age_sel_age[1, vac]*coun_rou_cov[vac]
    for(aging in 1:nrow(cou_age_sel_age)){
      cou_vac_age_sel[aging,aging + vac - 1] = vac_ppl
    }
  }
  cou_vac_age_sel = cou_vac_age_sel[,1:length(coun_rou_cov)]
  
  #Vaccinated population:
  vac_pop = cou_age_sel_age - cou_vac_age_sel
  vac_pop[vac_pop < 0] <- 0
  rou_pop_age_sel_sum = rbind(rou_pop_age_sel_sum, colSums(vac_pop))
  rou_pop_age_sel = c(rou_pop_age_sel, list(vac_pop))
}
rou_pop_age_sel_sum =  data.frame(country_ls, rou_pop_age_sel_sum)
names(rou_pop_age_sel) = country_ls

#Campaign + Routine scenario:
cam_pop_age_sel = list() # No sum data
cam_pop_age_sel_sum = c() #sum data
cam_data = coun_cam_all[coun_cam_all$activity_type == "campaign",]
cam_data$target[is.na(cam_data$target)] = 0
rou_data = coun_cam_all[coun_cam_all$activity_type == "routine",]
for(cam in 1:length(country_ls)){
  ##Vaccinate routine first:
  #the coverage data:
  coun_rou_cov_cam = rou_data$coverage[rou_data$country == country_ls[cam]]
  coun_rou_cov_cam[is.na(coun_rou_cov_cam)] = 0
  coun_cam_cov = cam_data[cam_data$country == country_ls[cam],c("coverage", "year")]
  coun_cam_cov = coun_cam_cov[!is.na(coun_cam_cov$coverage),]
  
  #population age group < 15 years old of the country:
  cou_age_sel_age = naive_pop_age_sel[naive_pop_age_sel$country %in% country_ls[cam],4:length(naive_pop_age_sel)]
  
  #generate the routine data - vaccinate the first age group, let the vaccinated ppl age:
  cou_vac_age_sel = data.frame(matrix(0, nrow = nrow(cou_age_sel_age), ncol = ncol(cou_age_sel_age),
                                      dimnames = list(age_sel, colnames(cou_age_sel_age))))
  for(vac_rou in 1:ncol(cou_age_sel_age)){
    vac_ppl = cou_age_sel_age[1, vac_rou]*coun_rou_cov_cam[vac_rou]
    for(aging in 1:nrow(cou_age_sel_age)){
      cou_vac_age_sel[aging,aging + vac_rou - 1] = vac_ppl
    }
  }
  cou_vac_age_sel = cou_vac_age_sel[,1:length(coun_rou_cov_cam)]
  
  #Vaccinated population - routine:
  vac_pop_rou = cou_age_sel_age - cou_vac_age_sel
  vac_pop_rou[vac_pop_rou < 0] <- 0
  
  ##Vaccinate campaign later:
  cou_age_sel_age_prop = sweep(cou_age_sel_age, 2, colSums(cou_age_sel_age), FUN = "/")
  index_year_cam_vac = match(paste("X",coun_cam_cov$year, sep = "") , colnames(cou_age_sel_age_prop))
  age_group_cam = 1:15
  vac_pop_cam = vac_pop_rou
  for(vac_cam in 1:length(index_year_cam_vac)){
    cam_vac_age_sel = data.frame(matrix(0, nrow = nrow(cou_age_sel_age), ncol = ncol(cou_age_sel_age),
                                        dimnames = list(age_sel, colnames(cou_age_sel_age))))
    cam_year_cov = coun_cam_cov[vac_cam,1]
    index_cam_vac = index_year_cam_vac[vac_cam]
    
    #Aging process:
    for(aging in 0:(nrow(cou_age_sel_age)-1)){
      if(aging == 0) age_group_cam_aging = 1:nrow(cou_age_sel_age)
      else age_group_cam_aging = head(age_group_cam, - aging)
      cam_vac_age_sel[age_group_cam_aging + aging,index_cam_vac + aging] = vac_pop_cam[age_group_cam_aging,index_cam_vac]*cam_year_cov
    }
    cam_vac_age_sel = cam_vac_age_sel[, 1:length(coun_rou_cov_cam)]
    vac_pop_cam = vac_pop_cam - cam_vac_age_sel
  }
  vac_pop_cam[vac_pop_cam < 0] = 0
  cam_pop_age_sel_sum = rbind(cam_pop_age_sel_sum, colSums(vac_pop_cam))
  cam_pop_age_sel = c(cam_pop_age_sel, list(vac_pop_cam))
}
cam_pop_age_sel_sum =  data.frame(country_ls, cam_pop_age_sel_sum)
names(cam_pop_age_sel) = country_ls

#Saving data:
saving_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_clean/Montagu_data/test3"
write.xlsx(naive_pop, paste(saving_path, "naive_pop_1950_2100.xlsx", sep = "/"), row.names= F)
write.xlsx(naive_pop_age_sel_sum, paste(saving_path, "naive_pop__sum_age_sel_2000_2044.xlsx", sep = "/"), row.names= F)
write.xlsx(rou_pop_age_sel_sum, paste(saving_path, "rou_pop__sum_age_sel_2000_2044.xlsx", sep = "/"), row.names= F)
write.xlsx(cam_pop_age_sel_sum, paste(saving_path, "cam_pop__sum_age_sel_2000_2044.xlsx", sep = "/"), row.names= F)
save(rou_pop_age_sel, file = paste(saving_path, "rou_pop_age_sel_2000_2044.RData", sep = "/"))
save(cam_pop_age_sel, file = paste(saving_path, "cam_pop_age_sel_2000_2044.RData", sep = "/"))
write.xlsx(coun_cam_all, paste(saving_path, "coun_cam_all_2000_2044.xlsx", sep = "/"), row.names= F)
write.xlsx(coun_rou_all, paste(saving_path, "coun_rou_all_2000_2044.xlsx", sep = "/"), row.names= F)
