#Montagu data - test version 1 - PAK and IND => generate naive population, scenario population
#Montagu data - test version 2 - PAK and IND => generate naive population, scenario population

#Function to generate routine population: format: coun_name: "IND", "PAK", ...; naive_pop_less_15: table by age group and year only take <15 years old population.
rou_scen_pop_gen <- function(coun_name, naive_pop_less_15){
  folder_path = "D:/OUCRU/Hannah/JEV/JEV_model/data_JE_raw/Montagu-data" 
  cov_rou_gavi = read.csv(paste(folder_path, "/coverage_201708test-2_je-routine-gavi.csv", sep = ""))
  year_sel = 2000:2044
  #Country coverage data routine:
  coun_rou_all = cov_rou_gavi[cov_rou_gavi$year %in% year_sel
                              ,c("country", "year", "age_first", "age_last", "target", "coverage")]
  
  #Routine scenario:

    #the coverage data:
    
    if(any(coun_rou_all$country == coun_name)){
  
      coun_rou_cov = coun_rou_all$coverage[coun_rou_all$country == coun_name]
      coun_rou_cov[is.na(coun_rou_cov)] = 0
      
      #generate the routine data - vaccinate the first age group, let the vaccinated ppl age:
      cou_vac_less_15 = data.frame(matrix(0, nrow = nrow(naive_pop_less_15), ncol = ncol(naive_pop_less_15),
                                          dimnames = list(0:14, colnames(naive_pop_less_15))))
      for(vac in 1:ncol(naive_pop_less_15)){
        vac_ppl = naive_pop_less_15[1, vac]*coun_rou_cov[vac]
        for(aging in 1:nrow(naive_pop_less_15)){
          cou_vac_less_15[aging,aging + vac - 1] = vac_ppl
        }
      }
      cou_vac_less_15 = cou_vac_less_15[,1:length(coun_rou_cov)]
      
      #Vaccinated population:
      vac_pop = naive_pop_less_15 - cou_vac_less_15
      vac_pop[vac_pop < 0] <- 0
  
    return(vac_pop)
    } else { 
      print("Coverage data not found!") 
      return(naive_pop_less_15) }
}

#Function to generate Cam + routine population: format: coun_name: "IDN", "PAK", ...; naive_pop_less_15: table by age group and year only take <15 years old population.
cam_scen_pop_gen <-function(coun_name, naive_pop_less_15){
  folder_path = "D:/OUCRU/Hannah/JEV/JEV_model/data_JE_raw/Montagu-data" 
  cov_cam_gavi = read.csv(paste(folder_path, "/coverage_201708test-2_je-campaign-gavi.csv", sep = ""))
  cov_rou_gavi = read.csv(paste(folder_path, "/coverage_201708test-2_je-routine-gavi.csv", sep = ""))
  year_sel = 2000:2044
  
  #Country coverage data routine:
  coun_rou_all = cov_rou_gavi[cov_rou_gavi$year %in% year_sel
                              ,c("country", "year", "age_first", "age_last", "target", "coverage")]

  #Country coverage data campaign:
  coun_cam_all = cov_cam_gavi[cov_cam_gavi$year %in% year_sel,
                              c("country","activity_type", "year", "age_first", "age_last", "target", "coverage")]
  
  cam_data = coun_cam_all[coun_cam_all$activity_type == "campaign",]
  cam_data$target[is.na(cam_data$target)] = 0
  rou_data = coun_cam_all[coun_cam_all$activity_type == "routine",]
  
  if(any(c(rou_data$country == coun_name,cam_data$country == coun_name))){
    
      ##Vaccinate routine first:
      #the coverage data:
      coun_rou_cov_cam = rou_data$coverage[rou_data$country == coun_name]
      coun_rou_cov_cam[is.na(coun_rou_cov_cam)] = 0
      coun_cam_cov = cam_data[cam_data$country == coun_name,c("coverage", "year")]
      coun_cam_cov = coun_cam_cov[!is.na(coun_cam_cov$coverage),]
      
      #generate the routine data - vaccinate the first age group, let the vaccinated ppl age:
      cou_vac_less_15 = data.frame(matrix(0, nrow = nrow(naive_pop_less_15), ncol = ncol(naive_pop_less_15),
                                          dimnames = list(0:14, colnames(naive_pop_less_15))))
      for(vac_rou in 1:ncol(naive_pop_less_15)){
        vac_ppl = naive_pop_less_15[1, vac_rou]*coun_rou_cov_cam[vac_rou]
        for(aging in 1:nrow(naive_pop_less_15)){
          cou_vac_less_15[aging,aging + vac_rou - 1] = vac_ppl
        }
      }
      cou_vac_less_15 = cou_vac_less_15[,1:length(coun_rou_cov_cam)]
      
      #Vaccinated population - routine:
      vac_pop_rou = naive_pop_less_15 - cou_vac_less_15
      vac_pop_rou[vac_pop_rou < 0] <- 0
      
      ##Vaccinate campaign later:
      naive_pop_less_15_prop = sweep(naive_pop_less_15, 2, colSums(naive_pop_less_15), FUN = "/")
      index_year_cam_vac = match(paste("X",coun_cam_cov$year, sep = "") , colnames(naive_pop_less_15_prop))
      age_group_cam = 1:15
      vac_pop_cam = vac_pop_rou
      for(vac_cam in 1:length(index_year_cam_vac)){
        cam_vac_less_15 = data.frame(matrix(0, nrow = nrow(naive_pop_less_15), ncol = ncol(naive_pop_less_15),
                                            dimnames = list(0:14, colnames(naive_pop_less_15))))
        cam_year_cov = coun_cam_cov[vac_cam,1]
        index_cam_vac = index_year_cam_vac[vac_cam]
        
        #Aging process:
        for(aging in 0:(nrow(naive_pop_less_15)-1)){
          if(aging == 0) age_group_cam_aging = 1:nrow(naive_pop_less_15)
          else age_group_cam_aging = head(age_group_cam, - aging)
          cam_vac_less_15[age_group_cam_aging + aging,index_cam_vac + aging] = vac_pop_cam[age_group_cam_aging,index_cam_vac]*cam_year_cov
        }
        cam_vac_less_15 = cam_vac_less_15[, 1:length(coun_rou_cov_cam)]
        vac_pop_cam = vac_pop_cam - cam_vac_less_15
      }
      vac_pop_cam[vac_pop_cam < 0] = 0
      
      return(vac_pop_cam)
  } else {
    print("Coverage data not found!") 
    return(naive_pop_less_15)
    }
}