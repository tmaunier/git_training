#Montagu data - test version 1 - PAK and IND => generate naive population, scenario population
#Montagu data - test version 2 - PAK and IND => generate naive population, scenario population
##Montagu data - test version 3 - 15 countries => generate naive population, scenario population of 100 age groups
##Montagu data - test version 3 - 16 countries => generate naive population, scenario population of 100 age groups, 
#campaign assumption: age group - target*coverage (age groups are getting from data of naive_pop_age_sel (maybe local of national))

#get coverage data
folder_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_raw/Montagu-data" 
cov_cam_gavi = read.csv(paste(folder_path, "/coverage_201710gavi-2_je-campaign-gavi.csv", sep = ""))
cov_rou_gavi = read.csv(paste(folder_path, "/coverage_201710gavi-2_je-routine-gavi.csv", sep = ""))


##Function to generate routine population: format: coun_name: "IND", "PAK", ...; naive_pop_age_sel: table by age group and year with age group selected population.
rou_scen_pop_gen <- function(coun_name, naive_pop_age_sel, year_sel){
  age_sel = 0:99
  naive_pop_age_sel_data = naive_pop_age_sel[,paste("X",year_sel, sep = "")]
  naive_pop_age_sel_outline = naive_pop_age_sel[,c("country", "age_from", "age_to")]
  #Country coverage data routine:
  coun_rou_all = cov_rou_gavi[cov_rou_gavi$year %in% year_sel
                              ,c("country_code", "year", "age_first", "age_last", "target", "coverage")]
  
  #Routine scenario:
    #the coverage data:
    if(any(coun_rou_all$country_code == coun_name)){
      
      coun_rou_cov_year = coun_rou_all[coun_rou_all$country_code == coun_name, c("coverage", "year")]
      coun_rou_cov_year[is.na(coun_rou_cov_year)] = 0
      coun_rou_cov = coun_rou_cov_year$coverage[match(year_sel, coun_rou_cov_year$year)]
      coun_rou_cov[is.na(coun_rou_cov)] = 0
      
      #generate the routine data - vaccinate the first age group, let the vaccinated ppl age:
      cou_vac_age_sel = data.frame(matrix(0, nrow = nrow(naive_pop_age_sel_data), ncol = ncol(naive_pop_age_sel_data),
                                          dimnames = list(age_sel, colnames(naive_pop_age_sel_data))))
      for(vac in 1:ncol(naive_pop_age_sel_data)){
        vac_ppl = naive_pop_age_sel_data[1, vac]*coun_rou_cov[vac]
        for(aging in 1:nrow(naive_pop_age_sel_data)){
          cou_vac_age_sel[aging,aging + vac - 1] = vac_ppl
        }
      }
      cou_vac_age_sel = cou_vac_age_sel[,1:length(coun_rou_cov)]
      
      #Vaccinated population:
      vac_pop = naive_pop_age_sel_data - cou_vac_age_sel
      vac_pop[vac_pop < 0] <- 0
  
    return(cbind(naive_pop_age_sel_outline, vac_pop))
    } else { 
      print("Coverage data not found!") 
      return(naive_pop_age_sel) }
}

#Function to generate Cam + routine population: format: coun_name: "IDN", "PAK", ...; naive_pop_age_sel: table by age group and year with age group selected population.
cam_scen_pop_gen <-function(coun_name, naive_pop_age_sel, year_sel){
  age_sel = 0:99
  naive_pop_age_sel_data = naive_pop_age_sel[,paste("X",year_sel, sep = "")]
  naive_pop_age_sel_outline = naive_pop_age_sel[,c("country", "age_from", "age_to")]
  #Country coverage data routine:
  coun_rou_all = cov_rou_gavi[cov_rou_gavi$year %in% year_sel
                              ,c("country_code", "year", "age_first", "age_last", "target", "coverage")]

  #Country coverage data campaign:
  coun_cam_all = cov_cam_gavi[cov_cam_gavi$year %in% year_sel,
                              c("country_code","activity_type", "year", "age_first", "age_last", "target", "coverage")]
  
  cam_data = coun_cam_all[coun_cam_all$activity_type == "campaign",]
  cam_data$target[is.na(cam_data$target)] = 0
  rou_data = coun_cam_all[coun_cam_all$activity_type == "routine",]
  
  if(any(c(rou_data$country_code == coun_name,cam_data$country_code == coun_name))){

      ##Vaccinate routine first:
      #the coverage data:
      coun_rou_cov_cam_year = coun_rou_all[coun_rou_all$country_code == coun_name, c("coverage", "year")]
      coun_rou_cov_cam_year[is.na(coun_rou_cov_cam_year)] = 0
      coun_rou_cov_cam = coun_rou_cov_cam_year$coverage[match(year_sel, coun_rou_cov_cam_year$year)]
      coun_rou_cov_cam[is.na(coun_rou_cov_cam)] = 0

      coun_cam_cov_year = cam_data[cam_data$country_code == coun_name,c("target","coverage", "year")]
      coun_cam_cov_year = coun_cam_cov_year[!is.na(coun_cam_cov_year$coverage),]
      coun_cam_cov = coun_cam_cov_year[match(year_sel, coun_cam_cov_year$year),]
      coun_cam_cov$year = year_sel
      coun_cam_cov[is.na(coun_cam_cov)] = 0
      
      #generate the routine data - vaccinate the first age group, let the vaccinated ppl age:
      cou_vac_age_sel = data.frame(matrix(0, nrow = nrow(naive_pop_age_sel_data), ncol = ncol(naive_pop_age_sel_data),
                                          dimnames = list(age_sel, colnames(naive_pop_age_sel_data))))
      for(vac_rou in 1:ncol(naive_pop_age_sel_data)){
        vac_ppl = naive_pop_age_sel_data[1, vac_rou]*coun_rou_cov_cam[vac_rou]
        for(aging in 1:nrow(naive_pop_age_sel_data)){
          cou_vac_age_sel[aging,aging + vac_rou - 1] = vac_ppl
        }
      }
      cou_vac_age_sel = cou_vac_age_sel[,1:length(coun_rou_cov_cam)]
      
      #Vaccinated population - routine:
      vac_pop_rou = naive_pop_age_sel_data - cou_vac_age_sel
      vac_pop_rou[vac_pop_rou < 0] <- 0
      
      ##Vaccinate campaign later:
      age_group_cam = 1:15
      naive_pop_age_group_cam_sel_prop = sweep(naive_pop_age_sel_data[age_group_cam,], 2, colSums(naive_pop_age_sel_data[age_group_cam,]), FUN = "/")
      index_year_cam_vac = match(paste("X",coun_cam_cov$year, sep = "") , colnames(naive_pop_age_group_cam_sel_prop))
    
      target_number = as.numeric(as.vector(coun_cam_cov$target))*unlist(coun_cam_cov$coverage)
      vac_pop_cam = vac_pop_rou
      for(vac_cam in 1:length(index_year_cam_vac)){
        cam_vac_age_sel = data.frame(matrix(0, nrow = nrow(naive_pop_age_sel_data), ncol = ncol(naive_pop_age_sel_data),
                                            dimnames = list(age_sel, colnames(naive_pop_age_sel_data))))
        index_cam_vac = index_year_cam_vac[vac_cam]
        
        #Aging process:
        for(aging in 0:(nrow(naive_pop_age_sel_data)-1)){
          cam_vac_age_sel[age_group_cam + aging,index_cam_vac + aging] = naive_pop_age_group_cam_sel_prop[,index_cam_vac]*target_number[vac_cam]
        }
        cam_vac_age_sel = cam_vac_age_sel[1:nrow(naive_pop_age_sel_data), 1:length(coun_rou_cov_cam)]
        vac_pop_cam = vac_pop_cam - cam_vac_age_sel
      }
      vac_pop_cam[vac_pop_cam < 0] = 0
      
      return(cbind(naive_pop_age_sel_outline,vac_pop_cam))
  } else {
    print("Coverage data not found!") 
    return(naive_pop_age_sel)
    }
}