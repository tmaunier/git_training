#Calculate the susceptible population based on vaccination coverage:
#loading library:
library(xlsx)
pop_data_path = "D:/OUCRU/Hannah/JEV/JEV_model/data_JE_raw/"
output_data_path = 'D:/OUCRU/Hannah/JEV/JEV_model/data_JE_clean/pop_data/'
data_name = "vaccine_data.xlsx"
Rou_Best = read.xlsx(paste(pop_data_path,data_name, sep = ""), sheetName = "Routine_Coverage_BestEstimate")
Rou_NG = read.xlsx(paste(pop_data_path,data_name, sep = ""), sheetName = "Routine_Coverage_NoGavi")
Cam_Best = read.xlsx(paste(pop_data_path,data_name, sep = ""), sheetName = "Campaign_Coverage_BestEstimate")
Cam_NG = read.xlsx(paste(pop_data_path,data_name, sep = ""), sheetName = "Campaign_Coverage_NoGavi")
demo_population = read.xlsx2(paste(output_data_path,'IDB_0_100_2000_2045_24.xlsx', sep = ""),1)
demo_population[,4:49] <- mapply(demo_population[,4:49], FUN = f <- function(x) as.numeric(as.character(x)))

#List of country need:
ISO_country_list = as.character(unique(demo_population$ISO))
num_age_group = length(unique(demo_population$Age_group))
#Routine_Coverage_BestEstimate
Rou_Best_func <- function(){
Rou_Best_demo_2000_2045 = demo_population[0,]
for(ii in 1:length(ISO_country_list)){
  country = ISO_country_list[ii]
  Rou_Best_country = Rou_Best[Rou_Best$ISO == country,]
  demo_country = demo_population[demo_population$ISO == country,]
  if (nrow(Rou_Best_country) != 0){
    index_2000_year_vac = grep("X2000", colnames(Rou_Best_country))
    index_2000_year_pop = grep("X2000", colnames(demo_country))
    vaccinated_infant = as.numeric(Rou_Best_country[index_2000_year_vac:(index_2000_year_vac + 45)]*demo_country[1,index_2000_year_pop:(index_2000_year_pop + 45)])
    vaccinated_pop_2000_2045 = matrix(ncol = length(vaccinated_infant), nrow = nrow(demo_country))
    seq_sel = min(nrow(demo_country), ncol(demo_country) - 3)
    for (i in (1:seq_sel-1)){
      vaccinated_pop_2000_2045[i,] = c(rep(NA,i),vaccinated_infant[1:(seq_sel - i)])
    }
    vaccinated_pop_2000_2045[is.na(vaccinated_pop_2000_2045)] <-0
    sus_pop_2000_2045 = demo_country[,4:49] - vaccinated_pop_2000_2045
    sus_pop_2000_2045[sus_pop_2000_2045 < 0] <-0
    sus_pop_2000_2045 = cbind(demo_country[,1:3],sus_pop_2000_2045)
  } else {
    sus_pop_2000_2045 = demo_country
  }
  Rou_Best_demo_2000_2045 = rbind(Rou_Best_demo_2000_2045, sus_pop_2000_2045)
}
return(Rou_Best_demo_2000_2045)
}
Rou_Best_demo_2000_2045 <- Rou_Best_func()
write.xlsx(Rou_Best_demo_2000_2045, paste(output_data_path, "Rou_Best_demo_2000_2045_",num_age_group,".xlsx", sep = ""),
           row.names=F)

#Routine_Coverage_NoGavi
Rou_NG_func <- function(){
  Rou_NG_demo_2000_2045 = demo_population[0,]
  for(ii in 1:length(ISO_country_list)){
    country = ISO_country_list[ii]
    Rou_NG_country = Rou_NG[Rou_NG$ISO == country,]
    demo_country = demo_population[demo_population$ISO == country,]
    if (nrow(Rou_NG_country) != 0){
      index_2000_year_vac = grep("X2000", colnames(Rou_NG_country))
      index_2000_year_pop = grep("X2000", colnames(demo_country))
      vaccinated_infant = as.numeric(Rou_NG_country[index_2000_year_vac:(index_2000_year_vac + 45)]*demo_country[1,index_2000_year_pop:(index_2000_year_pop + 45)])
      vaccinated_pop_2000_2045 = matrix(ncol = length(vaccinated_infant), nrow = nrow(demo_country))
      seq_sel = min(nrow(demo_country), ncol(demo_country) - 3)
      for (i in (1:seq_sel-1)){
        vaccinated_pop_2000_2045[i,] = c(rep(NA,i),vaccinated_infant[1:(seq_sel - i)])
      }
      vaccinated_pop_2000_2045[is.na(vaccinated_pop_2000_2045)] <-0
      sus_pop_2000_2045 = demo_country[,4:49] - vaccinated_pop_2000_2045
      sus_pop_2000_2045[sus_pop_2000_2045 < 0] <-0
      sus_pop_2000_2045 = cbind(demo_country[,1:3],sus_pop_2000_2045)
    } else {
      sus_pop_2000_2045 = demo_country
    }
    Rou_NG_demo_2000_2045 = rbind(Rou_NG_demo_2000_2045, sus_pop_2000_2045)
  }
  return(Rou_NG_demo_2000_2045)
}
Rou_NG_demo_2000_2045 <- Rou_NG_func()
write.xlsx(Rou_NG_demo_2000_2045, paste(output_data_path, "Rou_NG_demo_2000_2045_",num_age_group,".xlsx", sep = ""),
           row.names=F)

#Campaign_Coverage_BestEstimate
Cam_Best_func <- function(){
  Cam_Best_demo_2000_2045 = demo_population[0,]
  for(ii in 1:length(ISO_country_list)){
    country = ISO_country_list[ii]
    Cam_Best_country = Cam_Best[Cam_Best$ISO == country,]
    Cam_Best_country = Cam_Best_country[order(Cam_Best[Cam_Best$ISO == country,"Year"]),]
    demo_country = demo_population[demo_population$ISO == country,]
    n_Cam_vac = nrow(Cam_Best_country)
    if (n_Cam_vac != 0){
      index_2000_year_pop = grep("X2000", colnames(demo_country))
      Age_target = 1:nrow(demo_country)
      Cam_Year = Cam_Best_country$Year
      percent_vac = Cam_Best_country$Percent.target.vaccinated
      percent_vac[is.na(percent_vac)] <- 0
      demo_temp = demo_country[,4:49]
      for(year in 1:length(Cam_Year)){
        vaccinated_Cam_Year_total = matrix(0, nrow = length(Age_target), ncol = ncol(demo_country))
        vac_year_index = grep(paste("X",Cam_Year[year], sep = ""), colnames(demo_temp))
        vaccinated_Cam_Year =  demo_temp[,vac_year_index]*percent_vac[year]
        for (i in Age_target){
          tryCatch({ #ignore error message
          vaccinated_Cam_Year_total[,(vac_year_index - 1 + i)] = c(rep(0,length.out = i - 1),vaccinated_Cam_Year[1:(length(Age_target) - i + 1)])
          }, error=function(e){})
        }
        demo_temp = demo_temp - vaccinated_Cam_Year_total
      }
      demo_temp[demo_temp < 0] <- 0
      sus_pop_2000_2045 = cbind(demo_country[,1:3],demo_temp)
    } else {
      sus_pop_2000_2045 = demo_country
    }
    Cam_Best_demo_2000_2045 = rbind(Cam_Best_demo_2000_2045, sus_pop_2000_2045)
  }
  return(Cam_Best_demo_2000_2045)
}
Cam_Best_demo_2000_2045 <- Cam_Best_func()
write.xlsx(Cam_Best_demo_2000_2045, paste(output_data_path, "Cam_Best_demo_2000_2045_",num_age_group,".xlsx", sep = ""),
           row.names=F)

#Campaign_Coverage_NoGavi
Cam_NG_func <- function(){
  Cam_NG_demo_2000_2045 = demo_population[0,]
  for(ii in 1:length(ISO_country_list)){
    country = ISO_country_list[ii]
    Cam_NG_country = Cam_NG[Cam_NG$ISO == country,]
    Cam_NG_country = Cam_NG_country[order(Cam_NG[Cam_NG$ISO == country,"Year"]),]
    demo_country = demo_population[demo_population$ISO == country,]
    n_Cam_vac = nrow(Cam_NG_country)
    if (n_Cam_vac != 0){
      index_2000_year_pop = grep("X2000", colnames(demo_country))
      Age_target = 1:nrow(demo_country)
      Cam_Year = Cam_NG_country$Year
      percent_vac = Cam_NG_country$Percent.target.vaccinated
      percent_vac[is.na(percent_vac)] <- 0
      demo_temp = demo_country[,4:49]
      for(year in 1:length(Cam_Year)){
        vaccinated_Cam_Year_total = matrix(0, nrow = length(Age_target), ncol = ncol(demo_country))
        vac_year_index = grep(paste("X",Cam_Year[year], sep = ""), colnames(demo_temp))
        vaccinated_Cam_Year =  demo_temp[,vac_year_index]*percent_vac[year]
        for (i in Age_target){
          tryCatch({
          vaccinated_Cam_Year_total[,(vac_year_index - 1 + i)] = c(rep(0,length.out = i - 1),vaccinated_Cam_Year[1:(length(Age_target) - i + 1)])
          }, error=function(e){})
        }
        demo_temp = demo_temp - vaccinated_Cam_Year_total
      }
      demo_temp[demo_temp < 0] <- 0
      sus_pop_2000_2045 = cbind(demo_country[,1:3],demo_temp)
    } else {
      sus_pop_2000_2045 = demo_country
    }
    Cam_NG_demo_2000_2045 = rbind(Cam_NG_demo_2000_2045, sus_pop_2000_2045)
  }
  return(Cam_NG_demo_2000_2045)
}
Cam_NG_demo_2000_2045 <- Cam_NG_func()
write.xlsx(Cam_NG_demo_2000_2045, paste(output_data_path, "Cam_NG_demo_2000_2045_",num_age_group,".xlsx", sep = ""),
           row.names=F)
