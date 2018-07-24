#Code to test the template:
rm(list = ls())
library(ggplot2)
library(dplyr)
library(magrittr)
library(data.table)

#loading template:
Det_template_save_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/result_model/first_model/Det_template/"
#Det_template_save_path = "C:/Users/DELL-PC/Dropbox/JEV_GAVI_Gates_model/result_model/first_model/Det_template/"
Det_No_vac = fread(file = paste(Det_template_save_path, "test3_Det_No_vac.csv", sep = ""))
Det_Cam_vac = fread(file = paste(Det_template_save_path, "test3_Det_Cam_vac.csv", sep = ""))
Det_Rou_vac = fread(file = paste(Det_template_save_path, "test3_Det_Rou_vac.csv", sep = ""))

Sto_template_save_path = "D:/OUCRU/Hannah/JEV/JEV_model/result_model/first_model/Sto_template/stochastic_burden_template_JE-OUCRU-Clapham"
Sto_No_vac = fread(file = paste(Sto_template_save_path, "-No-vaccination.csv", sep = ""))
Sto_Cam_vac = fread(file = paste(Sto_template_save_path, "-Campaign.csv", sep = ""))
Sto_Rou_vac = fread(file = paste(Sto_template_save_path, "-Routine.csv", sep = ""))

choosen_col = c("cohort_size","cases",  "deaths", "dalys")
#choosen_col = c("cohort_size","cases",  "deaths", "dalys", "diff_vs_cam", "diff_vs_rou") #choose column to plot
#function to plot Determinstic model
Det_plotting <-function(Det_temp){
  list_plot = list()
  for(i in 1:length(choosen_col)){
    temp_plot = ggplot(Det_temp, aes_string(x = "year", y = choosen_col[i], group = "age"))+
      geom_line(aes(color = age))+
      facet_wrap(~country_name, scales = "free") +
      ggtitle(choosen_col[i])
    list_plot = c(list_plot, list(temp_plot))
  }
  names(list_plot) = choosen_col
  return(list_plot)
}
#function to plot Stochastic model
Sto_plotting <-function(Sto_temp, Det_temp){
  list_plot = list()
  for(i in 1:length(choosen_col)){
    temp_plot = ggplot(Sto_temp, aes_string(x = "year", y = choosen_col[i], group = "age"))+
      geom_line(data = Det_temp, mapping = aes_string(x = "year", y = choosen_col[i], 
                                                      group = "age", alpha = 0.1), color = "grey")+
      geom_line(aes(color = age))+
      facet_grid(country ~ run_id, scales = "free") +
      ggtitle(choosen_col[i])+
      theme(panel.background = element_blank(), axis.text.x = element_text(angle = -90))
    list_plot = c(list_plot, list(temp_plot))
  }
  return(list_plot)
}

#####PLOTTING:
Det_No_vac_plot = Det_plotting(Det_No_vac)
Det_Cam_vac_plot = Det_plotting(Det_Cam_vac)
Det_Rou_vac_plot = Det_plotting(Det_Rou_vac)

Sto_No_vac_plot = Sto_plotting(Sto_No_vac, Det_No_vac)
Sto_Cam_vac_plot = Sto_plotting(Sto_Cam_vac, Det_Cam_vac)
Sto_Rou_vac_plot = Sto_plotting(Sto_Rou_vac, Det_Rou_vac)

plot_save_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/result_model/first_model/Testing_template/"

save_list_plot_func <- function(plot_list){
  for(i in 1:length(plot_list)){
    ggsave(filename = paste("test3_",deparse(substitute(plot_list)),"_", names(plot_list)[i],".png", sep = ""), plot = plot_list[[i]], 
           path = plot_save_path, units = "mm", width = 1600/5, height = 840/5)
  }
}

save_list_plot_func(Det_No_vac_plot)
save_list_plot_func(Det_Cam_vac_plot)
save_list_plot_func(Det_Rou_vac_plot)

save_list_plot_func(Sto_No_vac_plot)
save_list_plot_func(Sto_Cam_vac_plot)
save_list_plot_func(Sto_Rou_vac_plot)

##################
#checking template:
#total cases in every age groups:
choosen_col = c("cohort_size","deaths", "cases","dalys")
sum_by_ages_Det_No_vac <- aggregate(Det_No_vac[,choosen_col_ages], by = list(year = Det_No_vac$year, country_name = Det_No_vac$country_name), sum)
sum_by_ages_Det_No_vac$age = 1
sum_by_ages_Det_No_vac_plot <- Det_plotting(sum_by_ages_Det_No_vac)
save_list_plot_func(sum_by_ages_Det_No_vac_plot)

sum_by_ages_Det_Cam_vac <- aggregate(Det_Cam_vac[,choosen_col_ages], by = list(year = Det_Cam_vac$year, country_name = Det_Cam_vac$country_name), sum)
sum_by_ages_Det_Cam_vac$age = 1
sum_by_ages_Det_Cam_vac_plot <- Det_plotting(sum_by_ages_Det_Cam_vac)
save_list_plot_func(sum_by_ages_Det_Cam_vac_plot)

sum_by_ages_Det_Rou_vac <- aggregate(Det_Rou_vac[,choosen_col_ages], by = list(year = Det_Rou_vac$year, country_name = Det_Rou_vac$country_name), sum)
sum_by_ages_Det_Rou_vac$age = 1
sum_by_ages_Det_Rou_vac_plot <- Det_plotting(sum_by_ages_Det_Rou_vac)
save_list_plot_func(sum_by_ages_Det_Rou_vac_plot)
#total cases in every age groups globally:
sum_global_Det_No_vac <- aggregate(Det_No_vac[,choosen_col_ages], by = list(year = Det_No_vac$year), sum)
sum_global_Det_No_vac$age = 1
sum_global_Det_No_vac$country_name = "World"
sum_global_Det_No_vac_plot <- Det_plotting(sum_global_Det_No_vac)
save_list_plot_func(sum_global_Det_No_vac_plot)

sum_global_Det_Cam_vac <- aggregate(Det_Cam_vac[,choosen_col_ages], by = list(year = Det_Cam_vac$year), sum)
sum_global_Det_Cam_vac$age = 1
sum_global_Det_Cam_vac$country_name = "World"
sum_global_Det_Cam_vac_plot <- Det_plotting(sum_global_Det_Cam_vac)
save_list_plot_func(sum_global_Det_Cam_vac_plot)

sum_global_Det_Rou_vac <- aggregate(Det_Rou_vac[,choosen_col_ages], by = list(year = Det_Rou_vac$year), sum)
sum_global_Det_Rou_vac$age = 1
sum_global_Det_Rou_vac$country_name = "World"
sum_global_Det_Rou_vac_plot <- Det_plotting(sum_global_Det_Rou_vac)
save_list_plot_func(sum_global_Det_Rou_vac_plot)

#cases averted vs vaccinated number:
Det_No_vac$cases_aversed_cam = Det_No_vac[,"cases"] - Det_Cam_vac[,"cases"]
Det_No_vac$cases_aversed_rou = Det_No_vac[,"cases"] - Det_Rou_vac[,"cases"]

vaccinated_number = fread("C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/result_model/first_model/Pop_scen.csv")
Det_No_vac$vac_cam_number = vaccinated_number$No_vac - vaccinated_number$Cam
Det_No_vac$vac_rou_number = vaccinated_number$No_vac - vaccinated_number$Rou

choosen_col = c("cohort_size","cases", "cases_aversed_cam", "cases_aversed_rou", "vac_cam_number", "vac_rou_number") #choose column to plot
sum_by_ages_cases_aversed_vac_num = aggregate(Det_No_vac[,choosen_col], by = list(year = Det_No_vac$year, country_name = Det_No_vac$country_name), sum)
sum_by_ages_cases_aversed_vac_num$age = 1
#calculate ratio between cases aversed and 
sum_by_ages_cases_aversed_vac_num$ratio_case_vs_cam = sum_by_ages_cases_aversed_vac_num$cases_aversed_cam/sum_by_ages_cases_aversed_vac_num$vac_cam_number
sum_by_ages_cases_aversed_vac_num$ratio_case_vs_cam[sum_by_ages_cases_aversed_vac_num$ratio_case_vs_cam == -Inf | sum_by_ages_cases_aversed_vac_num$ratio_case_vs_cam == Inf] = NaN
sum_by_ages_cases_aversed_vac_num$ratio_case_vs_rou = sum_by_ages_cases_aversed_vac_num$cases_aversed_rou/sum_by_ages_cases_aversed_vac_num$vac_rou_number
sum_by_ages_cases_aversed_vac_num$ratio_case_vs_rou[sum_by_ages_cases_aversed_vac_num$ratio_case_vs_rou == -Inf | sum_by_ages_cases_aversed_vac_num$ratio_case_vs_rou == Inf] = NaN
#save plot
choosen_col = c("cohort_size","cases", "cases_aversed_cam", "cases_aversed_rou", "vac_cam_number", "vac_rou_number", "ratio_case_vs_cam", "ratio_case_vs_rou")
sum_by_ages_cases_aversed_vac_num_plot <- Det_plotting(sum_by_ages_cases_aversed_vac_num)
save_list_plot_func(sum_by_ages_cases_aversed_vac_num_plot)

# large_diff_cam = Det_No_vac[cases_averted_cam < -1,]
# large_diff_cam$diff_vs_cam = cases_averted_cam[cases_averted_cam < -1]
# large_diff_rou = Det_No_vac[cases_averted_rou < -1,]
# large_diff_rou$diff_vs_rou = cases_averted_rou[cases_averted_rou < -1]

#####New Sto template with 200 run ids group by country
sel_run_id = 6:10
country_list = Det_No_vac$country %>% as.character %>% unique
Sto_template_save_path = "D:/OUCRU/Hannah/JEV/test3_16coun_200runid_every_countries/"
template_file_name = "stochastic_burden_template_JE-OUCRU-Clapham_"

Sto_No_vac = c()
Sto_Cam_vac = c()
Sto_Rou_vac = c()

for(run in sel_run_id){
  print(run)
  Sto_No_vac_temp <- fread(paste(Sto_template_save_path, template_file_name, "je-routine-no-vaccination_",run,".csv", sep = ""))
  Sto_Cam_vac_temp <- fread(paste(Sto_template_save_path, template_file_name, "je-campaign-gavi_",run,".csv", sep = ""))
  Sto_Rou_vac_temp <- fread(paste(Sto_template_save_path, template_file_name, "je-routine-gavi_",run,".csv", sep = ""))
  
  Sto_No_vac = c(Sto_No_vac, list(Sto_No_vac_temp))
  Sto_Cam_vac = c(Sto_Cam_vac, list(Sto_Cam_vac_temp))
  Sto_Rou_vac = c(Sto_Rou_vac, list(Sto_Rou_vac_temp))
}

#combine data for all countries:
Sto_No_vac = Reduce(rbind, Sto_No_vac) 
Sto_Cam_vac = Reduce(rbind, Sto_Cam_vac) 
Sto_Rou_vac = Reduce(rbind, Sto_Rou_vac) 

Sto_No_vac_plot = Sto_plotting(Sto_No_vac, Det_No_vac)
Sto_Cam_vac_plot = Sto_plotting(Sto_Cam_vac, Det_Cam_vac)
Sto_Rou_vac_plot = Sto_plotting(Sto_Rou_vac, Det_Rou_vac)

save_list_plot_func(Sto_No_vac_plot)
save_list_plot_func(Sto_Cam_vac_plot)
save_list_plot_func(Sto_Rou_vac_plot)

#####New Sto template with 200 run ids group by run id: compare between Sto and Det model
sel_run_id = 1:200
country_list = Det_No_vac$country %>% as.character %>% unique
Sto_template_save_path = "D:/OUCRU/Hannah/JEV/test3_16coun_200runid_every_countries/"
template_file_name = "stochastic_burden_template_JE-OUCRU-Clapham_"

Sto_No_vac = data.frame(matrix(0, ncol = 4, nrow = nrow(Det_No_vac), dimnames = list(NULL, c("cohort_size","deaths", "cases", "dalys"))))
Sto_Cam_vac = Sto_No_vac
Sto_Rou_vac = Sto_No_vac

for(run in sel_run_id){
  print(run)
  Sto_No_vac_temp <- fread(paste(Sto_template_save_path, template_file_name, "je-routine-no-vaccination_",run,".csv", sep = ""))
  Sto_Cam_vac_temp <- fread(paste(Sto_template_save_path, template_file_name, "je-campaign-gavi_",run,".csv", sep = ""))
  Sto_Rou_vac_temp <- fread(paste(Sto_template_save_path, template_file_name, "je-routine-gavi_",run,".csv", sep = ""))
  
  Sto_No_vac = Sto_No_vac + Sto_No_vac_temp[,colnames(Sto_No_vac)]
  Sto_Cam_vac = Sto_Cam_vac + Sto_Cam_vac_temp[,colnames(Sto_No_vac)]
  Sto_Rou_vac = Sto_Rou_vac + Sto_Rou_vac_temp[,colnames(Sto_No_vac)]
}

Sto_No_vac = Sto_No_vac/length(sel_run_id)
Sto_Cam_vac = Sto_Cam_vac/length(sel_run_id)
Sto_Rou_vac = Sto_Rou_vac/length(sel_run_id)

choosen_col = c("cohort_size","deaths", "cases","dalys")
average_Sto_No_vac_plot = cbind(Sto_No_vac_temp[,c("year", "age", "country", "country_name", "run_id")], Sto_No_vac) %>% Sto_plotting(Det_No_vac)
average_Sto_Cam_vac_plot = cbind(Sto_Cam_vac_temp[,c("year", "age", "country", "country_name", "run_id")], Sto_Cam_vac) %>% Sto_plotting(Det_Cam_vac)
average_Sto_Rou_vac_plot = cbind(Sto_Rou_vac_temp[,c("year", "age", "country", "country_name", "run_id")], Sto_Rou_vac) %>% Sto_plotting(Det_Rou_vac)

choosen_col = c("cohort_size","deaths", "cases","dalys")
Det_vs_average_Sto_No <- cbind(select(Det_No_vac, -one_of(choosen_col)), Sto_No_vac-Det_No_vac[,choosen_col]) %>% magrittr::extract(choosen_col) %>% aggregate(by = list(year = Det_No_vac$year, country_name = Det_No_vac$country_name), sum)
Det_vs_average_Sto_No$age = 1
Det_vs_average_Sto_No_plot <- Det_plotting(Det_vs_average_Sto_No)
Det_vs_average_Sto_Rou_plot <- cbind(select(Det_No_vac, -one_of(choosen_col)), Sto_Rou_vac/Det_Rou_vac[,choosen_col]) %>% Det_plotting

tail(Sto_No_vac[which(Det_No_vac$country == "TLS"),])
tail(Det_No_vac[which(Det_No_vac$country == "TLS"),])

#####New Sto template with 200 run ids group by run id: compare between Sto and Det model: plot all run id and the det model with total cases in every ages
sel_run_id = 1:200
Sto_template_save_path = "D:/OUCRU/Hannah/JEV/test3_16coun_200runid_every_countries_rounded/"
template_file_name = "stochastic_burden_template_JE-OUCRU-Clapham_"

l_age_sum_Sto_No_vac = c()
l_age_sum_Sto_Cam_vac = c()
l_age_sum_Sto_Rou_vac = c()

choosen_col = c("run_id","cohort_size","deaths", "cases","dalys")
for(run in sel_run_id){
  print(run)
  Sto_No_vac_temp <- fread(paste(Sto_template_save_path, template_file_name, "je-routine-no-vaccination_",run,".csv", sep = ""))
  Sto_Cam_vac_temp <- fread(paste(Sto_template_save_path, template_file_name, "je-campaign-gavi_",run,".csv", sep = ""))
  Sto_Rou_vac_temp <- fread(paste(Sto_template_save_path, template_file_name, "je-routine-gavi_",run,".csv", sep = ""))
  
  sum_by_ages_Sto_No_vac_temp <- aggregate(Sto_No_vac_temp[,choosen_col, with = F], by = list(year = Sto_No_vac_temp$year, country_name = Sto_No_vac_temp$country_name, run_id = Sto_No_vac_temp$run_id), sum)
  sum_by_ages_Sto_Cam_vac_temp <- aggregate(Sto_Cam_vac_temp[,choosen_col, with = F], by = list(year = Sto_Cam_vac_temp$year, country_name = Sto_Cam_vac_temp$country_name, run_id = Sto_No_vac_temp$run_id), sum)
  sum_by_ages_Sto_Rou_vac_temp <- aggregate(Sto_Rou_vac_temp[,choosen_col, with = F], by = list(year = Sto_Rou_vac_temp$year, country_name = Sto_Rou_vac_temp$country_name, run_id = Sto_No_vac_temp$run_id), sum)

  l_age_sum_Sto_No_vac = c(l_age_sum_Sto_No_vac, list(sum_by_ages_Sto_No_vac_temp))
  l_age_sum_Sto_Cam_vac = c(l_age_sum_Sto_Cam_vac, list(sum_by_ages_Sto_Cam_vac_temp))
  l_age_sum_Sto_Rou_vac = c(l_age_sum_Sto_Rou_vac, list(sum_by_ages_Sto_Rou_vac_temp))
}

age_sum_Sto_No_vac = Reduce(rbind, l_age_sum_Sto_No_vac)
age_sum_Sto_Cam_vac = Reduce(rbind, l_age_sum_Sto_Cam_vac)
age_sum_Sto_Rou_vac = Reduce(rbind, l_age_sum_Sto_Rou_vac)

choosen_col = c("cohort_size","deaths", "cases","dalys")
sum_by_ages_Det_No_vac <- aggregate(Det_No_vac[,choosen_col, with = F], by = list(year = Det_No_vac$year, country_name = Det_No_vac$country_name), sum)
sum_by_ages_Det_Cam_vac <- aggregate(Det_Cam_vac[,choosen_col, with = F], by = list(year = Det_Cam_vac$year, country_name = Det_Cam_vac$country_name), sum)
sum_by_ages_Det_Rou_vac <- aggregate(Det_Rou_vac[,choosen_col, with = F], by = list(year = Det_Rou_vac$year, country_name = Det_Rou_vac$country_name), sum)

#function to plot Stochastic model
Sto_plotting_age_sum <-function(Sto_temp, Det_temp){
  list_plot = list()
  for(i in 1:length(choosen_col)){
    temp_plot = ggplot(Sto_temp, aes_string(x = "year", y = choosen_col[i]))+
      geom_line(aes(group = run_id), color = "blue", alpha = 0.1)+
      geom_line(data = Det_temp, mapping = aes_string(x = "year", y = choosen_col[i]), color = "red", size = 1)+
      facet_wrap(~country_name, scales = "free") +
      ggtitle(choosen_col[i])+
      theme(panel.background = element_blank(), axis.text.x = element_text(angle = -90))
    list_plot = c(list_plot, list(temp_plot))
  }
  return(list_plot)
}

age_sum_Sto_plotting_No_vac <- Sto_plotting_age_sum(age_sum_Sto_No_vac, sum_by_ages_Det_No_vac)
names(age_sum_Sto_plotting_No_vac) = choosen_col
age_sum_Sto_plotting_Cam_vac <- Sto_plotting_age_sum(age_sum_Sto_Cam_vac, sum_by_ages_Det_Cam_vac)
names(age_sum_Sto_plotting_Cam_vac) = choosen_col
age_sum_Sto_plotting_Rou_vac <- Sto_plotting_age_sum(age_sum_Sto_Rou_vac, sum_by_ages_Det_Rou_vac)
names(age_sum_Sto_plotting_Rou_vac) = choosen_col

plot_save_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/result_model/first_model/Testing_template/"
save_list_plot_func(age_sum_Sto_plotting_No_vac)
save_list_plot_func(age_sum_Sto_plotting_Cam_vac)
save_list_plot_func(age_sum_Sto_plotting_Rou_vac)

#####New Sto template with 200 run ids group by run id: compare between Sto and Det model: plot all run id and the det model with cases in a specific age
sel_run_id = 1:200
age_sel = c(0, 1, 49, 50, 98, 99)
Sto_template_save_path = "D:/OUCRU/Hannah/JEV/test3_16coun_200runid_every_countries_rounded/"
template_file_name = "stochastic_burden_template_JE-OUCRU-Clapham_"

l_age_i_Sto_No_vac = c()
l_age_i_Sto_Cam_vac = c()
l_age_i_Sto_Rou_vac = c()

choosen_col = c("run_id","cohort_size","deaths", "cases","dalys")
for(run in sel_run_id){
  print(run)
  Sto_No_vac_temp <- fread(paste(Sto_template_save_path, template_file_name, "je-routine-no-vaccination_",run,".csv", sep = ""), showProgress = F)
  Sto_Cam_vac_temp <- fread(paste(Sto_template_save_path, template_file_name, "je-campaign-gavi_",run,".csv", sep = ""), showProgress = F)
  Sto_Rou_vac_temp <- fread(paste(Sto_template_save_path, template_file_name, "je-routine-gavi_",run,".csv", sep = ""), showProgress = F)
  
  l_age_i_Sto_No_vac = c(l_age_i_Sto_No_vac, list(filter(Sto_No_vac_temp, age %in% age_sel)))
  l_age_i_Sto_Cam_vac = c(l_age_i_Sto_Cam_vac, list(filter(Sto_Cam_vac_temp, age %in% age_sel)))
  l_age_i_Sto_Rou_vac = c(l_age_i_Sto_Rou_vac, list(filter(Sto_Rou_vac_temp, age %in% age_sel)))
}

age_i_Sto_No_vac = Reduce(rbind, l_age_i_Sto_No_vac)
age_i_Sto_Cam_vac = Reduce(rbind, l_age_i_Sto_Cam_vac)
age_i_Sto_Rou_vac = Reduce(rbind, l_age_i_Sto_Rou_vac)

i_age_Det_No_vac <- filter(Det_No_vac, age %in% age_sel)
i_age_Det_Cam_vac <- filter(Det_Cam_vac, age %in% age_sel)
i_age_Det_Rou_vac <- filter(Det_Rou_vac, age %in% age_sel)

choosen_col = c("cohort_size","deaths", "cases","dalys")
Sto_plotting_i_age <-function(Sto_temp, Det_temp, i_age){
  Sto_temp = filter(Sto_temp, age == i_age)
  Det_temp = filter(Det_temp, age == i_age)
  list_plot = list()
  for(i in 1:length(choosen_col)){
    temp_plot = ggplot(Sto_temp, aes_string(x = "year", y = choosen_col[i]))+
      geom_line(aes(group = run_id), color = "blue", alpha = 0.1)+
      geom_line(data = Det_temp, mapping = aes_string(x = "year", y = choosen_col[i]), color = "red", size = 1)+
      facet_wrap(~country_name, scales = "free") +
      ggtitle(paste(choosen_col[i], i_age, sep = "_"))+
      theme(panel.background = element_blank(), axis.text.x = element_text(angle = -90))
    list_plot = c(list_plot, list(temp_plot))
  }
  return(list_plot)
}

for(index_age in 1:length(age_sel)){
  age_i_Sto_plotting_No_vac <- Sto_plotting_i_age(age_i_Sto_No_vac, i_age_Det_No_vac, age_sel[index_age])
  names(age_i_Sto_plotting_No_vac) = paste(choosen_col,age_sel[index_age], sep = "_")
  age_i_Sto_plotting_Cam_vac <- Sto_plotting_i_age(age_i_Sto_Cam_vac, i_age_Det_Cam_vac, age_sel[index_age])
  names(age_i_Sto_plotting_Cam_vac) = paste(choosen_col,age_sel[index_age], sep = "_")
  age_i_Sto_plotting_Rou_vac <- Sto_plotting_i_age(age_i_Sto_Rou_vac, i_age_Det_Rou_vac, age_sel[index_age])
  names(age_i_Sto_plotting_Rou_vac) = paste(choosen_col,age_sel[index_age], sep = "_")
  
  plot_save_path = paste("C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/result_model/first_model/Testing_template/test3(16 countries)/Det_vs_Sto model/Age_specific_Sto_vs_Det_rounded/",age_sel[index_age], sep = "")
  save_list_plot_func(age_i_Sto_plotting_No_vac)
  save_list_plot_func(age_i_Sto_plotting_Cam_vac)
  save_list_plot_func(age_i_Sto_plotting_Rou_vac)
}

###################Round function:
#round function:
sel_run_id = 1:200
Sto_template_save_path = "D:/OUCRU/Hannah/JEV/test3_16coun_200runid_every_countries/"
Sto_template_write_path = "D:/OUCRU/Hannah/JEV/test3_16coun_200runid_every_countries_rounded/"
template_file_name = "stochastic_burden_template_JE-OUCRU-Clapham_"

for(run in sel_run_id){
  print(run)
  Sto_No_vac_temp <- fread(paste(Sto_template_save_path, template_file_name, "je-routine-no-vaccination_",run,".csv", sep = ""))
  Sto_Cam_vac_temp <- fread(paste(Sto_template_save_path, template_file_name, "je-campaign-gavi_",run,".csv", sep = ""))
  Sto_Rou_vac_temp <- fread(paste(Sto_template_save_path, template_file_name, "je-routine-gavi_",run,".csv", sep = ""))
  
  Sto_No_vac_no_BTN <- filter(Sto_No_vac_temp, country != "BTN")
  Sto_No_vac_BTN <- filter(Sto_No_vac_temp, country == "BTN")
  Sto_Cam_vac_no_BTN <- filter(Sto_Cam_vac_temp, country != "BTN")
  Sto_Cam_vac_BTN <- filter(Sto_Cam_vac_temp, country == "BTN")
  Sto_Rou_vac_no_BTN <- filter(Sto_Rou_vac_temp, country != "BTN")
  Sto_Rou_vac_BTN <- filter(Sto_Rou_vac_temp, country == "BTN")
  
  Sto_No_vac_temp <- mutate(Sto_No_vac_no_BTN, deaths = round(deaths), cases = round(cases), dalys = round(dalys)) %>% rbind(Sto_No_vac_BTN, stringsAsFactors = F)
  Sto_No_vac_temp <- Sto_No_vac_temp[order(Sto_No_vac_temp$country_name, method = "radix"),]
  
  Sto_Cam_vac_temp <- mutate(Sto_Cam_vac_no_BTN, deaths = round(deaths), cases = round(cases), dalys = round(dalys)) %>% rbind(Sto_Cam_vac_BTN, stringsAsFactors = F)
  Sto_Cam_vac_temp <- Sto_Cam_vac_temp[order(Sto_Cam_vac_temp$country_name, method = "radix"),]
  
  Sto_Rou_vac_temp <- mutate(Sto_Rou_vac_no_BTN, deaths = round(deaths), cases = round(cases), dalys = round(dalys)) %>% rbind(Sto_Rou_vac_BTN, stringsAsFactors = F)
  Sto_Rou_vac_temp <- Sto_Rou_vac_temp[order(Sto_Rou_vac_temp$country_name, method = "radix"),]
  
  fwrite(Sto_No_vac_temp, file = paste(Sto_template_write_path, template_file_name, "je-routine-no-vaccination_",run,".csv", sep = ""))
  fwrite(Sto_Cam_vac_temp, file = paste(Sto_template_write_path, template_file_name, "je-campaign-gavi_",run,".csv", sep = ""))
  fwrite(Sto_Rou_vac_temp, file = paste(Sto_template_write_path, template_file_name, "je-routine-gavi_",run,".csv", sep = ""))
}
