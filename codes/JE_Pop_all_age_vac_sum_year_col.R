rm(list = ls())
###age group population interpolation => generate the Pop_all_age_year_sum column in country incidence data:
#library loading:
library(xlsx)
library(dplyr)
#data loading:
#Pop data load:
output_data_path = 'C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_clean/Montagu_data/test3/' #change this one
naive_pop = read.xlsx2(paste(output_data_path,'naive_pop_1950_2100.xlsx', sep = ""), 1, colClasses=NA)
data_inci_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_raw/cases_sero_data/" #change this one
year_sel = 1950:2100
age_group_sel = 0:99
#UN data for korea, taiwan, japan
UN_KOR_TWN_JAP_pop = read.xlsx2("C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_raw/UN-world-bank-data/UN_KOR_TWN_JAP_pop - Copy.xlsx",1, colClasses= NA) #change this one
Taiwan_data = UN_KOR_TWN_JAP_pop[1:65,paste("X", 0:99, sep = "")]
Japan_data = UN_KOR_TWN_JAP_pop[1:65 + 66,paste("X", 0:99, sep = "")]
Korean_data = UN_KOR_TWN_JAP_pop[1:65 + 66*2,paste("X", 0:99, sep = "")]
t(cbind(Taiwan_data, Japan_data, Korean_data))
convert_UN_data = data.frame(country = rep(c("TWN", "JPN", "KOR"), each = 100), age_from = rep(0:99,3), age_to = rep(0:99,3), 
                             t(cbind(Taiwan_data, Japan_data, Korean_data)))
colnames(convert_UN_data) = c("country", "age_from", "age_to", paste("X", 1951:2015, sep = ""))

#function to generate subnation population over time by the national population data and subnation population of a certain time:
sub_pop_gen <-function(sub_pop_year, year, naive_pop_country){
  naive_pop_country = select(naive_pop_country, starts_with("X"))
  naive_pop_prop_year = colSums(naive_pop_country)/colSums(naive_pop_country)[year]
  naive_pop_sub_year = naive_pop_prop_year*sub_pop_year
  naive_pop_prop_age = sweep(naive_pop_country, 2, colSums(naive_pop_country), FUN = "/")
  naive_pop_sub_age_year = sweep(naive_pop_prop_age, 2, naive_pop_sub_year, FUN = "*")
  return(naive_pop_sub_age_year)
}
#function to calculate vaccinated population from the other data:
pop_vac_gen <-function(vac_cam_data, sub_pop_gen){
  temp_sub_naive_pop = sub_pop_gen
  for(cam in 1:ncol(vac_cam_data)){
    temp_vac_camp_data = data.frame(matrix(0,ncol = ncol(naive_pop_less_15), nrow = nrow(naive_pop_less_15), 
                                           dimnames = list(0:14, colnames(naive_pop_less_15)))) #temp data to generate 1 campaign of vaccination.
    vac_1cam_pop = naive_pop_prop*vac_cam_data[1,cam] #interpolate subnation pop by the national pop
    sel_year_cam = names(vac_cam_data)[cam]
    ind_sel_year_pop = grep(sel_year_cam, colnames(sub_pop_gen))
    age_group = 1:15
    vac_pop_year = temp_sub_naive_pop[,ind_sel_year_pop]
    #Aging process:
    for(aging in 0:(nrow(naive_pop_less_15)-1)){
      if(aging == 0) {age_group_cam_aging = 1:nrow(naive_pop_less_15)}
      else {age_group_cam_aging = head(age_group_cam, - aging)}
      temp_vac_camp_data[age_group_cam_aging + aging,ind_sel_year_pop + aging] = vac_pop_year[age_group_cam_aging + aging]
    }
    temp_vac_camp_data = temp_vac_camp_data[, 1:ncol(naive_pop_less_15)]
    temp_sub_naive_pop = temp_sub_naive_pop - temp_vac_camp_data*vac_cam_data[2,cam]
  }
  temp_sub_naive_pop[temp_sub_naive_pop < 0] = 0
  vac_cam_pop = temp_sub_naive_pop
}
#function to calculate Pop_all_age_year_sum column:
pop_age_year_sum_gen <- function(country_study, country_pop){
  #Calculate the vector of years the study conducted and proportion of population that under that study period:
  temp_date_data = unlist(strsplit(as.character(country_study$Year[1]),"_"))
  start_year = as.numeric(substring(temp_date_data[1],4,7))
  start_year_prop = 1 - as.numeric(substring(temp_date_data[1],1,2))/12.0
  stop_year = as.numeric(substring(temp_date_data[2],4,7))
  stop_year_prop = as.numeric(substring(temp_date_data[2],1,2))/12.0
  year_study_vector = rep(1,stop_year - start_year + 1)
  if(length(year_study_vector) == 1){
    year_study_vector = stop_year_prop - (1 - start_year_prop)
  } else {
    year_study_vector[1] = start_year_prop
    year_study_vector[stop_year - start_year + 1] = stop_year_prop
  }
  names(year_study_vector) = paste("X",start_year:stop_year, sep = "")
  
  #sum population of years that effect by the study:
  effect_pop = country_pop[,names(year_study_vector)] * matrix(rep(year_study_vector, nrow(country_pop)), ncol = length(year_study_vector), byrow = T)
  effect_pop_sum = rowSums(effect_pop)
  
  #Calculate the age_l and age_u:
  age_group_split = strsplit(as.character(country_study$Age_group),"-")
  age_l = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[1])}))
  age_u = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[2])}))
  
  age_range_total = data.frame(character(),character())
  for(age in 1:length(age_l)){
    age_range = age_u[age] - age_l[age] + 1
    age_info = c(age_group_split[[age]][1], age_range)
    age_range_total = rbind(age_range_total, age_info, stringsAsFactors = F)
  }
  
  age_start = as.numeric(age_range_total[,1]) + 1
  
  pop_gen_case = c()
  for(age_s in 1:length(age_start)){
    age_stop = as.numeric(age_range_total[age_s, 2]) - 1
    sel_col_index = age_start[age_s]:(age_start[age_s] + age_stop)
    
    pop_gen_sel = as.numeric(as.character(effect_pop_sum))[sel_col_index]
    
    if(age_stop == 0){
      pop_gen_case = c(pop_gen_case,pop_gen_sel)
    } else {
      pop_gen_case = c(pop_gen_case,sum(pop_gen_sel))
    }
  }
  return(pop_gen_case)
}

#pop data:
IDM_data = read.xlsx2("C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_clean/pop_data/IDB_0_100_1981_2045_24.xlsx", 1, #change this one
                      colClasses = c(rep("character", 3), rep("numeric", 45)))
age_group_split = strsplit(as.character(IDM_data$Age_group),"-")
age_from = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[1])}))
age_to = age_from
IDM_data_reformat = data.frame(IDM_data$ISO, age_from, age_to, IDM_data[,4:48], stringsAsFactors = F)
colnames(IDM_data_reformat)[1] = "country"

#################################
#India:

#load population data of the country:
India_pop = naive_pop[naive_pop$country == "IND",]
#load age group from incidence data:
data_inci_India = read.xlsx(paste(data_inci_path, "India_JE.xlsx", sep = ""), 1)

#Gorakhpur (Uttar Pradesh) district: First study
#data from: http://www.census2011.co.in/census/district/559-gorakhpur.html 
Go_ce2011 = 4440895
India_Go_pop_gen = sub_pop_gen(Go_ce2011, "X2001", India_pop)
India_Gorakhpur_campaign_vac_pop = matrix(0, ncol = ncol(India_Go_pop_gen), nrow = nrow(India_Go_pop_gen), dimnames = list(NULL,colnames(India_Go_pop_gen)))
for(age in 0:(nrow(India_Go_pop_gen)-1)){
  temp_India_Gorakhpur_campaign_vac_pop = matrix(0, ncol = ncol(India_Go_pop_gen), nrow = nrow(India_Go_pop_gen), dimnames = list(NULL,colnames(India_Go_pop_gen)))
  temp_India_Gorakhpur_campaign_vac_pop[(1:15 + age),match("X2006", colnames(India_Go_pop_gen)) + age] = 1349047/15
  India_Gorakhpur_campaign_vac_pop = India_Gorakhpur_campaign_vac_pop + temp_India_Gorakhpur_campaign_vac_pop
}
India_Go_pop_gen = India_Go_pop_gen - India_Gorakhpur_campaign_vac_pop
India_Go_pop_gen[India_Go_pop_gen < 0] = 0
#calculate column:
India_subnation_index_1 = c(data_inci_India$Reference == unique(data_inci_India$Reference)[1])
India_study_1 = data_inci_India[India_subnation_index_1,]
India_subnation_index_1_1 = c(India_study_1$subnation == unique(India_study_1$subnation)[1])
India_study_1_1 = India_study_1[India_subnation_index_1_1,]
India_pop_all_age_year_sum_1_1 = pop_age_year_sum_gen(India_study_1_1, India_Go_pop_gen)
India_subnation_index_1_2 = c(India_study_1$subnation == unique(India_study_1$subnation)[2])
India_study_1_2 = India_study_1[India_subnation_index_1_2,]
India_pop_all_age_year_sum_1_2 = pop_age_year_sum_gen(India_study_1_2, India_Go_pop_gen)
India_pop_all_age_year_sum_1 = c(India_pop_all_age_year_sum_1_1, India_pop_all_age_year_sum_1_2)

#Uttar Pradesh: Second study => before any vaccination program
Utt_cen_all_year = sub_pop_gen(199812341, "X2011", India_pop)
India_subnation_index_2 = c(data_inci_India$Reference == unique(data_inci_India$Reference)[2])
India_study_2 = data_inci_India[India_subnation_index_2,]
India_pop_all_age_year_sum_2 = pop_age_year_sum_gen(India_study_2, Utt_cen_all_year)

#Get vaccination pop in India from the gavi data => for Third and Fifth study:
#get the function of generating campaign and routine population:
source("C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/code_JEV_model/JEV_Montagu_data - Scen_pop_generate_function.R") #change this one
India_pop_cam_vac = cam_scen_pop_gen("IND", India_pop, 1950:2100)

#Tamil Nadu: Third study
TaNa_cen_all_year = select(India_pop_cam_vac, starts_with("X"))*72147030/sum(India_pop$X2011)
India_subnation_index_3 = c(data_inci_India$Reference == unique(data_inci_India$Reference)[3])
India_study_3 = data_inci_India[India_subnation_index_3,]
#cohort study:
India_pop_all_age_year_sum_3 = c()
for(i in 1:4){
  India_subnation_index_3_i = c(India_study_3$subnation == unique(India_study_3$subnation)[i])
  India_study_3_i = India_study_3[India_subnation_index_3_i,]
  India_pop_all_age_year_sum_3_temp = pop_age_year_sum_gen(India_study_3_i, TaNa_cen_all_year)
  India_pop_all_age_year_sum_3 = c(India_pop_all_age_year_sum_3, India_pop_all_age_year_sum_3_temp)
}

#Assam: Forth Study => before any vaccination program
#population of Assam:
Assam_pop_2011 = 31205576 #data from http://www.census2011.co.in/census/state/assam.html
India_Assam_pop_gen = sub_pop_gen(Assam_pop_2011, "X2011", India_pop)
India_subnation_index_4 = c(data_inci_India$Reference == unique(data_inci_India$Reference)[4])
India_study_4 = data_inci_India[India_subnation_index_4,]
#cohort study:
India_pop_all_age_year_sum_4 = c()
for(i in 1:4){
  India_subnation_index_4_i = c(India_study_4$subnation == unique(India_study_4$subnation)[i])
  India_study_4_i = India_study_4[India_subnation_index_4_i,]
  India_pop_all_age_year_sum_4_temp = pop_age_year_sum_gen(India_study_4_i, India_Assam_pop_gen)
  India_pop_all_age_year_sum_4 = c(India_pop_all_age_year_sum_4, India_pop_all_age_year_sum_4_temp)
}

#Fifth Study: 7 districts Medical College Hospital Assam
# 7 upper districts of Assam: Jorhat, Dibrugarh, Dhemaji, Golaghat, Lakhimpur, Sivasagar, and Tinsukia
#+ Arunachal Pradesh + Nagaland
Assam_7districts_2011 = 1092256 + 1326335 + 686133 + 1066888 + 1042137 + 1151050 + 1327929 #data from https://www.mapsofindia.com/maps/assam/assam-district.htm
Assam_7districts_2011 = Assam_7districts_2011 + 1383727 + 1978502#data from : http://www.census2011.co.in/census/state/
Assam_7districts_cen_all_year = select(India_pop_cam_vac, starts_with("X"))*Assam_7districts_2011/sum(India_pop$X2011)
India_subnation_index_5 = c(data_inci_India$Reference == unique(data_inci_India$Reference)[5])
India_study_5 = data_inci_India[India_subnation_index_5,]
India_pop_all_age_year_sum_5 = pop_age_year_sum_gen(India_study_5, India_Assam_pop_gen)

#Sixfth Study: northern dítricts of West Bengal:
#northern dítricts of West Bengal: Jalpaiguri, Darjeeling, Dakshin Dinajpur, Uttar Dinajpur, and Cooch Behar; pop data from https://en.wikipedia.org/wiki/West_Bengal
Northern_dis_West_Bengal_2017 = c(3872846, 1846823, 1676276, 3007134, 2819086)
Northern_dis_West_Bengal_cam_vac_2013 = c(0.7893, 0.8038, 0.761, 0.6789, 0)
Northern_dis_West_Bengal_rou_vac_2013 = c(0.90, 0.90, 0.90, 0.90, 0)
Northern_dis_West_Bengal_cen_all_year = lapply(Northern_dis_West_Bengal_2017, function(x)sub_pop_gen(x, "X2017",India_pop))
#Campaign vac in 2 to 15 years old, Rou vac in 1 year old:
Northern_dis_West_Bengal_after_vac_2013_2017 = list()
for(dist in 1:length(Northern_dis_West_Bengal_2017)){
  Northern_dis_West_Bengal_dist_pop = Northern_dis_West_Bengal_cen_all_year[[dist]]
  Northern_dis_West_Bengal_year_vac = paste("X", 2013:2017, sep = "")
  Northern_dis_West_Bengal_cam_pop = Northern_dis_West_Bengal_cam_vac_2013[dist]*Northern_dis_West_Bengal_dist_pop[1:15,"X2013"]
  Northern_dis_West_Bengal_rou_pop = Northern_dis_West_Bengal_rou_vac_2013[dist]*Northern_dis_West_Bengal_dist_pop[1,Northern_dis_West_Bengal_year_vac]
  Northern_dis_West_Bengal_dist_after_vac = Northern_dis_West_Bengal_dist_pop
  #aging process:
  for(i in 1:length(Northern_dis_West_Bengal_year_vac)){
    Northern_dis_West_Bengal_dist_after_vac[1:15 - 1 + i, Northern_dis_West_Bengal_year_vac[i]] = Northern_dis_West_Bengal_dist_pop[1:15, Northern_dis_West_Bengal_year_vac[1]] - Northern_dis_West_Bengal_cam_pop
    Northern_dis_West_Bengal_dist_after_vac[1:i, Northern_dis_West_Bengal_year_vac[i]] <- unlist(Northern_dis_West_Bengal_dist_after_vac[1:i, Northern_dis_West_Bengal_year_vac[i]] - Northern_dis_West_Bengal_rou_pop[1:i])
  }
  Northern_dis_West_Bengal_dist_after_vac[Northern_dis_West_Bengal_dist_after_vac < 0] = 0
  Northern_dis_West_Bengal_after_vac_2013_2017[[dist]] = Northern_dis_West_Bengal_dist_after_vac
}
Northern_dis_West_Bengal_all_after_vac <- Reduce("+", Northern_dis_West_Bengal_after_vac_2013_2017)
India_subnation_index_6 = c(data_inci_India$Reference == unique(data_inci_India$Reference)[6])
India_study_6 = data_inci_India[India_subnation_index_6,]
India_pop_all_age_year_sum_6 = pop_age_year_sum_gen(India_study_6, Northern_dis_West_Bengal_all_after_vac)

#Seventh Study: Uttar Pradesh after vaccination: campaign vaccination in 2006 -> 2009
Utt_cen_all_year = sub_pop_gen(199812341, "X2011", India_pop)
#Campaign vaccination from 2006-2009, 1-15 years old. data from "A review of Japanese Encephalitis in Uttar Pradesh, India"
Utt_cam_vac_2006_2009 <- c(6836506, 9499157, 10708393, 7831079)
names(Utt_cam_vac_2006_2009) <- paste("X", 2006:2009, sep = "")
Utt_cen_after_cam_vac <- Utt_cen_all_year
for(utt_cam_index in 1:length(Utt_cam_vac_2006_2009)){
  Utt_year_cam_vac_to_2017 <- paste("X", 2006:2017, sep = "")[utt_cam_index:length(2006:2017)]
  for(utt_cam in 1:length(Utt_year_cam_vac_to_2017)){
    Utt_cen_after_cam_vac[1:15 - 1 + utt_cam, Utt_year_cam_vac_to_2017[utt_cam]] = Utt_cen_all_year[1:15 - 1 + utt_cam, Utt_year_cam_vac_to_2017[utt_cam]] - Utt_cam_vac_2006_2009[utt_cam_index]/15
  }
}
India_subnation_index_7 = c(data_inci_India$Reference == unique(data_inci_India$Reference)[7])
India_study_7 = data_inci_India[India_subnation_index_7,]
India_pop_all_age_year_sum_7 = c()
for(i in 1:4){
  India_subnation_index_7_i = c(India_study_7$subnation == unique(India_study_7$subnation)[i])
  India_study_7_i = India_study_7[India_subnation_index_7_i,]
  India_pop_all_age_year_sum_7_temp = pop_age_year_sum_gen(India_study_7_i, India_pop)
  India_pop_all_age_year_sum_7 = c(India_pop_all_age_year_sum_7, India_pop_all_age_year_sum_7_temp)
}

#Eight study: Kushinagar in Uttar Pradesh
Kushinagar_cen_all_year = sub_pop_gen(3560830, "X2011", India_pop) #data from https://en.wikipedia.org/wiki/Kushinagar_district
#Campaign vaccinate 1 085 055 children in 2006, 1-15 years old. data from "A review of Japanese Encephalitis in Uttar Pradesh, India"
Kushinagar_cen_after_cam_vac <- Kushinagar_cen_all_year
Kushinagar_year_cam_vac_to_2017 <- paste("X", 2006:2017, sep = "")
  for(Kushinagar_cam in 1:length(Kushinagar_year_cam_vac_to_2017)){
    Kushinagar_cen_after_cam_vac[1:15 - 1 + Kushinagar_cam, Kushinagar_year_cam_vac_to_2017[Kushinagar_cam]] = Kushinagar_cen_all_year[1:15 - 1 + Kushinagar_cam, Kushinagar_year_cam_vac_to_2017[Kushinagar_cam]] - 1085055/sum(Kushinagar_cen_all_year[1:15,"X2006"])*Kushinagar_cen_all_year[1:15,"X2006"]
  }
Kushinagar_cen_after_cam_vac[Kushinagar_cen_after_cam_vac < 0] <- 0
India_subnation_index_8 = c(data_inci_India$Reference == unique(data_inci_India$Reference)[8])
India_study_8 = data_inci_India[India_subnation_index_8,]
India_pop_all_age_year_sum_8 = pop_age_year_sum_gen(India_study_8, Kushinagar_cen_after_cam_vac)

#Nineth study: Gorakhpur district 2011, 2013
India_subnation_index_9 = c(data_inci_India$Reference == unique(data_inci_India$Reference)[9])
India_study_9 = data_inci_India[India_subnation_index_9,]
India_pop_all_age_year_sum_9 = pop_age_year_sum_gen(India_study_9, India_Go_pop_gen)

#Tenth study: Nothern district Uttar Pradesh 2011-2012 - 34 districts. Data from https://www.mapsofindia.com/maps/uttarpradesh/uttar-pradesh-district.htm
Northern_Utt_2011 <- c(4589838, 4483992, 3260699, 3405559, 3108367, 4092845, 4021243, 3487731, 3433919, 2470996, 3797117, 2397888, 2464464, 2559297, 4440895)
Northern_Utt_pop_gen = sum(Northern_Utt_2011)/sum(Utt_cen_all_year[,"X2011"])*Utt_cen_after_cam_vac
India_subnation_index_10 = c(data_inci_India$Reference == unique(data_inci_India$Reference)[10])
India_study_10 = data_inci_India[India_subnation_index_10,]
India_pop_all_age_year_sum_10 = pop_age_year_sum_gen(India_study_10, Northern_Utt_pop_gen)

#Eleventh study: Gorakhpur district 2005
India_subnation_index_11 = c(data_inci_India$Reference == unique(data_inci_India$Reference)[11])
India_study_11 = data_inci_India[India_subnation_index_11,]
India_pop_all_age_year_sum_11 = pop_age_year_sum_gen(India_study_11, India_Go_pop_gen)

#Twelfth study: Cuddalore district 2002-2003
Cuddalore_2011 <- 2605914 #data from: https://en.wikipedia.org/wiki/Cuddalore_district
Cuddalore_pop_gen = sub_pop_gen(Cuddalore_2011, "X2011", India_pop)
India_subnation_index_12 = c(data_inci_India$Reference == unique(data_inci_India$Reference)[12])
India_study_12 = data_inci_India[India_subnation_index_12,]
India_pop_all_age_year_sum_12 = pop_age_year_sum_gen(India_study_12, Cuddalore_pop_gen)

#Thirdteenth study: Pondicherry 2001
Pondicherry_2011 <- 654392 #data from: https://en.wikipedia.org/wiki/Pondicherry_district
Pondicherry_pop_gen = sub_pop_gen(Pondicherry_2011, "X2011", India_pop)
India_subnation_index_13 = c(data_inci_India$Reference == unique(data_inci_India$Reference)[13])
India_study_13 = data_inci_India[India_subnation_index_13,]
India_pop_all_age_year_sum_13 = pop_age_year_sum_gen(India_study_13, Pondicherry_pop_gen)

#Fourteenth study: Tamil Nadu 1986-1990
Tamil_Nadu_1987 <- 120000 #cathment area based on the paper
Tamil_Nadu_pop_gen = sub_pop_gen(Tamil_Nadu_1987, "X1987", India_pop)
India_subnation_index_14 = c(data_inci_India$Reference == unique(data_inci_India$Reference)[14])
India_study_14 = data_inci_India[India_subnation_index_14,]
India_pop_all_age_year_sum_14 = pop_age_year_sum_gen(India_study_14, Tamil_Nadu_pop_gen)

#Fifteenth study: Dhemaji district
Dhemaji_cen_all_year = sub_pop_gen(688077, "X2011", India_pop) #data from the paper
#Campaign vaccinate in 2008, 1-15 years old. data from the paper
Dhemaji_cen_after_cam_vac <- Dhemaji_cen_all_year
Dhemaji_year_cam_vac_to_2017 <- paste("X", 2008:2017, sep = "")
for(Dhemaji_cam in 1:length(Dhemaji_year_cam_vac_to_2017)){
  Dhemaji_cen_after_cam_vac[1:15 - 1 + Dhemaji_cam, Dhemaji_year_cam_vac_to_2017[Dhemaji_cam]] = 
    Dhemaji_cen_all_year[1:15 - 1 + Dhemaji_cam, Dhemaji_year_cam_vac_to_2017[Dhemaji_cam]] - Dhemaji_cen_all_year[1:15,"X2008"]*0.9
}
Dhemaji_cen_after_cam_vac[Dhemaji_cen_after_cam_vac < 0] <- 0
India_subnation_index_15 = c(data_inci_India$Reference == unique(data_inci_India$Reference)[15])
India_study_15 = data_inci_India[India_subnation_index_15,]
India_pop_all_age_year_sum_15 = pop_age_year_sum_gen(India_study_15, Dhemaji_cen_after_cam_vac)

#Sixteenth study: Bellary 1995-1997 and neighbor: Raichur, Chitradurga, Kurnool, and Anantapur:
Bellary_1996 <- 2452595 + 1928812 + 1659456 + 4046601 + 4083315# data from this data: https://en.wikipedia.org/wiki/List_of_districts_of_Karnataka ; https://en.wikipedia.org/wiki/Kurnool_district ; https://en.wikipedia.org/wiki/Anantapur_district 
Bellary_pop_gen = sub_pop_gen(Bellary_1996, "X1996", India_pop)
India_subnation_index_16 = c(data_inci_India$Reference == unique(data_inci_India$Reference)[16])
India_study_16 = data_inci_India[India_subnation_index_16,]
India_pop_all_age_year_sum_16 = pop_age_year_sum_gen(India_study_16, Bellary_pop_gen)

#create the Pop_all_age_year_sum column:
data_inci_India$Pop_all_age_year_sum = c(India_pop_all_age_year_sum_1, India_pop_all_age_year_sum_2, India_pop_all_age_year_sum_3, 
                                         India_pop_all_age_year_sum_4, India_pop_all_age_year_sum_5, India_pop_all_age_year_sum_6,
                                         India_pop_all_age_year_sum_7, India_pop_all_age_year_sum_8, India_pop_all_age_year_sum_9,
                                         India_pop_all_age_year_sum_10, India_pop_all_age_year_sum_11, India_pop_all_age_year_sum_12,
                                         India_pop_all_age_year_sum_13, India_pop_all_age_year_sum_14, India_pop_all_age_year_sum_15, India_pop_all_age_year_sum_16)
write.xlsx(data_inci_India, paste(data_inci_path, "India_JE.xlsx", sep = ""), row.names = F)


##################################
#Japan:
#Japan census data:
Japan_pop = IDM_data_reformat[IDM_data_reformat$country == "JAP", ]
Japan_pop[,paste("X", 1976:1980, sep = "")] = Japan_pop$X1981
Japan_pop_gen = Japan_pop[,c("country", "age_from", "age_to",sort(colnames(Japan_pop[,paste("X", 1976:2025, sep = "")])))]
#japan routine vaccination from 1976 to 2005 with about 50% coverage:
Japan_National_routine_vac_pop = matrix(0, ncol = ncol(Japan_pop_gen), nrow = nrow(Japan_pop_gen), dimnames = list(NULL,colnames(Japan_pop_gen)))
Japan_routine_vac_age_0_from_1976_to_2005 =c(rep(0.5,length(1976:1994)),rep(0.8,length(1995:2005)))*Japan_pop_gen[1,paste("X", 1976:2005, sep = "")]
Jap_year_boundary = 0
for(age in 0:(nrow(Japan_pop_gen)-1)){
  temp_Japan_National_routine_vac_pop = matrix(0, ncol = ncol(Japan_pop_gen), nrow = nrow(Japan_pop_gen), dimnames = list(NULL,colnames(Japan_pop_gen)))
  if (age <= 20){
    temp_Japan_National_routine_vac_pop[(1 + age),paste("X", 1976:2005 + age, sep = "")] = unlist(Japan_routine_vac_age_0_from_1976_to_2005)
  } else {
    Jap_year_boundary = Jap_year_boundary + 1
    temp_Japan_National_routine_vac_pop[(1 + age),head(paste("X", 1976:2005 + age, sep = ""),- Jap_year_boundary)] = head(unlist(Japan_routine_vac_age_0_from_1976_to_2005), - Jap_year_boundary)
    
  }
  #temp_Japan_National_routine_vac_pop = temp_Japan_National_routine_vac_pop[,1:ncol(Japan_National_routine_vac_pop)]
  Japan_National_routine_vac_pop = Japan_National_routine_vac_pop + temp_Japan_National_routine_vac_pop
}
Japan_pop_gen = Japan_pop_gen - Japan_National_routine_vac_pop
Japan_pop_gen[Japan_pop_gen < 0] = 0


#get age group data from incidence data:
data_inci_Jap = read.xlsx(paste(data_inci_path, "Japan_JE.xlsx", sep = ""), 1)

#first study: National
Jap_subnation_index_1 = c(data_inci_Jap$subnation == unique(data_inci_Jap$subnation)[1])
Jap_study_1 = data_inci_Jap[Jap_subnation_index_1,]
Jap_pop_all_age_year_sum_1 = pop_age_year_sum_gen(Jap_study_1, Japan_pop_gen)

#Second study: Chugoku
Chugoku_pop_2010 = 7563428 #data from: https://en.wikipedia.org/wiki/Ch%C5%ABgoku_region
Chugoku_pop_gen = sub_pop_gen(Chugoku_pop_2010, "X2010", Japan_pop)
Jap_subnation_index_2 = c(data_inci_Jap$subnation == unique(data_inci_Jap$subnation)[2])
Jap_study_2 = data_inci_Jap[Jap_subnation_index_2,]
Jap_pop_all_age_year_sum_2 = pop_age_year_sum_gen(Jap_study_2, Japan_pop)

#Third study: National
Jap_subnation_index_3 = c(data_inci_Jap$subnation == unique(data_inci_Jap$subnation)[3])
Jap_study_3 = data_inci_Jap[Jap_subnation_index_3,]
Jap_pop_all_age_year_sum_3 = pop_age_year_sum_gen(Jap_study_3, Japan_pop_gen)

#Forth study: National
Jap_subnation_index_4 = c(data_inci_Jap$subnation == unique(data_inci_Jap$subnation)[4])
Jap_study_4 = data_inci_Jap[Jap_subnation_index_4,]
Jap_pop_all_age_year_sum_4 = pop_age_year_sum_gen(Jap_study_4, Japan_pop_gen)

write.xlsx(data_inci_Jap, paste(data_inci_path, "Japan_JE.xlsx", sep = ""), row.names = F)

##################################
#Taiwan:
Taiwan_pop = convert_UN_data[convert_UN_data$country == "TWN", ]
Taiwan_pop = Taiwan_pop[, paste("X",1951:2015, sep = "")]
Taiwan_pop_gen = Taiwan_pop

#Vaccination program in Taiwan: campaign vaccination in 1968 for 0-3 year old, routine vaccination followed that time
#campaign vaccination
Taiwan_National_campaign_vac_pop = data.frame(matrix(0, ncol = ncol(Taiwan_pop_gen), nrow = nrow(Taiwan_pop_gen), dimnames = list(NULL,colnames(Taiwan_pop_gen))))
for(age in 0:(nrow(Taiwan_pop_gen)-1)){
  temp_Taiwan_National_campaign_vac_pop = data.frame(matrix(0, ncol = ncol(Taiwan_pop_gen), nrow = nrow(Taiwan_pop_gen), dimnames = list(NULL,colnames(Taiwan_pop_gen))))
  temp_Taiwan_National_campaign_vac_pop[(1:4 + age),paste("X",1968 + age, sep = "")] = Taiwan_pop_gen[1:4,"X1968"]
  Taiwan_National_campaign_vac_pop = Taiwan_National_campaign_vac_pop + temp_Taiwan_National_campaign_vac_pop
}
Taiwan_pop_gen = Taiwan_pop_gen - Taiwan_National_campaign_vac_pop
Taiwan_pop_gen[Taiwan_pop_gen < 0] = 0

#routine vaccination:
Taiwan_National_routine_vac_pop = matrix(0, ncol = ncol(Taiwan_pop_gen), nrow = nrow(Taiwan_pop_gen), dimnames = list(NULL,colnames(Taiwan_pop_gen)))
Taiwan_routine_vac_age_0_from_1969_to_2015 = rep(0.99, length(1969:2015))*Taiwan_pop_gen[1,paste("X", 1969:2015, sep = "")]
Taiwan_rou_year_boundary = 0
for(age in 0:(nrow(Taiwan_pop_gen)-1)){
  temp_Taiwan_National_routine_vac_pop = matrix(0, ncol = ncol(Taiwan_pop_gen), nrow = nrow(Taiwan_pop_gen), dimnames = list(NULL,colnames(Taiwan_pop_gen)))
  if(age == 0){
    temp_Taiwan_National_routine_vac_pop[1,match("X1969",colnames(Taiwan_pop_gen)):match("X2015",colnames(Taiwan_pop_gen))] = unlist(Taiwan_routine_vac_age_0_from_1969_to_2015)
  } else {
    temp_Taiwan_National_routine_vac_pop[(1 + age),head(match("X1969",colnames(Taiwan_pop_gen)):match("X2015",colnames(Taiwan_pop_gen)) + age, - Taiwan_rou_year_boundary)] = head(unlist(Taiwan_routine_vac_age_0_from_1969_to_2015), -Taiwan_rou_year_boundary)
  }
  
  Taiwan_rou_year_boundary = Taiwan_rou_year_boundary + 1
  Taiwan_National_routine_vac_pop = Taiwan_National_routine_vac_pop + temp_Taiwan_National_routine_vac_pop
}
Taiwan_pop_gen = Taiwan_pop_gen - Taiwan_National_routine_vac_pop
Taiwan_pop_gen[Taiwan_pop_gen < 0] = 0

#get age group data from incidence data:
data_inci_Taiwan = read.xlsx(paste(data_inci_path, "Taiwan_JE.xlsx", sep = ""), 1)

#First study: cohort study:
Taiwan_pop_all_age_year_sum_1 = c()
Taiwan_subnation_index_1 = c(data_inci_Taiwan$Reference == unique(data_inci_Taiwan$Reference)[1])
Taiwan_study_1 = data_inci_Taiwan[Taiwan_subnation_index_1,] 
for(i in 1:12){
  Taiwan_subnation_index_1_i = c(Taiwan_study_1$subnation == unique(Taiwan_study_1$subnation)[i])
  Taiwan_study_1_i = Taiwan_study_1[Taiwan_subnation_index_1_i,]
  Taiwan_pop_all_age_year_sum_1_temp = pop_age_year_sum_gen(Taiwan_study_1_i, Taiwan_pop_gen)
  Taiwan_pop_all_age_year_sum_1 = c(Taiwan_pop_all_age_year_sum_1, Taiwan_pop_all_age_year_sum_1_temp)
}

#Second study: cohort study:
Taiwan_pop_all_age_year_sum_2 = c()
Taiwan_subnation_index_2 = c(data_inci_Taiwan$Reference == unique(data_inci_Taiwan$Reference)[2])
Taiwan_study_2 = data_inci_Taiwan[Taiwan_subnation_index_2,] 
for(i in 1:length(unique(Taiwan_study_2$subnation))){
  Taiwan_subnation_index_2_i = c(Taiwan_study_2$subnation == unique(Taiwan_study_2$subnation)[i])
  Taiwan_study_2_i = Taiwan_study_2[Taiwan_subnation_index_2_i,]
  Taiwan_pop_all_age_year_sum_2_temp = pop_age_year_sum_gen(Taiwan_study_2_i, Taiwan_pop_gen)
  Taiwan_pop_all_age_year_sum_2 = c(Taiwan_pop_all_age_year_sum_2, Taiwan_pop_all_age_year_sum_2_temp)
}

#Third study: cohort study:
Taiwan_pop_all_age_year_sum_3 = c()
Taiwan_subnation_index_3 = c(data_inci_Taiwan$Reference == unique(data_inci_Taiwan$Reference)[3])
Taiwan_study_3 = data_inci_Taiwan[Taiwan_subnation_index_3,] 
for(i in 1:length(unique(Taiwan_study_3$subnation))){
  Taiwan_subnation_index_3_i = c(Taiwan_study_3$subnation == unique(Taiwan_study_3$subnation)[i])
  Taiwan_study_3_i = Taiwan_study_3[Taiwan_subnation_index_3_i,]
  Taiwan_pop_all_age_year_sum_3_temp = pop_age_year_sum_gen(Taiwan_study_3_i, Taiwan_pop_gen)
  Taiwan_pop_all_age_year_sum_3 = c(Taiwan_pop_all_age_year_sum_3, Taiwan_pop_all_age_year_sum_3_temp)
}

#create the Pop_all_age_year_sum column:
data_inci_Taiwan$Pop_all_age_year_sum = c(Taiwan_pop_all_age_year_sum_1, Taiwan_pop_all_age_year_sum_2, Taiwan_pop_all_age_year_sum_3)
write.xlsx(data_inci_Taiwan, paste(data_inci_path, "Taiwan_JE.xlsx", sep = ""), row.names = F)

##################################
#Cambodia:
#Cambodia census data:
Cam_pop = naive_pop[naive_pop$country == "KHM", ]
#add vaccination information: vaccinate 12500 for each of 00-00 and 01-01 age group in 2009 and 2010:
Cam_vac_pop = matrix(0, ncol = ncol(Cam_pop), nrow = nrow(Cam_pop), dimnames = list(NULL,colnames(Cam_pop)))
for(age in 0:(nrow(Cam_pop)-1)){
  temp_Cam_vac_pop = matrix(0, ncol = ncol(Cam_pop), nrow = nrow(Cam_pop), dimnames = list(NULL,colnames(Cam_pop)))
  temp_Cam_vac_pop[(1:2 + age),paste("X",2009:2010 + age, sep = "")] = 12500
  temp_Cam_vac_pop = temp_Cam_vac_pop[,1:ncol(Cam_vac_pop)]
  Cam_vac_pop = Cam_vac_pop + temp_Cam_vac_pop
}
Cam_pop_vaccinated = Cam_pop - Cam_vac_pop

#get age group data from incidence data:
data_inci_Cam = read.xlsx(paste(data_inci_path, "Cambodia_JE.xlsx", sep = ""), 1)
#first study: before the vaccination program:
Cam_subnation_index_1 = c(data_inci_Cam$subnation == unique(data_inci_Cam$subnation)[1])
Cam_study_1 = data_inci_Cam[Cam_subnation_index_1,]
Cam_pop_all_age_year_sum_1 = pop_age_year_sum_gen(Cam_study_1, Cam_pop_vaccinated)
#second study:
#population of siemreap: No vaccination program
Siemreap_pop_2017 = 139458 #data from: http://worldpopulationreview.com/countries/cambodia-population/cities/
Cam_Siemreap_pop_gen = sub_pop_gen(Siemreap_pop_2017, "X2017", Cam_pop)
Cam_subnation_index_2 = c(data_inci_Cam$subnation == unique(data_inci_Cam$subnation)[2])
Cam_study_2 = data_inci_Cam[Cam_subnation_index_2,]
Cam_pop_all_age_year_sum_2 = pop_age_year_sum_gen(Cam_study_2, Cam_Siemreap_pop_gen)
#third study: National scale, after vaccination program
Cam_subnation_index_3 = c(data_inci_Cam$subnation == unique(data_inci_Cam$subnation)[3])
Cam_study_3 = data_inci_Cam[Cam_subnation_index_3,]
Cam_pop_all_age_year_sum_3 = pop_age_year_sum_gen(Cam_study_3, Cam_pop_vaccinated)
#fourth study: grouping ages from the first study, before the vaccination program
Cam_subnation_index_4 = c(data_inci_Cam$subnation == unique(data_inci_Cam$subnation)[4])
Cam_study_4 = data_inci_Cam[Cam_subnation_index_4,]
Cam_pop_all_age_year_sum_4 = pop_age_year_sum_gen(Cam_study_4, Cam_pop)

#create the Pop_all_age_year_sum column:
data_inci_Cam$Pop_all_age_year_sum = c(Cam_pop_all_age_year_sum_1, Cam_pop_all_age_year_sum_2, Cam_pop_all_age_year_sum_3, Cam_pop_all_age_year_sum_4)
write.xlsx(data_inci_Cam, paste(data_inci_path, "Cambodia_JE.xlsx", sep = ""), row.names = F)

##################################
#Indonesia:
#Indonesia census data:
Indo_pop = naive_pop[naive_pop$country == "IDN", ]

#get age group data from incidence data:
data_inci_Indo = read.xlsx(paste(data_inci_path, "Indonesia_JE.xlsx", sep = ""), 1)

#first study:
Indo_subnation_index_1 = c(data_inci_Indo$subnation == unique(data_inci_Indo$subnation)[1])
Indo_study_1 = data_inci_Indo[Indo_subnation_index_1,]
Indo_pop_all_age_year_sum_1 = pop_age_year_sum_gen(Indo_study_1, Indo_pop)

#second study: 
#6 provinces of the dataset: East Java, East Nusa Tenggara (West Timor)(2010 only), Papua, West Kalimantan, West Nusa Tenggara (Lombok Island), West Sumatera
Indo_6_prov_pop_2014 = 38529481 + 1662056 + 3486432 + 4073304 + 3311044 + 	5131882
Indo_6_prov_pop_gen = sub_pop_gen(Indo_6_prov_pop_2014, "X2014", Indo_pop)
Indo_subnation_index_2 = c(data_inci_Indo$subnation == unique(data_inci_Indo$subnation)[2])
Indo_study_2 = data_inci_Indo[Indo_subnation_index_2,]
Indo_pop_all_age_year_sum_2 = pop_age_year_sum_gen(Indo_study_2, Indo_6_prov_pop_gen)

#create the Pop_all_age_year_sum column:
data_inci_Indo$Pop_all_age_year_sum = c(Indo_pop_all_age_year_sum_1, Indo_pop_all_age_year_sum_2)
write.xlsx(data_inci_Indo, paste(data_inci_path, "Indonesia_JE.xlsx", sep = ""), row.names = F)

###################################
#Laos:
#Laos census data:
Laos_pop = naive_pop[naive_pop$country == "LAO", ]
#population of Vientian:
Vientian_pop_2017 = 196731 #data from: http://worldpopulationreview.com/countries/laos-population/cities/
Laos_Vientian_pop_gen = sub_pop_gen(Vientian_pop_2017, "X2017", Laos_pop)
#get age group data from incidence data:
data_inci_Laos = read.xlsx(paste(data_inci_path, "Laos_JE.xlsx", sep = ""), 1)
Laos_subnation_index_1 = c(data_inci_Laos$subnation == unique(data_inci_Laos$subnation)[1])
Laos_study_1 = data_inci_Laos[Laos_subnation_index_1,]
Laos_pop_all_age_year_sum_1 = pop_age_year_sum_gen(Laos_study_1, Laos_Vientian_pop_gen)

#create the Pop_all_age_year_sum column:
data_inci_Laos$Pop_all_age_year_sum = c(Laos_pop_all_age_year_sum_1)
write.xlsx(data_inci_Laos, paste(data_inci_path, "Laos_JE.xlsx", sep = ""), row.names = F)

##################################
#South Korea:
#South Korea census data:
S_Korean_pop = convert_UN_data[convert_UN_data$country == "KOR", ]
S_Korean_pop = S_Korean_pop[, paste("X",1951:2015, sep = "")]
S_Korean_pop$X2016 = S_Korean_pop$X2015
S_Korean_pop_gen = S_Korean_pop

#Vaccination program in S_Korean: from 1981 to 1998 - campaign vaccination for 3 to 15 years old, distribution based on the graph:
#campaign vaccination
S_Korean_National_campaign_vac_pop = data.frame(matrix(0, ncol = ncol(S_Korean_pop_gen), nrow = nrow(S_Korean_pop_gen), dimnames = list(NULL,colnames(S_Korean_pop_gen))))
S_Korean_National_campaign_distribution_1981_1998 = c(1461738,1607554,2582605,1996664,12227829,11385806,11385215,11396884,11018276,10017717,10017192,
                                                      8821511,8162448,4186280,4185655,4209421,4184538,4184044)
S_Korean_National_campaign_distribution_1981_1998_matrix = matrix(rep(S_Korean_National_campaign_distribution_1981_1998/length(4:16),length(4:16)), nrow = length(4:16), byrow = T)
S_Korean_cam_year_boundary = 0
for(age in 0:(nrow(S_Korean_pop_gen)-1)){
  temp_S_Korean_National_campaign_vac_pop = data.frame(matrix(0, ncol = ncol(S_Korean_pop_gen), nrow = nrow(S_Korean_pop_gen), dimnames = list(NULL,colnames(S_Korean_pop_gen))))
  year_index_vac = match(paste("X", 1981:1998, sep = ""),colnames(S_Korean_pop_gen)) + age
  if(age<18){
    temp_S_Korean_National_campaign_vac_pop[(4:16 + age),year_index_vac] = S_Korean_National_campaign_distribution_1981_1998_matrix
  } else {
    S_Korean_cam_year_boundary = S_Korean_cam_year_boundary + 1
    temp_S_Korean_National_campaign_vac_pop[(4:16 + age),head(year_index_vac, - S_Korean_cam_year_boundary)] = S_Korean_National_campaign_distribution_1981_1998_matrix[,head(1:ncol(S_Korean_National_campaign_distribution_1981_1998_matrix), - S_Korean_cam_year_boundary)]
  }
  S_Korean_National_campaign_vac_pop = S_Korean_National_campaign_vac_pop + temp_S_Korean_National_campaign_vac_pop
}
S_Korean_pop_gen = S_Korean_pop_gen - S_Korean_National_campaign_vac_pop
S_Korean_pop_gen[S_Korean_pop_gen < 0] = 0

#routine vaccination:
#Since 2000, public health centers have provided a three-dose primary series of inactivated JE vaccine starting at 12 months of age
S_Korean_National_routine_vac_pop = matrix(0, ncol = ncol(S_Korean_pop_gen), nrow = nrow(S_Korean_pop_gen), dimnames = list(NULL,colnames(S_Korean_pop_gen)))
S_Korean_routine_vac_age_0_from_2000_to_2016 = rep(0.99, length(2000:2016))*S_Korean_pop_gen[1,paste("X", 2000:2016, sep = "")]
S_Korean_rou_year_boundary = 0
for(age in 0:(nrow(S_Korean_pop_gen)-1)){
  temp_S_Korean_National_routine_vac_pop = matrix(0, ncol = ncol(S_Korean_pop_gen), nrow = nrow(S_Korean_pop_gen), dimnames = list(NULL,colnames(S_Korean_pop_gen)))
  if(age == 0){
    temp_S_Korean_National_routine_vac_pop[1,match("X2000",colnames(S_Korean_pop_gen)):match("X2016",colnames(S_Korean_pop_gen))] = unlist(S_Korean_routine_vac_age_0_from_2000_to_2016)
  } else {
    temp_S_Korean_National_routine_vac_pop[(1 + age),head(match("X2000",colnames(S_Korean_pop_gen)):match("X2016",colnames(S_Korean_pop_gen)) + age, - S_Korean_rou_year_boundary)] = head(unlist(S_Korean_routine_vac_age_0_from_2000_to_2016), -S_Korean_rou_year_boundary)
  }
  
  S_Korean_rou_year_boundary = S_Korean_rou_year_boundary + 1
  S_Korean_National_routine_vac_pop = S_Korean_National_routine_vac_pop + temp_S_Korean_National_routine_vac_pop
}
S_Korean_pop_gen = S_Korean_pop_gen - S_Korean_National_routine_vac_pop
S_Korean_pop_gen[S_Korean_pop_gen < 0] = 0

#get age group data from incidence data:
data_inci_Kor = read.xlsx(paste(data_inci_path, "South_Korea_JE.xlsx", sep = ""), 1)

#First study: cohort study:
Kor_subnation_index_1 = c(data_inci_Kor$Reference == unique(data_inci_Kor$Reference)[1])
Kor_study_1 = data_inci_Kor[Kor_subnation_index_1,]
Kor_pop_all_age_year_sum_1 = c()
for(i in 1:14){
  Kor_subnation_index_1_i = c(Kor_study_1$subnation == unique(Kor_study_1$subnation)[i])
  Kor_study_1_i = Kor_study_1[Kor_subnation_index_1_i,]
  Kor_pop_all_age_year_sum_1_temp = pop_age_year_sum_gen(Kor_study_1_i, S_Korean_pop_gen)
  Kor_pop_all_age_year_sum_1 = c(Kor_pop_all_age_year_sum_1, Kor_pop_all_age_year_sum_1_temp)
}

#Second study: 3 hospitals in Seoul
Seoul_hos_pop_2016 = 25514000 #https://en.wikipedia.org/wiki/Seoul_Capital_Area
Seoul_pop_gen = Seoul_hos_pop_2016/sum(S_Korean_pop$X2015)*S_Korean_pop_gen
Kor_subnation_index_2 = c(data_inci_Kor$Reference == unique(data_inci_Kor$Reference)[2])
Kor_study_2 = data_inci_Kor[Kor_subnation_index_2,]
Kor_pop_all_age_year_sum_2 = pop_age_year_sum_gen(Kor_study_2, Seoul_pop_gen)

#Third study: national
Kor_subnation_index_3 = c(data_inci_Kor$Reference == unique(data_inci_Kor$Reference)[3])
Kor_study_3 = data_inci_Kor[Kor_subnation_index_3,]
Kor_pop_all_age_year_sum_3 = pop_age_year_sum_gen(Kor_study_3, S_Korean_pop_gen)

#create the Pop_all_age_year_sum column:
data_inci_Kor$Pop_all_age_year_sum = c(Kor_pop_all_age_year_sum_1, Kor_pop_all_age_year_sum_2)
write.xlsx(data_inci_Kor, paste(data_inci_path, "South_Korea_JE.xlsx", sep = ""), row.names = F)

##################################
#Nepal:
#Nepal census data:
Nepal_pop = naive_pop[naive_pop$country == "NPL",]
Nepal_pop = Nepal_pop[,paste("X", year_sel, sep = "")]
#get age group data from incidence data:
data_inci_Nepal = read.xlsx(paste(data_inci_path, "Nepal_JE.xlsx", sep = ""), 1)

#first study:
Nepal_subnation_index_1 = c(data_inci_Nepal$subnation == unique(data_inci_Nepal$subnation)[1])
Nepal_study_1 = data_inci_Nepal[Nepal_subnation_index_1,]
#Kath population:
Kathmandu_valley_pop_2005 = 793000 #data from: http://worldpopulationreview.com/world-cities/kathmandu-population/
Nepal_Kath_pop_gen = sub_pop_gen(Kathmandu_valley_pop_2005, "X2005", Nepal_pop)
#pop age gen
Nepal_pop_all_age_year_sum_1 = pop_age_year_sum_gen(Nepal_study_1, Nepal_Kath_pop_gen)

#Second study:incomplete
Nepal_subnation_index_2 = c(data_inci_Nepal$subnation == unique(data_inci_Nepal$subnation)[2])
Nepal_study_2 = data_inci_Nepal[Nepal_subnation_index_2,]
#non_Kath population:
non_Kathmandu_valley_pop_2001 = 1070000 - Kathmandu_valley_pop_2005 #assumption from Campbell 2011:
Nepal_non_Kath_pop_gen = sub_pop_gen(non_Kathmandu_valley_pop_2001, "X2001", Nepal_pop)
#pop age gen
Nepal_pop_all_age_year_sum_2 = pop_age_year_sum_gen(Nepal_study_2, Nepal_non_Kath_pop_gen)

#Third and forth studies: pop data included in the paper
Nepal_subnation_index_3 = c(data_inci_Nepal$subnation == unique(data_inci_Nepal$subnation)[3])
Nepal_study_3 = data_inci_Nepal[Nepal_subnation_index_3,]
Nepal_W_Terai_pop = sum(Nepal_study_3$Pop_all_age)/sum(Nepal_pop[1:45,"X2005"])*Nepal_pop
#add vaccination information from the paper:
#vaccinate 224000/15 for each of age group from 1 to 15 in 1999 and 
#and vaccinate 37811.2/2 for each of age group from 1 to 10 in 2000 and 2001
Nepal_W_Terai_vac_pop = matrix(0, ncol = ncol(Nepal_W_Terai_pop), nrow = nrow(Nepal_W_Terai_pop), dimnames = list(NULL,colnames(Nepal_W_Terai_pop)))
for(age in 0:(nrow(Nepal_W_Terai_pop)-1)){
  temp_Nepal_W_Terai_vac_pop = matrix(0, ncol = ncol(Nepal_W_Terai_pop), nrow = nrow(Nepal_W_Terai_pop), dimnames = list(NULL,colnames(Nepal_W_Terai_pop)))
  temp_Nepal_W_Terai_vac_pop[(1:15 + age),(match("X1999",colnames(temp_Nepal_W_Terai_vac_pop)) + age)] =  224000/15 #vaccination in 1999
  temp_Nepal_W_Terai_vac_pop[(1:10 + age),(match(c("X2000","X2001"),colnames(temp_Nepal_W_Terai_vac_pop)) + age)] = 37811.2/2 #vaccination in 2000, 2001
  temp_Nepal_W_Terai_vac_pop = temp_Nepal_W_Terai_vac_pop[,1:ncol(Nepal_W_Terai_vac_pop)]
  Nepal_W_Terai_vac_pop = Nepal_W_Terai_vac_pop + temp_Nepal_W_Terai_vac_pop
}
Nepal_W_Terai_pop = Nepal_W_Terai_pop - Nepal_W_Terai_vac_pop
Nepal_W_Terai_pop[Nepal_W_Terai_pop < 0] = 0

Nepal_pop_all_age_year_sum_3 = pop_age_year_sum_gen(Nepal_study_3, Nepal_W_Terai_pop)

Nepal_subnation_index_4 = c(data_inci_Nepal$subnation == unique(data_inci_Nepal$subnation)[4])
Nepal_study_4 = data_inci_Nepal[Nepal_subnation_index_4,] 
Nepal_pop_all_age_year_sum_4 = Nepal_study_4$Pop_all_age*2

#Fifth Study: => vaccination information from GAVI file: 
#previous vaccination from terai:
Nepal_pop_all_vac = Nepal_pop - Nepal_W_Terai_vac_pop
#campaign vaccination in 2005, 2006, 2008, 2016
#Routine vaccination from 2009 to 2017
#campaign vaccination
Nepal_National_campaign_vac_pop = matrix(0, ncol = ncol(Nepal_pop_all_vac), nrow = nrow(Nepal_pop_all_vac), dimnames = list(NULL,colnames(Nepal_pop_all_vac)))
for(age in 0:(nrow(Nepal_pop_all_vac)-1)){
  temp_Nepal_National_campaign_vac_pop = matrix(0, ncol = ncol(Nepal_pop_all_vac), nrow = nrow(Nepal_pop_all_vac), dimnames = list(NULL,colnames(Nepal_pop_all_vac)))
  temp_Nepal_National_campaign_vac_pop[(1:15 + age),(match("X2005",colnames(temp_Nepal_W_Terai_vac_pop)) + age)] = 745099*0.817538340542666/15 #vaccination in 2005
  temp_Nepal_National_campaign_vac_pop[(1:15 + age),(match("X2006",colnames(temp_Nepal_W_Terai_vac_pop)) + age)] = 3319229*0.764243744556341/15 #vaccination in 2006
  temp_Nepal_National_campaign_vac_pop[(1:15 + age),(match("X2008",colnames(temp_Nepal_W_Terai_vac_pop)) + age)] = 1967051*0.96225415609458/15 #vaccination in 2008
  temp_Nepal_National_campaign_vac_pop[(1:15 + age),(match("X2016",colnames(temp_Nepal_W_Terai_vac_pop)) + age)] = 3251435*1.058/15 #vaccination in 2016
  temp_Nepal_National_campaign_vac_pop = temp_Nepal_National_campaign_vac_pop[,1:ncol(Nepal_National_campaign_vac_pop)]
  Nepal_National_campaign_vac_pop = Nepal_National_campaign_vac_pop + temp_Nepal_National_campaign_vac_pop
}
Nepal_pop_all_vac = Nepal_pop_all_vac - Nepal_National_campaign_vac_pop
Nepal_pop_all_vac[Nepal_pop_all_vac < 0] = 0
#routine vaccination:
Nepal_National_routine_vac_pop = matrix(0, ncol = ncol(Nepal_pop_all_vac), nrow = nrow(Nepal_pop_all_vac), dimnames = list(NULL,colnames(Nepal_pop_all_vac)))
routine_vac_age_0_from_2009_to_2017 = c(0.21, 0.47, 0.54, 0.62, 0.72, 0.78, 0.697407, 0.71, 0.87)*Nepal_pop[1,63:71]
for(age in 0:(nrow(Nepal_pop_all_vac)-1)){
  temp_Nepal_National_routine_vac_pop = matrix(0, ncol = ncol(Nepal_pop_all_vac), nrow = nrow(Nepal_pop_all_vac), dimnames = list(NULL,colnames(Nepal_pop_all_vac)))
  temp_Nepal_National_routine_vac_pop[(1 + age),(match("X2009",colnames(temp_Nepal_W_Terai_vac_pop)):match("X2017",colnames(temp_Nepal_W_Terai_vac_pop)) + age)] = unlist(routine_vac_age_0_from_2009_to_2017)
  temp_Nepal_National_routine_vac_pop = temp_Nepal_National_routine_vac_pop[,1:ncol(Nepal_National_routine_vac_pop)]
  Nepal_National_routine_vac_pop = Nepal_National_routine_vac_pop + temp_Nepal_National_routine_vac_pop
}
Nepal_pop_all_vac = Nepal_pop_all_vac - Nepal_National_routine_vac_pop
Nepal_pop_all_vac[Nepal_pop_all_vac < 0] = 0

Nepal_subnation_index_5 = c(data_inci_Nepal$subnation == unique(data_inci_Nepal$subnation)[5])
Nepal_study_5 = data_inci_Nepal[Nepal_subnation_index_5,]
Nepal_pop_all_age_year_sum_5 = pop_age_year_sum_gen(Nepal_study_5, Nepal_pop_all_vac)

#Sixth Study: => non-western regions, no information of vaccination
Nepal_subnation_index_6 = c(data_inci_Nepal$subnation == unique(data_inci_Nepal$subnation)[6])
Nepal_study_6 = data_inci_Nepal[Nepal_subnation_index_6,]
#pop from 6 districts: Chitwan, Makwanpur, Bara, Parsa, Nawalparasi, Rupandehi
Nepal_6district_2011 = 579984 + 420477 + 687708 + 601017 + 643508 + 880196
Nepal_6district_pop_gen = Nepal_pop_all_vac*(Nepal_6district_2011/sum(Nepal_pop[,"X2011"]))
#data from: https://data.humdata.org/dataset/nepal-census-2011-district-profiles-demography
Nepal_pop_all_age_year_sum_6 = pop_age_year_sum_gen(Nepal_study_6, Nepal_6district_pop_gen)

#Seventh Study: => vaccination information:
Nepal_subnation_index_7 = c(data_inci_Nepal$subnation == unique(data_inci_Nepal$subnation)[7])
Nepal_study_7 = data_inci_Nepal[Nepal_subnation_index_7,]
Nepal_pop_all_age_year_sum_7 = pop_age_year_sum_gen(Nepal_study_7, Nepal_pop_all_vac)

#Eighth Study: (combine of 1st and 2nd study)
Nepal_pop_all_age_year_sum_8 = Nepal_pop_all_age_year_sum_1 + Nepal_pop_all_age_year_sum_2

#Nineth Study: (combine of 3rd and 4th study)
Nepal_pop_all_age_year_sum_9 = Nepal_pop_all_age_year_sum_3 + Nepal_pop_all_age_year_sum_4

#Tenth study: South and eastern regions of Nepal:
BPKIHS_south_eastern_hos_2005 = 10383478 #based on non-Terai region on the paper in study 3rd and 4th
Nepal_south_eastern_pop_gen = sub_pop_gen(BPKIHS_south_eastern_hos_2005, "X2005", Nepal_pop)
Nepal_subnation_index_10 = c(data_inci_Nepal$Reference == unique(data_inci_Nepal$Reference)[9])
Nepal_study_10 = data_inci_Nepal[Nepal_subnation_index_10,]
Nepal_pop_all_age_year_sum_10 = pop_age_year_sum_gen(Nepal_study_10, Nepal_south_eastern_pop_gen)

#Eleventh study: 
Nepal_subnation_index_11 = c(data_inci_Nepal$Reference == unique(data_inci_Nepal$Reference)[10])
Nepal_study_11 = data_inci_Nepal[Nepal_subnation_index_11,]
#cohort study:
Nepal_pop_all_age_year_sum_11 = c()
for(i in 1:5){
  Nepal_subnation_index_11_cohort = c(Nepal_study_11$subnation == unique(Nepal_study_11$subnation)[i])
  Nepal_study_11_cohort = Nepal_study_11[Nepal_subnation_index_11_cohort,]
  Nepal_pop_all_age_year_sum_11_temp = pop_age_year_sum_gen(Nepal_study_11_cohort, Nepal_pop_all_vac)
  Nepal_pop_all_age_year_sum_11 = c(Nepal_pop_all_age_year_sum_11, Nepal_pop_all_age_year_sum_11_temp)
}

#Twefth study: Kosi zone
Kosi_zone_2011 = 2335047 #https://en.wikipedia.org/wiki/Kosi_Zone
Nepal_Kosi_zone_pop_gen = sub_pop_gen(Kosi_zone_2011, "X2011", Nepal_pop)
Nepal_subnation_index_12 = c(data_inci_Nepal$Reference == unique(data_inci_Nepal$Reference)[11])
Nepal_study_12 = data_inci_Nepal[Nepal_subnation_index_12,]
Nepal_pop_all_age_year_sum_12 = pop_age_year_sum_gen(Nepal_study_12, Nepal_Kosi_zone_pop_gen)

#Thirdteenth study: Banke, Bardia and Kailali districts
dists_West_Terai_2011 = 491313 +  426576 + 775709 #https://en.wikipedia.org/wiki/Banke_District ; https://en.wikipedia.org/wiki/Bardiya_District ; 
Nepal_dists_West_Terai_pop_gen = sub_pop_gen(dists_West_Terai_2011, "X2011", Nepal_pop)
Nepal_subnation_index_13 = c(data_inci_Nepal$Reference == unique(data_inci_Nepal$Reference)[12])
Nepal_study_13 = data_inci_Nepal[Nepal_subnation_index_13,]
Nepal_pop_all_age_year_sum_13 = pop_age_year_sum_gen(Nepal_study_13, Nepal_dists_West_Terai_pop_gen)

#create the Pop_all_age_year_sum column:
data_inci_Nepal$Pop_all_age_year_sum = c(Nepal_pop_all_age_year_sum_1, Nepal_pop_all_age_year_sum_2, Nepal_pop_all_age_year_sum_3, Nepal_pop_all_age_year_sum_4, 
                                         Nepal_pop_all_age_year_sum_5, Nepal_pop_all_age_year_sum_6, Nepal_pop_all_age_year_sum_7, Nepal_pop_all_age_year_sum_8, 
                                         Nepal_pop_all_age_year_sum_9, Nepal_pop_all_age_year_sum_10, Nepal_pop_all_age_year_sum_11, Nepal_pop_all_age_year_sum_12, Nepal_pop_all_age_year_sum_13)
write.xlsx(data_inci_Nepal, paste(data_inci_path, "Nepal_JE.xlsx", sep = ""), row.names = F)

##################################
#Taiwan:
#Taiwan census data:
Taiwan_pop = convert_UN_data[convert_UN_data$country == "TWN", ]
Taiwan_pop = Taiwan_pop[, paste("X",1951:2015, sep = "")]

#Vaccination program in Taiwan: from 1968 to 1979 - campaign vaccination for all < 3 years old, from 1980 onward - routine vaccination for older 15m.
#campaign vaccination
Taiwan_National_campaign_vac_pop = data.frame(matrix(0, ncol = ncol(Taiwan_pop), nrow = nrow(Taiwan_pop), dimnames = list(NULL,colnames(Taiwan_pop))))
Taiwan_cam_year_boundary = 0
for(age in 0:(nrow(Taiwan_pop)-1)){
  temp_Taiwan_National_campaign_vac_pop = data.frame(matrix(0, ncol = ncol(Taiwan_pop), nrow = nrow(Taiwan_pop), dimnames = list(NULL,colnames(Taiwan_pop))))
  year_index_vac = match(paste("X", 1968:1979, sep = ""),colnames(Taiwan_pop)) + age
  if(age<37){
    temp_Taiwan_National_campaign_vac_pop[(1:3 + age),year_index_vac] = Taiwan_pop[1:3,year_index_vac - age]
  } else {
    Taiwan_cam_year_boundary = Taiwan_cam_year_boundary + 1
    temp_Taiwan_National_campaign_vac_pop[(1:3 + age),head(year_index_vac, - Taiwan_cam_year_boundary)] = Taiwan_pop[1:3,head(year_index_vac - age, - Taiwan_cam_year_boundary)]
  }
  temp_Taiwan_National_campaign_vac_pop = temp_Taiwan_National_campaign_vac_pop[,1:ncol(Taiwan_National_campaign_vac_pop)]
  Taiwan_National_campaign_vac_pop = Taiwan_National_campaign_vac_pop + temp_Taiwan_National_campaign_vac_pop
}
Taiwan_pop = Taiwan_pop - Taiwan_National_campaign_vac_pop
Taiwan_pop[Taiwan_pop < 0] = 0
#routine vaccination:
Taiwan_National_routine_vac_pop = matrix(0, ncol = ncol(Taiwan_pop), nrow = nrow(Taiwan_pop), dimnames = list(NULL,colnames(Taiwan_pop)))
Taiwan_routine_vac_age_0_from_1980_to_2015 = rep(0.99, length(1980:2015))*Taiwan_pop[1,paste("X", 1980:2015, sep = "")]
Taiwan_rou_year_boundary = 0
for(age in 0:(nrow(Taiwan_pop)-1)){
  temp_Taiwan_National_routine_vac_pop = matrix(0, ncol = ncol(Taiwan_pop), nrow = nrow(Taiwan_pop), dimnames = list(NULL,colnames(Taiwan_pop)))
  if(age == 0){
    temp_Taiwan_National_routine_vac_pop[1,match("X1980",colnames(Taiwan_pop)):match("X2015",colnames(Taiwan_pop))] = unlist(Taiwan_routine_vac_age_0_from_1980_to_2015)
  } else {
    temp_Taiwan_National_routine_vac_pop[(1 + age),head(match("X1980",colnames(Taiwan_pop)):match("X2015",colnames(Taiwan_pop)) + age, - Taiwan_rou_year_boundary)] = head(unlist(Taiwan_routine_vac_age_0_from_1980_to_2015), -Taiwan_rou_year_boundary)
  }
  
  Taiwan_rou_year_boundary = Taiwan_rou_year_boundary + 1
  Taiwan_National_routine_vac_pop = Taiwan_National_routine_vac_pop + temp_Taiwan_National_routine_vac_pop
}
Taiwan_pop = Taiwan_pop - Taiwan_National_routine_vac_pop
Taiwan_pop[Taiwan_pop < 0] = 0

#get age group data from incidence data:
data_inci_Taiwan = read.xlsx(paste(data_inci_path, "Taiwan_JE.xlsx", sep = ""), 1)

#cohort study:
Taiwan_pop_all_age_year_sum_1 = c()
for(i in 1:14){
  Taiwan_subnation_index = c(data_inci_Taiwan$subnation == unique(data_inci_Taiwan$subnation)[i])
  Taiwan_study = data_inci_Taiwan[Taiwan_subnation_index,]
  Taiwan_pop_all_age_year_sum_1_temp = pop_age_year_sum_gen(Taiwan_study, Taiwan_pop)
  Taiwan_pop_all_age_year_sum_1 = c(Taiwan_pop_all_age_year_sum_1, Taiwan_pop_all_age_year_sum_1_temp)
}

#create the Pop_all_age_year_sum column:
data_inci_Taiwan$Pop_all_age_year_sum = Taiwan_pop_all_age_year_sum_1
write.xlsx(data_inci_Taiwan, paste(data_inci_path, "Taiwan_JE.xlsx", sep = ""), row.names = F)

###################################
#Vietnam:
#Vietnam census data:
Vietnam_pop = naive_pop[naive_pop$country == "VNM", ]
Vietnam_pop_gen = Vietnam_pop[, paste("X",1950:2044, sep = "")]

#get age group data from incidence data:
data_inci_Vietnam = read.xlsx(paste(data_inci_path, "Vietnam_JE.xlsx", sep = ""), 1)

#First study:
#population of 12 districts: Hue, KH, DL, B Phuoc, D Thap, A Giang, K Giang, Can Tho, Soc Trang, Bac Lieu, Ca Mau, Tra Vinh
Vietnam_12dist_pop_2011 = (1103.1 + 1174.1 + 1771.8 + 905.3 + 1673.2 + 2151 + 1200.3 + 1303.7 + 873.3 +  1214.9)*1000 #data from: http://www.gso.gov.vn/default_en.aspx?tabid=467&ItemID=12941
Vietnam_12dist_pop_gen = sub_pop_gen(Vietnam_12dist_pop_2011, "X2011", Vietnam_pop_gen)
#get age group
Vietnam_subnation_index_1 = c(data_inci_Vietnam$subnation == unique(data_inci_Vietnam$subnation)[1])
Vietnam_study_1 = data_inci_Vietnam[Vietnam_subnation_index_1,]
#get pop:
Vietnam_pop_all_gen_1 = pop_age_year_sum_gen(Vietnam_study_1, Vietnam_12dist_pop_gen)

#Second study:
#get population from 5 provinces:Bac giang, Hai Duong, Hai Phong, Thai Binh, Thanh Hoa.
Vietnam_5provinces_pop_2011 = (1574.3 + 1718.9 + 1878.5 + 1786 + 3412.6)*1000
Vietnam_5provinces_pop_gen = sub_pop_gen(Vietnam_5provinces_pop_2011, "X2011", Vietnam_pop_gen)
#get age group
Vietnam_subnation_index_2 = c(data_inci_Vietnam$subnation == unique(data_inci_Vietnam$subnation)[2])
Vietnam_study_2 = data_inci_Vietnam[Vietnam_subnation_index_2,]
#gen pop:
Vietnam_pop_all_gen_2 = pop_age_year_sum_gen(Vietnam_study_2, Vietnam_5provinces_pop_gen)

#create the Pop_all_age_year_sum column:
data_inci_Vietnam$Pop_all_age_year_sum = c(Vietnam_pop_all_gen_1, Vietnam_pop_all_gen_2)
write.xlsx(data_inci_Vietnam, paste(data_inci_path, "Vietnam_JE.xlsx", sep = ""), row.names = F)

###################################
#Thailand:
#Thailand census data:
Thailand_pop = IDM_data_reformat[IDM_data_reformat$country == "THA", ]
Thailand_pop_gen = Thailand_pop[, paste("X",1981:2025, sep = "")]
#population of 2 districts: Bangkok, Haiyan
Thailand_2dist_pop_2011 = 8280925 + 158218 #data from: https://en.wikipedia.org/wiki/Provinces_of_Thailand
Thailand_2dist_pop_gen = sub_pop_gen(Thailand_2dist_pop_2011, "X2011", Thailand_pop_gen)

Thailand_2dist_routine_vac_pop = matrix(0, ncol = ncol(Thailand_2dist_pop_gen), nrow = nrow(Thailand_2dist_pop_gen), dimnames = list(NULL,colnames(Thailand_2dist_pop_gen)))
Thailand_routine_vac_age_0_from_2000_to_2004 = c(0.99, 0.99, 0.99, 0.99)*Thailand_2dist_pop_gen[1,20:24]
for(age in 0:(nrow(Thailand_2dist_pop_gen)-1)){
  temp_Thailand_2dist_routine_vac_pop = matrix(0, ncol = ncol(Thailand_2dist_pop_gen), nrow = nrow(Thailand_2dist_pop_gen), dimnames = list(NULL,colnames(Thailand_2dist_pop_gen)))
  temp_Thailand_2dist_routine_vac_pop[(1 + age),(match("X2000",colnames(temp_Thailand_2dist_routine_vac_pop)):match("X2004",colnames(temp_Thailand_2dist_routine_vac_pop)) + age)] = unlist(Thailand_routine_vac_age_0_from_2000_to_2004)
  temp_Thailand_2dist_routine_vac_pop = temp_Thailand_2dist_routine_vac_pop[,1:ncol(Thailand_2dist_routine_vac_pop)]
  Thailand_2dist_routine_vac_pop = Thailand_2dist_routine_vac_pop + temp_Thailand_2dist_routine_vac_pop
}
Thailand_2dist_pop_gen = Thailand_2dist_pop_gen - Thailand_2dist_routine_vac_pop
Thailand_2dist_pop_gen[Thailand_2dist_pop_gen < 0] = 0


#get age group data from incidence data:
data_inci_Thailand = read.xlsx(paste(data_inci_path, "Thailand_JE.xlsx", sep = ""), 1)

Thailand_subnation_index_1 = c(data_inci_Thailand$subnation == unique(data_inci_Thailand$subnation)[1])
Thailand_study_1 = data_inci_Thailand[Thailand_subnation_index_1,]

#create the Pop_all_age_year_sum column:
data_inci_Thailand$Pop_all_age_year_sum = pop_age_year_sum_gen(Thailand_study_1, Thailand_2dist_pop_gen)
write.xlsx(data_inci_Thailand, paste(data_inci_path, "Thailand_JE.xlsx", sep = ""), row.names = F)

###################################
#Bangladesh:
#Bangladesh census data:
Bangladesh_pop = naive_pop[naive_pop$country == "BGD", ]
Bangladesh_pop_gen = Bangladesh_pop[, paste("X",1950:2044, sep = "")]

#get age group data from incidence data:
data_inci_Bangladesh = read.xlsx(paste(data_inci_path, "Bangladesh_JE.xlsx", sep = ""), 1)

#first study:
#population of 4 divisions: Dhaka, Mymensingh, Rajshahi, Sylhet
Bangladesh_4div_pop_2011 = 46729000 + 11370000 + 18484858 + 9910219 #data from: https://en.wikipedia.org/wiki/Divisions_of_Bangladesh
Bangladesh_4div_pop_gen = sub_pop_gen(Bangladesh_4div_pop_2011, "X2011", Bangladesh_pop_gen)
#pop all age gen:
Bangladesh_subnation_index_1 = c(data_inci_Bangladesh$subnation == unique(data_inci_Bangladesh$subnation)[1])
Bangladesh_study_1 = data_inci_Bangladesh[Bangladesh_subnation_index_1,]
Bangladesh_pop_all_age_year_sum_1 = pop_age_year_sum_gen(Bangladesh_study_1, Bangladesh_4div_pop_gen)

#second study:
#catchment areas of Rajshahi Medical College Hospital: Rajshahi, Naogaon, Chapai-Nawabgonj districts
Bangladesh_3hos_pop_2011 = 46729000 + 11370000 + 18484858 + 9910219 #data from: https://en.wikipedia.org/wiki/Divisions_of_Bangladesh
Bangladesh_3hos_pop_gen = sub_pop_gen(Bangladesh_3hos_pop_2011, "X2011", Bangladesh_pop_gen)
#pop all age gen:
Bangladesh_subnation_index_2 = c(data_inci_Bangladesh$subnation == unique(data_inci_Bangladesh$subnation)[2])
Bangladesh_study_2 = data_inci_Bangladesh[Bangladesh_subnation_index_2,]
Bangladesh_pop_all_age_year_sum_2 = pop_age_year_sum_gen(Bangladesh_study_2, Bangladesh_3hos_pop_gen)

#create the Pop_all_age_year_sum column:
data_inci_Bangladesh$Pop_all_age_year_sum
write.xlsx(data_inci_Bangladesh, paste(data_inci_path, "Bangladesh_JE.xlsx", sep = ""), row.names = F)

###################################
#China:
#China census data:
China_pop = naive_pop[naive_pop$country == "CHN", paste("X",year_sel, sep = "")]
China_pop_gen = China_pop
#Routine vaccination in 2008:
China_National_routine_vac_pop = matrix(0, ncol = ncol(China_pop_gen), nrow = nrow(China_pop_gen), dimnames = list(NULL,colnames(China_pop_gen)))
China_routine_vac_age_0_from_2008_to_2015 = rep(0.99, length(2008:2015))*China_pop_gen[1,paste("X", 2008:2015, sep = "")]
China_rou_year_boundary = 0
for(age in 0:(nrow(China_pop_gen)-1)){
  temp_China_National_routine_vac_pop = matrix(0, ncol = ncol(China_pop_gen), nrow = nrow(China_pop_gen), dimnames = list(NULL,colnames(China_pop_gen)))
  if(age == 0){
    temp_China_National_routine_vac_pop[1,match("X2008",colnames(China_pop_gen)):match("X2015",colnames(China_pop_gen))] = unlist(China_routine_vac_age_0_from_2008_to_2015)
  } else {
    temp_China_National_routine_vac_pop[(1 + age),head(match("X2008",colnames(China_pop_gen)):match("X2015",colnames(China_pop_gen)) + age, - China_rou_year_boundary)] = head(unlist(China_routine_vac_age_0_from_2008_to_2015), -China_rou_year_boundary)
  }
  
  China_rou_year_boundary = China_rou_year_boundary + 1
  China_National_routine_vac_pop = China_National_routine_vac_pop + temp_China_National_routine_vac_pop
}
China_pop_gen = China_pop_gen - China_National_routine_vac_pop
China_pop_gen[China_pop_gen < 0] = 0

#get age group data from incidence data:
data_inci_China = read.xlsx(paste(data_inci_path, "China_JE.xlsx", sep = ""), 1)

#First study: Guizhou cohort: 2004-2009
#poplation of Guizhou
China_Guizhou_pop_2010 = 	34746468 #https://en.wikipedia.org/wiki/Guizhou
China_Guizhou_pop_gen = China_pop_gen*China_Guizhou_pop_2010/sum(China_pop$X2010)
#campaign vaccination in Guizhou from 2004 to 2008:
China_Guizhou_campaign_vac_pop = matrix(0, ncol = ncol(China_Guizhou_pop_gen), nrow = nrow(China_Guizhou_pop_gen), dimnames = list(NULL,colnames(China_Guizhou_pop_gen)))
for(age in 0:(nrow(China_Guizhou_pop_gen)-1)){
  temp_China_Guizhou_campaign_vac_pop = matrix(0, ncol = ncol(China_Guizhou_pop_gen), nrow = nrow(China_Guizhou_pop_gen), dimnames = list(NULL,colnames(China_Guizhou_pop_gen)))
  temp_China_Guizhou_campaign_vac_pop[(1:11 + age),(match("X2004",colnames(temp_China_Guizhou_campaign_vac_pop)) + age)] = 351045/11 #vaccination in 2004
  temp_China_Guizhou_campaign_vac_pop[(1:11+ age),(match("X2005",colnames(temp_China_Guizhou_campaign_vac_pop)) + age)] = 396968/11 #vaccination in 2005
  temp_China_Guizhou_campaign_vac_pop[(1:11 + age),(match("X2006",colnames(temp_China_Guizhou_campaign_vac_pop)) + age)] = 754608/11 #vaccination in 2006
  temp_China_Guizhou_campaign_vac_pop[(1:7 + age),(match("X2007",colnames(temp_China_Guizhou_campaign_vac_pop)) + age)] = 2201124/7 #vaccination in 2007
  temp_China_Guizhou_campaign_vac_pop[(1:11 + age),(match("X2008",colnames(temp_China_Guizhou_campaign_vac_pop)) + age)] = 1503052/11 #vaccination in 2008
  temp_China_Guizhou_campaign_vac_pop = temp_China_Guizhou_campaign_vac_pop[,1:ncol(China_Guizhou_campaign_vac_pop)]
  China_Guizhou_campaign_vac_pop = China_Guizhou_campaign_vac_pop + temp_China_Guizhou_campaign_vac_pop
}
China_Guizhou_pop_gen = China_Guizhou_pop_gen - China_Guizhou_campaign_vac_pop
China_Guizhou_pop_gen[China_Guizhou_pop_gen < 0] = 0
#pop all age gen:
China_subnation_index_1 = c(data_inci_China$Reference == unique(data_inci_China$Reference)[1])
China_study_1 = data_inci_China[China_subnation_index_1,]
#cohort study:
China_pop_all_age_year_sum_1 = c()
for(i in 1:6){
  China_subnation_index_1_cohort = c(China_study_1$subnation == unique(China_study_1$subnation)[i])
  China_study_1_cohort = China_study_1[China_subnation_index_1_cohort,]
  China_pop_all_age_year_sum_1_temp = pop_age_year_sum_gen(China_study_1_cohort, China_Guizhou_pop_gen)
  China_pop_all_age_year_sum_1 = c(China_pop_all_age_year_sum_1, China_pop_all_age_year_sum_1_temp)
}
data_inci_China[China_subnation_index_1,]$Pop_all_age_year_sum = China_pop_all_age_year_sum_1

#second study: Guizhou 2006
#pop all age gen:
China_subnation_index_2 = c(data_inci_China$Reference == unique(data_inci_China$Reference)[2])
China_study_2 = data_inci_China[China_subnation_index_2,]
China_pop_all_age_year_sum_2 = pop_age_year_sum_gen(China_study_2, China_Guizhou_pop_gen)
data_inci_China[China_subnation_index_2,]$Pop_all_age_year_sum = China_pop_all_age_year_sum_2

#Third study: 4 cities: Jinan, Yichang, Shijiazhuang, Guigang
#poplation of Jinan, Yichang, Shijiazhuang, Guigang
China_4cities_pop = list(pop = c(6814000, 4059686, 10701600, 4400000), 
                         year = c("X2010","X2010","X2015", "X2004"))
#https://en.wikipedia.org/wiki/Jinan https://en.wikipedia.org/wiki/Yichanghttps://en.wikipedia.org/wiki/Shijiazhuang  https://en.wikipedia.org/wiki/Guigang
#pop all age gen:
China_subnation_index_3 = c(data_inci_China$Reference == unique(data_inci_China$Reference)[3])
China_study_3 = data_inci_China[China_subnation_index_3,]
#cohort study:
China_pop_all_age_year_sum_3 = c()
for(i in 1:4){
  China_subnation_index_3_4cities = c(China_study_3$subnation == unique(China_study_3$subnation)[i])
  China_study_3_4cities = China_study_3[China_subnation_index_3_4cities,]
  China_4cities_pop_gen = China_pop_gen*China_4cities_pop$pop[i]/sum(China_pop[,China_4cities_pop$year[i]])
  
  China_pop_all_age_year_sum_3_temp = pop_age_year_sum_gen(China_study_3_4cities, China_4cities_pop_gen)
  China_pop_all_age_year_sum_3 = c(China_pop_all_age_year_sum_3, China_pop_all_age_year_sum_3_temp)
}
data_inci_China[China_subnation_index_3,]$Pop_all_age_year_sum = China_pop_all_age_year_sum_3

#forth study: cohort data
#Endemic poulation of China based on campbell 2011
China_endemic_national_pop_2010_2007 = (1025.7 + 276.6)*1000000
China_endemic_pop_gen = sub_pop_gen(China_endemic_national_pop_2010_2007,"X2010",China_pop)
#pop all age gen:
China_subnation_index_4 = c(data_inci_China$Reference == unique(data_inci_China$Reference)[4])
China_study_4 = data_inci_China[China_subnation_index_4,]
#cohort study:
China_pop_all_age_year_sum_4 = c()
for(i in 1:10){
  China_subnation_index_4_cohort = c(China_study_4$subnation == unique(China_study_4$subnation)[i])
  China_study_4_cohort = China_study_4[China_subnation_index_4_cohort,]
  China_pop_all_age_year_sum_4_temp = pop_age_year_sum_gen(China_study_4_cohort, China_endemic_pop_gen)
  China_pop_all_age_year_sum_4 = c(China_pop_all_age_year_sum_4, China_pop_all_age_year_sum_4_temp)
}
data_inci_China[China_subnation_index_4,]$Pop_all_age_year_sum = China_pop_all_age_year_sum_4

#Fifth study: cohort data of Longnan city
#Population of Longnan https://en.wikipedia.org/wiki/Longnan
China_Longnan_2010 = 	2567718
China_Longnan_pop_gen = China_pop_gen*China_Longnan_2010/sum(China_pop[,"X2010"])
#Vaccination program from 2004 to 2009
#pop all age gen:
China_subnation_index_5 = c(data_inci_China$Reference == unique(data_inci_China$Reference)[5])
China_study_5 = data_inci_China[China_subnation_index_5,]
#cohort study:
China_pop_all_age_year_sum_5 = c()
for(i in 1:7){
  China_subnation_index_5_cohort = c(China_study_5$subnation == unique(China_study_5$subnation)[i])
  China_study_5_cohort = China_study_5[China_subnation_index_5_cohort,]
  China_pop_all_age_year_sum_5_temp = pop_age_year_sum_gen(China_study_5_cohort, China_Longnan_pop_gen)
  China_pop_all_age_year_sum_5 = c(China_pop_all_age_year_sum_5, China_pop_all_age_year_sum_5_temp)
}
data_inci_China[China_subnation_index_5,]$Pop_all_age_year_sum = China_pop_all_age_year_sum_5

#sixth study: Guigang - Guangxi - continuous data => discretized age in excel file: pone.0144366.s001.xlsx
#poplation of Guiang city: https://en.wikipedia.org/wiki/Guigang
China_Guigang_Guangxi_pop_2005 = 	4400000
China_Guigang_Guangxi_pop_gen = China_pop_gen*China_Guigang_Guangxi_pop_2005/sum(China_pop[,"X2005"])
#pop all age gen:
China_subnation_index_6 = c(data_inci_China$Reference == unique(data_inci_China$Reference)[6])
China_study_6 = data_inci_China[China_subnation_index_6,]
China_pop_all_age_year_sum_6 = pop_age_year_sum_gen(China_study_6, China_Guigang_Guangxi_pop_gen)
data_inci_China[China_subnation_index_6,]$Pop_all_age_year_sum = China_pop_all_age_year_sum_6

#Seventh study: hospital in Baoji-Shaanxi province
#poplation of Guiang city: https://en.wikipedia.org/wiki/Guigang
China_Baoji_Shaanxi_pop_2010 = 	3716731
China_Baoji_Shaanxi_pop_gen = China_pop_gen*China_Baoji_Shaanxi_pop_2010/sum(China_pop[,"X2010"])
#pop all age gen:
China_subnation_index_7 = c(data_inci_China$Reference == unique(data_inci_China$Reference)[7])
China_study_7 = data_inci_China[China_subnation_index_7,]
China_pop_all_age_year_sum_7 = pop_age_year_sum_gen(China_study_7, China_Baoji_Shaanxi_pop_gen)
data_inci_China[China_subnation_index_7,]$Pop_all_age_year_sum = China_pop_all_age_year_sum_7

#Eighth study: hospital in Hebei province
#poplation of Hebei province: https://en.wikipedia.org/wiki/Hebei
China_Hebei_pop_2016 = 	74700500
China_Hebei_pop_gen = China_pop_gen*China_Hebei_pop_2016/sum(China_pop[,"X2016"])
#pop all age gen:
China_subnation_index_8 = c(data_inci_China$Reference == unique(data_inci_China$Reference)[8])
China_study_8 = data_inci_China[China_subnation_index_8,]
China_pop_all_age_year_sum_8 = pop_age_year_sum_gen(China_study_8, China_Hebei_pop_gen)
data_inci_China[China_subnation_index_8,]$Pop_all_age_year_sum = China_pop_all_age_year_sum_8

#Nineth study: Nyingchi District (3 counties Nyingchi, Mainling, Gongbo'gyamda)
#poplation of 3 counties Nyingchi District: https://en.wikipedia.org/wiki/Nyingchi
China_3_counties_Nyingchi_pop_2010 = 	54702 + 29929 + 22834
China_3_counties_Nyingchi_pop_gen = China_pop_gen*China_3_counties_Nyingchi_pop_2010/sum(China_pop[,"X2010"])
#pop all age gen:
China_subnation_index_9 = c(data_inci_China$Reference == unique(data_inci_China$Reference)[9])
China_study_9 = data_inci_China[China_subnation_index_9,]
China_pop_all_age_year_sum_9 = pop_age_year_sum_gen(China_study_9, China_3_counties_Nyingchi_pop_gen)
data_inci_China[China_subnation_index_9,]$Pop_all_age_year_sum = China_pop_all_age_year_sum_9

#Tenth study: Wuhan => too complicated vaccinatino program, only use the info of national routine vac in 2008
#poplation of Wuhan: https://en.wikipedia.org/wiki/Wuhan
China_Wuhan_pop_2015 = 	10607700
China_Wuhan_pop_gen = China_pop_gen*China_Wuhan_pop_2015/sum(China_pop[,"X2015"])
#pop all age gen:
China_subnation_index_10 = c(data_inci_China$Reference == unique(data_inci_China$Reference)[10])
China_study_10 = data_inci_China[China_subnation_index_10,]
China_pop_all_age_year_sum_10 = pop_age_year_sum_gen(China_study_10, China_Wuhan_pop_gen)
data_inci_China[China_subnation_index_10,]$Pop_all_age_year_sum = China_pop_all_age_year_sum_10

#Eleventh study: Mainland China
#poplation of Mainland China: https://en.wikipedia.org/wiki/Demographics_of_China
China_Mainland_China_pop_2016 = 	1403500365
China_Mainland_China_pop_gen = China_pop_gen*China_Mainland_China_pop_2016/sum(China_pop[,"X2016"])
#pop all age gen:
China_subnation_index_11 = c(data_inci_China$Reference == unique(data_inci_China$Reference)[11])
China_study_11 = data_inci_China[China_subnation_index_11,]
#cohort study:
China_pop_all_age_year_sum_11 = c()
for(i in 1:14){
  China_subnation_index_11_cohort = c(China_study_11$subnation == unique(China_study_11$subnation)[i])
  China_study_11_cohort = China_study_11[China_subnation_index_11_cohort,]
  China_pop_all_age_year_sum_11_temp = pop_age_year_sum_gen(China_study_11_cohort, China_Mainland_China_pop_gen)
  China_pop_all_age_year_sum_11 = c(China_pop_all_age_year_sum_11, China_pop_all_age_year_sum_11_temp)
}

data_inci_China[China_subnation_index_11,]$Pop_all_age_year_sum = China_pop_all_age_year_sum_11

write.xlsx(data_inci_China, paste(data_inci_path, "China_JE.xlsx", sep = ""), row.names = F)
###################################
#Philippines:
#Philippines census data:
Philippines_pop = naive_pop[naive_pop$country == "PHL", ]

#get age group data from incidence data:
data_inci_Philippines = read.xlsx(paste(data_inci_path, "Philippines_JE.xlsx", sep = ""), 1)

#first study:
Philippines_subnation_index_1 = c(data_inci_Philippines$subnation == unique(data_inci_Philippines$subnation)[1])
Philippines_study_1 = data_inci_Philippines[Philippines_subnation_index_1,]
Philippines_pop_all_age_year_sum_1 = pop_age_year_sum_gen(Philippines_study_1, Philippines_pop)

#second study:
Phil_manila_pop_2010 = 1652171 # data from: https://vi.wikipedia.org/wiki/Manila
Phil_manila_pop_gen = sub_pop_gen(Phil_manila_pop_2010, "X2010", Philippines_pop)
Philippines_subnation_index_2 = c(data_inci_Philippines$subnation == unique(data_inci_Philippines$subnation)[2])
Philippines_study_2 = data_inci_Philippines[Philippines_subnation_index_2,]
Philippines_pop_all_age_year_sum_2 = pop_age_year_sum_gen(Philippines_study_2, Phil_manila_pop_gen)

#third study: 
Phil_North_central_pop_2010 = 42822878 + 15528346	 # data from: https://en.wikipedia.org/wiki/Island_groups_of_the_Philippines
Phil_North_central_pop_gen = sub_pop_gen(Phil_North_central_pop_2010, "X2010", Philippines_pop)
Philippines_subnation_index_3 = c(data_inci_Philippines$subnation == unique(data_inci_Philippines$subnation)[3])
Philippines_study_3 = data_inci_Philippines[Philippines_subnation_index_3,]
Philippines_pop_all_age_year_sum_3 = pop_age_year_sum_gen(Philippines_study_3, Phil_North_central_pop_gen)

#create the Pop_all_age_year_sum column:
data_inci_Philippines$Pop_all_age_year_sum = c(Philippines_pop_all_age_year_sum_1, Philippines_pop_all_age_year_sum_2, Philippines_pop_all_age_year_sum_3)
write.xlsx(data_inci_Philippines, paste(data_inci_path, "Philippines_JE.xlsx", sep = ""), row.names = F)

