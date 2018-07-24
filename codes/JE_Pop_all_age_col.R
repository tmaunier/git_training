###age group population interpolation => generate the Pop_all_age column in country incidence data:
#library loading:
rm(list = ls())
library(xlsx)

#data loading:
#Pop data load:
output_data_path = 'C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_clean/Montagu_data/test3/'
naive_pop = read.xlsx2(paste(output_data_path,'naive_pop_1950_2100.xlsx', sep = ""), 1, colClasses=NA)
data_inci_path = "D:/OUCRU/Hannah/JEV/JEV_model/data_JE_raw/cases_sero_data/"
year_sel = 2000:2100

#function to generate subnation population over time by the national population data and subnation population of a certain time:
sub_pop_gen <-function(sub_pop_year, year, naive_pop_country){
  naive_pop_country
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
    age_group_cam = 1:15
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
#function to calculate Pop_all_age column:
pop_age_gen <- function(age_group, sub_cen_year){
  #Calculate the age_l and age_u:
  age_group_split = strsplit(as.character(age_group),"-")
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
    
    pop_gen_sel = as.numeric(as.character(sub_cen_year))[sel_col_index]
    
    if(age_stop == 0){
      pop_gen_case = c(pop_gen_case,pop_gen_sel)
    } else {
      pop_gen_case = c(pop_gen_case,sum(pop_gen_sel))
    }
  }
  return(pop_gen_case)
}

#pop data:
IDM_data = read.xlsx2("C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_clean/pop_data/IDB_0_100_2000_2045_24.xlsx", 1, 
                      colClasses = c(rep("character", 3), rep("numeric", 45)))
age_group_split = strsplit(as.character(IDM_data$Age_group),"-")
age_from = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[1])}))
age_to = age_from
IDM_data_reformat = data.frame(IDM_data$ISO, age_from, age_to, IDM_data[,4:48], stringsAsFactors = F)
colnames(IDM_data_reformat)[1] = "country"

#Age group select:
age_group_sel = 0:99
year_sel = 2000:2044
#################################
#India:

#load population data of the country:
choosen_country = as.character(data_inci_India$ISO)[1]
India_pop = naive_pop[naive_pop$country == choosen_country,]

#load age group from incidence data:
data_inci_India = read.xlsx(paste(data_inci_path, "India_JE.xlsx", sep = ""), 1)
India_subnation = unlist(unique(data_inci_India$subnation))

#Gorakhpur (Uttar Pradesh) district: First and second study
#data from: http://www.census2011.co.in/census/district/559-gorakhpur.html => linear interpolate => pop in that year
Go_ce2011 = 4440895
Go_ce2001 = 3769456 
Go_gap = (Go_ce2011 - Go_ce2001)/(2011-2001)
Go_ce2005 = Go_ce2001 + Go_gap*(2005-2001)
Go_ce2009 = Go_ce2001 + Go_gap*(2009-2001)

Age_group_2009 = data_inci_India[data_inci_India$subnation == India_subnation[1],]$Age_group
Age_group_2005 = data_inci_India[data_inci_India$subnation == India_subnation[2],]$Age_group

#Data of vaccination campaign:
Go_vac_cam_2006 = data.frame(matrix(c(1390307, 0.9703), ncol = 1, nrow = 2, dimnames = list(c("target", "coverage"), "X2006")))

Go_cen_all_year = sub_pop_gen(Go_ce2011, "X2011", India_pop)

Go_vac_cam_2006_pop = pop_vac_gen(Go_vac_cam_2006, Go_cen_all_year)

Go_pop_age_group_2009 = pop_age_gen(Age_group_2009, Go_vac_cam_2006_pop$X2009)
Go_pop_age_group_2005 = pop_age_gen(Age_group_2005, Go_vac_cam_2006_pop$X2005)


#Uttar Pradesh: Third study
Utt_2011 = 199812341
Utt_2001 = 166197921
Utt_gap = (Utt_2011-Utt_2001)/(2011-2001)
Utt_2005 = Utt_2001 + Utt_gap*(2005-2001)
#age group from data
Utt_Age_group_2005 = as.character(data_inci_India[data_inci_India$subnation == India_subnation[3],]$Age_group)

Utt_cen_all_year = sub_pop_gen(Utt_2011, "X2011", India_pop)

Utt_pop_age_group_2005 = pop_age_gen(Utt_Age_group_2005, Utt_cen_all_year$X2005)

#Tamil Nadu: Fourth study
TaNa_ce2011 = 72147030
#age group from data
TaNa_Age_group = as.character(data_inci_India[data_inci_India$subnation == subnation[4],]$Age_group)

TaNa_cen_all_year = sub_pop_gen(TaNa_ce2011, "X2011", India_pop)
TaNa_pop_age_group_2008 = pop_age_gen(TaNa_Age_group, TaNa_cen_all_year$X2008)

#Assam: Fifth Study
Assam_2011 = 31205576 #data from http://www.census2011.co.in/census/state/assam.html
#age group from data
Assam_Age_group = as.character(data_inci_India[data_inci_India$subnation == India_subnation[5],]$Age_group)
#generate population stratified by age groups:
Assam_cen_all_year = sub_pop_gen(Assam_2011, "X2011", India_pop)
Assam_pop_age_group_2008 = pop_age_gen(Assam_Age_group, Assam_cen_all_year$X2008)

#7 districts Medical College Hospital Assam: sixth study
# 7 upper districts of Assam: Jorhat, Dibrugarh, Dhemaji, Golaghat, Lakhimpur, Sivasagar, and Tinsukia
#+ Arunachal Pradesh + Nagaland
Assam_7districts_2011 = 1092256 + 1326335 + 686133 + 1066888 + 1042137 + 1151050 + 1327929 #data from https://www.mapsofindia.com/maps/assam/assam-district.htm
Assam_7districts_2011 = Assam_7districts_2011 + 1383727 + 1978502#data from : http://www.census2011.co.in/census/state/
#age group from data
Assam_7districts_Age_group = as.character(data_inci_India[data_inci_India$subnation == India_subnation[6],]$Age_group)
#generate population stratified by age groups:
Assam_7districts_cen_all_year = sub_pop_gen(Assam_7districts_2011, "X2011", India_pop)
Assam_7districts_pop_age_group_2012 = pop_age_gen(Assam_7districts_Age_group, Assam_7districts_cen_all_year$X2012)

##################################
#Japan:
IDM_data = read.xlsx2("C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_clean/pop_data/IDB_0_100_2000_2045_24.xlsx", 1, 
                      colClasses = c(rep("character", 3), rep("numeric", 45)))
age_group_split = strsplit(as.character(IDM_data$Age_group),"-")
age_from = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[1])}))
age_to = age_from
IDM_data_reformat = data.frame(IDM_data$ISO, age_from, age_to, IDM_data[,4:48], stringsAsFactors = F)
colnames(IDM_data_reformat)[1] = "country"

#Japan census data:
Jap_pop = IDM_data_reformat[IDM_data_reformat$country == "JAP", ]

#get age group data from incidence data:
data_inci_Jap = read.xlsx(paste(data_inci_path, "Japan_JE.xlsx", sep = ""), 1)
Jap_age_group = as.character(data_inci_Jap$Age_group)

#create the Pop_all_age column:
data_inci_Jap$Pop_all_age = pop_age_gen(Jap_age_group, Jap_pop$X2000)
write.xlsx(data_inci_Jap, paste(data_inci_path, "Japan_JE.xlsx", sep = ""), row.names = F)

##################################
#Cambodia:
#Cambodia census data:
Cam_pop = IDM_data_reformat[IDM_data_reformat$country == "KHM", ]
Cam_pop_gen = Cam_pop[, paste("X",year_sel, sep = "")]

#get age group data from incidence data:
data_inci_Cam = read.xlsx(paste(data_inci_path, "Cambodia_JE.xlsx", sep = ""), 1)
#first study:
Cam_subnation_index_1 = c(data_inci_Cam$subnation == unique(data_inci_Cam$subnation)[1])
Cam_age_group_1 = as.character(data_inci_Cam$Age_group[Cam_subnation_index_1])
Cam_pop_all_age_1 = pop_age_gen(Cam_age_group_1, Cam_pop$X2007)
#second study:
#population of siemreap:
Siemreap_pop_2017 = 139458 #data from: http://worldpopulationreview.com/countries/cambodia-population/cities/
Cam_Siemreap_pop_gen = sub_pop_gen(Siemreap_pop_2017, "X2017", Cam_pop_gen)
Cam_subnation_index_2 = c(data_inci_Cam$subnation == unique(data_inci_Cam$subnation)[2])
Cam_age_group_2 = as.character(data_inci_Cam$Age_group[Cam_subnation_index_2])
Cam_pop_all_age_2 = pop_age_gen(Cam_age_group_2, Cam_Siemreap_pop_gen$X2010)
#third study:
Cam_subnation_index_3 = c(data_inci_Cam$subnation == unique(data_inci_Cam$subnation)[3])
Cam_age_group_3 = as.character(data_inci_Cam$Age_group[Cam_subnation_index_3])
Cam_pop_all_age_3 = pop_age_gen(Cam_age_group_3, Cam_pop$X2012)
#fourth study: grouping ages from the first study
Cam_subnation_index_4 = c(data_inci_Cam$subnation == unique(data_inci_Cam$subnation)[4])
Cam_age_group_4 = as.character(data_inci_Cam$Age_group[Cam_subnation_index_4])
Cam_pop_all_age_4 = pop_age_gen(Cam_age_group_4, Cam_pop$X2007)

#create the Pop_all_age column:
data_inci_Cam$Pop_all_age = c(Cam_pop_all_age_1, Cam_pop_all_age_2, Cam_pop_all_age_3)
write.xlsx(data_inci_Cam, paste(data_inci_path, "Cambodia_JE.xlsx", sep = ""), row.names = F)

##################################
#Indonesia:
#Indonesia census data:
Indo_pop = IDM_data_reformat[IDM_data_reformat$country == "IDN", ]
Indo_pop_gen = Indo_pop[age_group_sel + 1, paste("X",year_sel, sep = "")]
#6 provinces of the dataset: East Java, East Nusa Tenggara (West Timor)(2010 only), Papua, West Kalimantan, West Nusa Tenggara (Lombok Island), West Sumatera
Indo_6_prov_pop_2014 = 38529481 + 1662056 + 3486432 + 4073304 + 3311044 + 	5131882
Indo_6_prov_pop_gen = sub_pop_gen(Indo_6_prov_pop_2014, "X2014", Indo_pop_gen)
#get age group data from incidence data:
data_inci_Indo = read.xlsx(paste(data_inci_path, "Indonesia_JE.xlsx", sep = ""), 1)
Indo_subnation_index = c(data_inci_Indo$subnation == unique(data_inci_Indo$subnation)[2])
Indo_age_group = as.character(data_inci_Indo$Age_group[Indo_subnation_index])

#create the Pop_all_age column:
data_inci_Indo$Pop_all_age = c(data_inci_Indo$Pop_all_age[!Indo_subnation_index],pop_age_gen(Indo_age_group, Indo_6_prov_pop_gen$X2005))
write.xlsx(data_inci_Indo, paste(data_inci_path, "Indonesia_JE.xlsx", sep = ""), row.names = F)

###################################
#Laos:
#Laos census data:
Laos_pop = IDM_data_reformat[IDM_data_reformat$country == "LAO", ]
Laos_pop_gen = Laos_pop[, paste("X",year_sel, sep = "")]
#population of Vientian:
Vientian_pop_2017 = 196731 #data from: http://worldpopulationreview.com/countries/laos-population/cities/
Laos_Vientian_pop_gen = sub_pop_gen(Vientian_pop_2017, "X2017", Laos_pop_gen)
#get age group data from incidence data:
data_inci_Laos = read.xlsx(paste(data_inci_path, "Laos_JE.xlsx", sep = ""), 1)
Laos_subnation_index = c(data_inci_Laos$subnation == unique(data_inci_Laos$subnation)[1])
Laos_age_group = as.character(data_inci_Laos$Age_group[Laos_subnation_index])

#create the Pop_all_age column:
data_inci_Laos$Pop_all_age = c(data_inci_Laos$Pop_all_age[!Laos_subnation_index],pop_age_gen(Laos_age_group, Laos_Vientian_pop_gen$X2005))
write.xlsx(data_inci_Laos, paste(data_inci_path, "Laos_JE.xlsx", sep = ""), row.names = F)

##################################
#South Korea:
#South Korea census data:
Kor_pop = IDM_data_reformat[IDM_data_reformat$country == "KOR", ]
Kor_pop_gen = Kor_pop[, paste("X",year_sel, sep = "")]

#get age group data from incidence data:
data_inci_Kor = read.xlsx(paste(data_inci_path, "South_Korea_JE.xlsx", sep = ""), 1)

#cohort study:
Kor_subnation_index_1 = c(data_inci_Kor$subnation == unique(data_inci_Kor$subnation)[1])
Kor_age_group_1 = as.character(data_inci_Kor$Age_group[Kor_subnation_index_1])

Kor_pop_all_age_1_0 = pop_age_gen(Kor_age_group_1, Kor_pop_gen$X2001)
Kor_pop_all_age_1_1 = pop_age_gen(Kor_age_group_1, Kor_pop_gen$X2002)
Kor_pop_all_age_1_2 = pop_age_gen(Kor_age_group_1, Kor_pop_gen$X2003)
Kor_pop_all_age_1_3 = pop_age_gen(Kor_age_group_1, Kor_pop_gen$X2004)
Kor_pop_all_age_1_4 = pop_age_gen(Kor_age_group_1, Kor_pop_gen$X2005)
Kor_pop_all_age_1_5 = pop_age_gen(Kor_age_group_1, Kor_pop_gen$X2006)
Kor_pop_all_age_1_6 = pop_age_gen(Kor_age_group_1, Kor_pop_gen$X2007)
Kor_pop_all_age_1_7 = pop_age_gen(Kor_age_group_1, Kor_pop_gen$X2008)
Kor_pop_all_age_1_8 = pop_age_gen(Kor_age_group_1, Kor_pop_gen$X2009)
Kor_pop_all_age_1_9 = pop_age_gen(Kor_age_group_1, Kor_pop_gen$X2010)
Kor_pop_all_age_1_10 = pop_age_gen(Kor_age_group_1, Kor_pop_gen$X2011)
Kor_pop_all_age_1_11 = pop_age_gen(Kor_age_group_1, Kor_pop_gen$X2012)
Kor_pop_all_age_1_12 = pop_age_gen(Kor_age_group_1, Kor_pop_gen$X2013)
Kor_pop_all_age_1_13 = pop_age_gen(Kor_age_group_1, Kor_pop_gen$X2014)

#create the Pop_all_age column:
data_inci_Kor$Pop_all_age = c(Kor_pop_all_age_1_0, Kor_pop_all_age_1_1, Kor_pop_all_age_1_2, Kor_pop_all_age_1_3, Kor_pop_all_age_1_4, Kor_pop_all_age_1_5, Kor_pop_all_age_1_6, Kor_pop_all_age_1_7, Kor_pop_all_age_1_8, Kor_pop_all_age_1_9, Kor_pop_all_age_1_10, Kor_pop_all_age_1_11, Kor_pop_all_age_1_12, Kor_pop_all_age_1_13)
write.xlsx(data_inci_Kor, paste(data_inci_path, "South_Korea_JE.xlsx", sep = ""), row.names = F)

##################################
#Nepal:
#Nepal census data:
Nepal_pop = naive_pop[naive_pop$country == "NPL",]
Nepal_pop_gen = Nepal_pop[, paste("X",year_sel, sep = "")]

#get age group data from incidence data:
data_inci_Nepal = read.xlsx(paste(data_inci_path, "Nepal_JE.xlsx", sep = ""), 1)

#first study:
Nepal_subnation_index_1 = c(data_inci_Nepal$subnation == unique(data_inci_Nepal$subnation)[1])
Nepal_age_group_1 = as.character(data_inci_Nepal$Age_group[Nepal_subnation_index_1])
#Kath population:
Kathmandu_valley_pop_2005 = 793000 #data from: http://worldpopulationreview.com/world-cities/kathmandu-population/
Nepal_Kath_pop_gen = sub_pop_gen(Kathmandu_valley_pop_2005, "X2005", Nepal_pop_gen)
#pop age gen
Nepal_pop_all_age_1 = pop_age_gen(Nepal_age_group_1, Nepal_Kath_pop_gen$X2007)

#Second study:incomplete
Nepal_subnation_index_2 = c(data_inci_Nepal$subnation == unique(data_inci_Nepal$subnation)[2])
Nepal_age_group_2 = as.character(data_inci_Nepal$Age_group[Nepal_subnation_index_2])
#non_Kath population:
non_Kathmandu_valley_pop_2005 = NA #Hard to get all data:
Nepal_non_Kath_pop_gen = sub_pop_gen(non_Kathmandu_valley_pop_2005, "X2005", Nepal_pop_gen)
#pop age gen
Nepal_pop_all_age_2 = pop_age_gen(Nepal_age_group_2, Nepal_non_Kath_pop_gen$X2007)

#Third and forth studies: pop data included in the paper

#Fifth Study:
Nepal_subnation_index_5 = c(data_inci_Nepal$subnation == unique(data_inci_Nepal$subnation)[5])
Nepal_age_group_5 = as.character(data_inci_Nepal$Age_group[Nepal_subnation_index_5])
Nepal_pop_all_age_5 = pop_age_gen(Nepal_age_group_5, Nepal_pop_gen$X2015)

#Sixth Study:
Nepal_subnation_index_6 = c(data_inci_Nepal$subnation == unique(data_inci_Nepal$subnation)[6])
Nepal_age_group_6 = as.character(data_inci_Nepal$Age_group[Nepal_subnation_index_6])
#pop from 6 districts: Chitwan, Makwanpur, Bara, Parsa, Nawalparasi, Rupandehi
Nepal_6district_2011 = 579984 + 420477 + 687708 + 601017 + 643508 + 880196
Nepal_6district_pop_gen = sub_pop_gen(Nepal_6district_2011, "X2011", Nepal_pop_gen)
#data from: https://data.humdata.org/dataset/nepal-census-2011-district-profiles-demography
Nepal_pop_all_age_6 = pop_age_gen(Nepal_age_group_6, Nepal_6district_pop_gen$X2011)

#Seventh Study:
Nepal_subnation_index_7 = c(data_inci_Nepal$subnation == unique(data_inci_Nepal$subnation)[7])
Nepal_age_group_7 = as.character(data_inci_Nepal$Age_group[Nepal_subnation_index_7])
Nepal_pop_all_age_7 = pop_age_gen(Nepal_age_group_7, Nepal_pop_gen$X2011)

#Eighth Study: (combine of 1st and 2nd study)
Nepal_subnation_index_8 = c(data_inci_Nepal$subnation == unique(data_inci_Nepal$subnation)[8])
Nepal_age_group_8 = as.character(data_inci_Nepal$Age_group[Nepal_subnation_index_8])
#pop from mountain and hill district 
Nepal_mountain_hill_2001 = 10700000 #data from Cambell 2011
Nepal_mountain_hill_pop_gen = sub_pop_gen(Nepal_mountain_hill_2001, "X2001", Nepal_pop_gen)
Nepal_pop_all_age_8 = pop_age_gen(Nepal_age_group_8, Nepal_mountain_hill_pop_gen$X2006)

#create the Pop_all_age column:
#data_inci_Nepal$Pop_all_age = c(Nepal_pop_all_age_1)
#write.xlsx(data_inci_Nepal, paste(data_inci_path, "Nepal_JE.xlsx", sep = ""), row.names = F)

##################################
#Taiwan:
#Taiwan census data:
Taiwan_pop = IDM_data_reformat[IDM_data_reformat$country == "TWN", ]
Taiwan_pop_gen = Taiwan_pop[, paste("X",year_sel, sep = "")]

#get age group data from incidence data:
data_inci_Taiwan = read.xlsx(paste(data_inci_path, "Taiwan_JE.xlsx", sep = ""), 1)
#cohort study:
Taiwan_subnation_index_1 = c(data_inci_Taiwan$subnation == unique(data_inci_Taiwan$subnation)[1])
Taiwan_age_group_1 = as.character(data_inci_Taiwan$Age_group[Taiwan_subnation_index_1])

Taiwan_pop_all_age_1_1 = pop_age_gen(Taiwan_age_group_1, Taiwan_pop_gen$X2002)
Taiwan_pop_all_age_1_2 = pop_age_gen(Taiwan_age_group_1, Taiwan_pop_gen$X2003)
Taiwan_pop_all_age_1_3 = pop_age_gen(Taiwan_age_group_1, Taiwan_pop_gen$X2004)
Taiwan_pop_all_age_1_4 = pop_age_gen(Taiwan_age_group_1, Taiwan_pop_gen$X2005)
Taiwan_pop_all_age_1_5 = pop_age_gen(Taiwan_age_group_1, Taiwan_pop_gen$X2006)
Taiwan_pop_all_age_1_6 = pop_age_gen(Taiwan_age_group_1, Taiwan_pop_gen$X2007)
Taiwan_pop_all_age_1_7 = pop_age_gen(Taiwan_age_group_1, Taiwan_pop_gen$X2008)
Taiwan_pop_all_age_1_8 = pop_age_gen(Taiwan_age_group_1, Taiwan_pop_gen$X2009)
Taiwan_pop_all_age_1_9 = pop_age_gen(Taiwan_age_group_1, Taiwan_pop_gen$X2010)
Taiwan_pop_all_age_1_10 = pop_age_gen(Taiwan_age_group_1, Taiwan_pop_gen$X2011)
Taiwan_pop_all_age_1_11 = pop_age_gen(Taiwan_age_group_1, Taiwan_pop_gen$X2012)

#create the Pop_all_age column:
data_inci_Taiwan$Pop_all_age = c(Taiwan_pop_all_age_1_1, Taiwan_pop_all_age_1_2, Taiwan_pop_all_age_1_3, Taiwan_pop_all_age_1_4, Taiwan_pop_all_age_1_5, Taiwan_pop_all_age_1_6, Taiwan_pop_all_age_1_7, Taiwan_pop_all_age_1_8, Taiwan_pop_all_age_1_9, Taiwan_pop_all_age_1_10, Taiwan_pop_all_age_1_11)
write.xlsx(data_inci_Taiwan, paste(data_inci_path, "Taiwan_JE.xlsx", sep = ""), row.names = F)

###################################
#Vietnam:
#Vietnam census data:
Vietnam_pop = IDM_data_reformat[IDM_data_reformat$country == "VNM", ]
Vietnam_pop_gen = Vietnam_pop[, paste("X",year_sel, sep = "")]

#get age group data from incidence data:
data_inci_Vietnam = read.xlsx(paste(data_inci_path, "Vietnam_JE.xlsx", sep = ""), 1)

#First study:
#population of 12 districts: Hue, KH, DL, B Phuoc, D Thap, A Giang, K Giang, Can Tho, Soc Trang, Bac Lieu, Ca Mau, Tra Vinh
Vietnam_12dist_pop_2011 = (1103.1 + 1174.1 + 1771.8 + 905.3 + 1673.2 + 2151 + 1200.3 + 1303.7 + 873.3 +  1214.9)*1000 #data from: http://www.gso.gov.vn/default_en.aspx?tabid=467&ItemID=12941
Vietnam_12dist_pop_gen = sub_pop_gen(Vietnam_12dist_pop_2011, "X2011", Vietnam_pop_gen)
#get age group
Vietnam_subnation_index_1 = c(data_inci_Vietnam$subnation == unique(data_inci_Vietnam$subnation)[1])
Vietnam_age_group_1 = as.character(data_inci_Vietnam$Age_group[Vietnam_subnation_index_1])
#get pop:
Vietnam_pop_all_gen_1 = pop_age_gen(Vietnam_age_group_1, Vietnam_12dist_pop_gen$X2008)

#Second study:
#get population from 5 provinces:Bac giang, Hai Duong, Hai Phong, Thai Binh, Thanh Hoa.
Vietnam_5provinces_pop_2011 = (1574.3 + 1718.9 + 1878.5 + 1786 + 3412.6)*1000
Vietnam_5provinces_pop_gen = sub_pop_gen(Vietnam_5provinces_pop_2011, "X2011", Vietnam_pop_gen)
#get age group
Vietnam_subnation_index_2 = c(data_inci_Vietnam$subnation == unique(data_inci_Vietnam$subnation)[2])
Vietnam_age_group_2 = as.character(data_inci_Vietnam$Age_group[Vietnam_subnation_index_2])
#gen pop:
Vietnam_pop_all_gen_2 = pop_age_gen(Vietnam_age_group_2, Vietnam_5provinces_pop_gen$X2005)

#create the Pop_all_age column:
data_inci_Vietnam$Pop_all_age = c(Vietnam_pop_all_gen_1, Vietnam_pop_all_gen_2)
write.xlsx(data_inci_Vietnam, paste(data_inci_path, "Vietnam_JE.xlsx", sep = ""), row.names = F)

###################################
#Thailand:
#Thailand census data:
Thailand_pop = IDM_data_reformat[IDM_data_reformat$country == "THA", ]
Thailand_pop_gen = Thailand_pop[, paste("X",year_sel, sep = "")]
#population of 2 districts: Bangkok, Haiyan
Thailand_2dist_pop_2011 = 8280925 + 158218 #data from: https://en.wikipedia.org/wiki/Provinces_of_Thailand
Thailand_2dist_pop_gen = sub_pop_gen(Thailand_2dist_pop_2011, "X2011", Thailand_pop_gen)
#get age group data from incidence data:
data_inci_Thailand = read.xlsx(paste(data_inci_path, "Thailand_JE.xlsx", sep = ""), 1)

Thailand_subnation_index_1 = c(data_inci_Thailand$subnation == unique(data_inci_Thailand$subnation)[1])
Thailand_age_group_1 = as.character(data_inci_Thailand$Age_group[Thailand_subnation_index_1])

#create the Pop_all_age column:
data_inci_Thailand$Pop_all_age = pop_age_gen(Thailand_age_group_1, Thailand_2dist_pop_gen$X2004)
write.xlsx(data_inci_Thailand, paste(data_inci_path, "Thailand_JE.xlsx", sep = ""), row.names = F)

###################################
#Bangladesh:
#Bangladesh census data:
Bangladesh_pop = IDM_data_reformat[IDM_data_reformat$country == "BGD", ]
Bangladesh_pop_gen = Bangladesh_pop[, paste("X",year_sel, sep = "")]

#get age group data from incidence data:
data_inci_Bangladesh = read.xlsx(paste(data_inci_path, "Bangladesh_JE.xlsx", sep = ""), 1)

#first study:
#population of 4 divisions: Dhaka, Mymensingh, Rajshahi, Sylhet
Bangladesh_4div_pop_2011 = 46729000 + 11370000 + 18484858 + 9910219 #data from: https://en.wikipedia.org/wiki/Divisions_of_Bangladesh
Bangladesh_4div_pop_gen = sub_pop_gen(Bangladesh_4div_pop_2011, "X2011", Bangladesh_pop_gen)
#pop all age gen:
Bangladesh_subnation_index_1 = c(data_inci_Bangladesh$subnation == unique(data_inci_Bangladesh$subnation)[1])
Bangladesh_age_group_1 = as.character(data_inci_Bangladesh$Age_group[Bangladesh_subnation_index_1])
Bangladesh_pop_all_age_1 = pop_age_gen(Bangladesh_age_group_1, Bangladesh_4div_pop_gen$X2004)

#second study:
#catchment areas of Rajshahi Medical College Hospital: Rajshahi, Naogaon, Chapai-Nawabgonj districts
Bangladesh_Raj_hos_pop_2011 = 46729000 + 11370000 + 18484858 + 9910219 #data from: https://en.wikipedia.org/wiki/Divisions_of_Bangladesh
Bangladesh_4div_pop_gen = sub_pop_gen(Bangladesh_4div_pop_2011, "X2011", Bangladesh_pop_gen)
#pop all age gen:
Bangladesh_subnation_index_1 = c(data_inci_Bangladesh$subnation == unique(data_inci_Bangladesh$subnation)[1])
Bangladesh_age_group_1 = as.character(data_inci_Bangladesh$Age_group[Bangladesh_subnation_index_1])
Bangladesh_pop_all_age_1 = pop_age_gen(Bangladesh_age_group_1, Bangladesh_4div_pop_gen$X2004)

#create the Pop_all_age column:
data_inci_Bangladesh$Pop_all_age
write.xlsx(data_inci_Bangladesh, paste(data_inci_path, "Bangladesh_JE.xlsx", sep = ""), row.names = F)