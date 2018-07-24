#Calculate FOI from incidence data, then generates cases, health burden overtime.
#Input: country_inci (format: colnames: "Year, subnation, ISO, Age_group, Case_Sero"); 
#       pop_demo (format: age group of 1);
#       vaccine_data ;
#library loading:
library(rstan)
library(bayesplot)
library(ggplot2)
library(xlsx)

#data loading:
#Pop data load:
output_data_path = 'D:/OUCRU/Hannah/JEV/JEV_model/data_JE_clean/pop_data/'
naive_pop = read.xlsx2(paste(output_data_path,'IDB_0_100_2000_2045_24.xlsx', sep = ""),1)
Rou_Best_demo_2000_2045 = read.xlsx2(paste(output_data_path, "Rou_Best_demo_2000_2045_101.xlsx", sep = ""),1)
Rou_NG_demo_2000_2045 = read.xlsx2(paste(output_data_path, "Rou_NG_demo_2000_2045_101.xlsx", sep = ""), 1)
Cam_Best_demo_2000_2045 = read.xlsx2(paste(output_data_path, "Cam_Best_demo_2000_2045_101.xlsx", sep = ""), 1)
Cam_NG_demo_2000_2045 = read.xlsx2(paste(output_data_path, "Cam_NG_demo_2000_2045_101.xlsx", sep = ""), 1)
#cases data load:
data_inci_path = "D:/OUCRU/Hannah/JEV/JEV_model/data_JE_raw/cases_sero_data/"
#List of country need:
ISO_country_list = as.character(unique(naive_pop$ISO))

ll_case_je="

data {
  
  int N;  //number of age groups
  vector[N] age_l; //lower bound of the age group
  vector[N] age_u; //upper bound of the age group
  vector[N] case_age; //number of cases in each age group.
  vector[N] pop_age; //population demo of the study year.
  int a; //number of row of pop_gen_case_less_15: age group
  int b; //number of col of pop_gen_case_less_15: year
  matrix[a,b] pop_gen_case_less_15; //population that <15 years old from 2000 to 2045 to generate case data. 
  vector[b] pop_less_15; //total population that < 15 years old
  int t_case; //total cases in every groups
  //real year_study; // the year taken to conduct the study.

}

parameters {

  real<lower = 0> lambda; //constant FOI
  real<lower = 0> rho;

}

transformed parameters {

  vector[N] i_age; //the expected proportion of infections in each age group
  vector[N] e_age; //The expected number of cases in each age group
  
  i_age = exp(-lambda*age_l) - exp(-lambda*(age_u + 1));
  e_age = (pop_age .* i_age)*rho;

}

model {

  real l_MN;
  //prior distribution
  lambda ~ normal(0, 1000);
  rho ~ uniform(0, 1);
  
  //MN likelihood function:
  l_MN = lgamma(t_case + 1) - sum(lgamma(case_age + 1)) + sum(case_age .* log(e_age/sum(e_age))); 
  
  //likelihood function, included poisson for total cases across all age group:
  target += l_MN + t_case*log(sum(e_age)) - sum(e_age) - lgamma(t_case + 1);

}

generated quantities {
  vector[b] c_less_15; //the generated cases that <15 years old from 2000 to 2045
  vector[b] ir_less_15; //the generated incidence rate that <15 years old from 2000 to 2045
  vector[b] mo_less_15; //the mortality that <15 years old from 2000 to 2045
  vector[b] dis_less_15; //the disability proportion that <15 years old from 2000 to 2045
  vector[b] DALY_less_15; //the DALY that <15 years old from 2000 to 2045 = YLL + YLDacute + YLDchronic

  for(i in 1:b){
    c_less_15[i] = rho*sum((1 - e()^(-lambda))*pop_gen_case_less_15[,i]);
  }
  
  ir_less_15 = c_less_15 ./ pop_less_15*100000;
  mo_less_15 = c_less_15*uniform_rng(0.1,0.3);
  dis_less_15 = (c_less_15 - mo_less_15)*uniform_rng(0.3, 0.5);
  DALY_less_15 = mo_less_15*72 + (c_less_15 - mo_less_15)*uniform_rng(0.012,0.024) + dis_less_15*uniform_rng(1.43,44.5);

}
"

stan_sero_model = stan_model(model_code = ll_case_je)

data_inci = read.xlsx(paste(data_inci_path, "India_JE.xlsx", sep = ""), 1)
subnation = unlist(unique(data_inci$subnation))
data_inci_subnation = data_inci[data_inci$subnation == subnation[1],]
#data_inci_subnation = data_inci_subnation[1:(nrow(data_inci_subnation) - 1), ] #change depend on country

#Calculate the time taken to conduct the study and when in average the study conducted:
temp_date_data = strsplit(as.character(data_inci_subnation$Year),"_")
date_study_taken = lapply(temp_date_data, FUN = f <- function(x){
  as.Date(paste("01/",x[2], sep = ""), "%d/%m/%Y") - 
  as.Date(paste("01/",x[1], sep = ""), "%d/%m/%Y")})
year_study_taken = unlist(date_study_taken)/365.25

year_conduct = lapply(temp_date_data, FUN = f <- function(x){
    data.table::year(as.Date(paste("01/",x[1], sep = ""), "%d/%m/%Y")) + 
    floor(0.5*(data.table::year(as.Date(paste("01/",x[2], sep = ""), "%d/%m/%Y")) - 
    data.table::year(as.Date(paste("01/",x[1], sep = ""), "%d/%m/%Y"))))
  })
year_conduct = year_conduct[[1]]

#Calculate the age_l and age_u:
age_group_split = strsplit(as.character(data_inci_subnation$Age_group),"-")
age_l = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[1])}))
age_u = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[2])}))

#Calculate for each country with different scenario:
stan_model_for_each_scenario <- function(source_pop_data, naive_pop_data ,scenario){
  #group age group depend on age group of the incidence data:
  age_range_total = data.frame(character(),character())
  for(age in 1:length(age_l)){
    age_range = age_u[age] - age_l[age] 
    age_info = c(age_group_split[[age]][1], age_range)
    age_range_total = rbind(age_range_total, age_info, stringsAsFactors = F)
  }
  
  choosen_country = as.character(data_inci_subnation$ISO)[1]
  pop_temp = source_pop_data[source_pop_data$ISO == choosen_country,] #change depend on country
  naive_pop_country = naive_pop_data[naive_pop_data$ISO == choosen_country,] #change depend on country 
  age_start = match(age_range_total$X.00., gsub("-.*$", "",x = pop_temp$Age_group))
    
  pop_gen_case = pop_temp[0,4:49]
  sum_pop = naive_pop_country[0,]
  for(age_s in 1:length(age_start)){
    age_stop = as.numeric(age_range_total[age_s, 2]) - 1
    sel_col_index = age_start[age_s]:(age_start[age_s] + age_stop)
    
    pop_gen_sel_col = sapply(pop_temp[sel_col_index,4:49], FUN = f <- function(x) as.numeric(as.character(x)))
    sum_pop_sel_col = sapply(naive_pop_country[sel_col_index,4:49], FUN = f <- function(x) as.numeric(as.character(x)))
    
    if(age_stop == 0){
      pop_gen_case = rbind(pop_gen_case,pop_gen_sel_col)
      sum_pop = rbind(sum_pop,sum_pop_sel_col)
    } else {
      pop_gen_case = rbind(pop_gen_case,colSums(pop_gen_sel_col))
      sum_pop = rbind(sum_pop,colSums(sum_pop_sel_col))
    }
  }
  colnames(pop_gen_case) = colnames(pop_temp[0,4:49])
  less_15_index = as.numeric(row.names(sum_pop))[as.numeric(row.names(sum_pop)) <= 15]
  sum_pop_less_15 = colSums(sum_pop[less_15_index,])
  
  #Get the demo data from the study year:
  pop_data_index = grep(paste("X", year_conduct, sep = ""), colnames(pop_gen_case))
  pop_data = pop_gen_case[,pop_data_index]
  #Population demo data that < 15:
  pop_gen_case_less_15 = pop_gen_case[less_15_index,]
  
  #MCMC:
  data_for_HMC = list(N = nrow(data_inci_subnation), a = nrow(pop_gen_case_less_15), b = ncol(pop_gen_case_less_15),
                      # year_study = year_study_taken, 
                      age_l = age_l, age_u = age_u, 
                      case_age = data_inci_subnation$Case_Sero, 
                      pop_age = pop_data, pop_gen_case_less_15 = pop_gen_case_less_15, pop_less_15 = sum_pop_less_15,
                      t_case = sum(data_inci_subnation$Case_Sero))
  
  stan_FOI_fit = sampling(object = stan_sero_model, data = data_for_HMC, 
                          chains = 4, iter = 4000)
  print(stan_FOI_fit)
  sum_data = data.frame(summary(stan_FOI_fit)$summary)
  
  #Cases estimated:
  case_proj_index = grep("c_less_15" ,rownames(sum_data))
  case_proj = sum_data[case_proj_index,]
  case_proj$year = 2000:2045
  case_proj$Scenario = scenario
  
  #incidence rate estimated:
  ir_proj_index = grep("ir_less_15" ,rownames(sum_data))
  ir_proj = sum_data[ir_proj_index,]
  ir_proj$year = 2000:2045
  ir_proj$Scenario = scenario
  
  #Mortality estimated:
  mo_proj_index = grep("mo_less_15" ,rownames(sum_data))
  mo_proj = sum_data[mo_proj_index,]
  mo_proj$year = 2000:2045
  mo_proj$Scenario = scenario
  
  #Disability number estimated:
  dis_proj_index = grep("dis_less_15" ,rownames(sum_data))
  dis_proj = sum_data[dis_proj_index,]
  dis_proj$year = 2000:2045
  dis_proj$Scenario = scenario
  
  #DALY estimated:
  DALY_proj_index = grep("DALY_less_15" ,rownames(sum_data))
  DALY_proj = sum_data[DALY_proj_index,]
  DALY_proj$year = 2000:2045
  DALY_proj$Scenario = scenario
  
  return(list(case_proj, ir_proj, mo_proj, dis_proj, DALY_proj))
}

total_naive_pop_less_15 = naive_pop

source_pop_data = naive_pop
naive_scen = stan_model_for_each_scenario(source_pop_data, total_naive_pop_less_15, "Naive")

source_pop_data = Rou_Best_demo_2000_2045
Rou_B_scen = stan_model_for_each_scenario(source_pop_data, total_naive_pop_less_15, "Routine Best")

source_pop_data = Rou_NG_demo_2000_2045
Rou_NG_scen = stan_model_for_each_scenario(source_pop_data, total_naive_pop_less_15, "Routine No GAVI")

source_pop_data = Cam_Best_demo_2000_2045
Cam_B_scen = stan_model_for_each_scenario(source_pop_data, total_naive_pop_less_15, "Campaign Best")

source_pop_data = Cam_NG_demo_2000_2045
Cam_NG_scen = stan_model_for_each_scenario(source_pop_data, naive_pop, "Campaign No GAVI")

case_proj_total = rbind(naive_scen[[1]], Rou_B_scen[[1]], Rou_NG_scen[[1]], Cam_B_scen[[1]], Cam_NG_scen[[1]])
ggplot(case_proj_total, aes(x = year, y = mean, group = Scenario)) +
  geom_boxplot(aes(ymax = X97.5., ymin = X2.5., 
                   lower = X25., upper = X75., middle = X50., fill = Scenario),
               stat = "identity")

inci_rate_proj_total = rbind(naive_scen[[2]], Rou_B_scen[[2]], Rou_NG_scen[[2]], Cam_B_scen[[2]], Cam_NG_scen[[2]])
ggplot(inci_rate_proj_total, aes(x = year, y = mean, group = Scenario)) +
  geom_boxplot(aes(ymax = X97.5., ymin = X2.5., 
                   lower = X25., upper = X75., middle = X50., fill = Scenario),
               stat = "identity")

morta_proj_total = rbind(naive_scen[[3]], Rou_B_scen[[3]], Rou_NG_scen[[3]], Cam_B_scen[[3]], Cam_NG_scen[[3]])
ggplot(morta_proj_total, aes(x = year, y = mean, group = Scenario)) +
  geom_boxplot(aes(ymax = X97.5., ymin = X2.5., 
                   lower = X25., upper = X75., middle = X50., fill = Scenario),
               stat = "identity")

dis_proj_total = rbind(naive_scen[[4]], Rou_B_scen[[4]], Rou_NG_scen[[4]], Cam_B_scen[[4]], Cam_NG_scen[[4]])
ggplot(dis_proj_total, aes(x = year, y = mean, group = Scenario)) +
  geom_boxplot(aes(ymax = X97.5., ymin = X2.5., 
                   lower = X25., upper = X75., middle = X50., fill = Scenario),
               stat = "identity")

DALY_proj_total = rbind(naive_scen[[5]], Rou_B_scen[[5]], Rou_NG_scen[[5]], Cam_B_scen[[5]], Cam_NG_scen[[5]])
ggplot(DALY_proj_total, aes(x = year, y = mean, group = Scenario)) +
  geom_boxplot(aes(ymax = X97.5., ymin = X2.5., 
                   lower = X25., upper = X75., middle = X50., fill = Scenario),
               stat = "identity")
