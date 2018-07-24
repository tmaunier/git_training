#Calculate FOI from incidence data, then generates cases, health burden overtime.
#Input: country_inci (format: colnames: "Year, subnation, ISO, Age_group, Case_Sero, Pop_all_age"); 
#       pop_demo (format: age group of 1);
#       vaccine_data ;
#library loading:
library(rstan)
library(bayesplot)
library(ggplot2)
library(xlsx)

#data loading:
#Pop data load:
output_data_path = 'D:/OUCRU/Hannah/JEV/JEV_model/data_JE_clean/Montagu_data/test1/'
naive_pop = read.xlsx2(paste(output_data_path,'naive_pop_2000_2045.xlsx', sep = ""),1)
naive_pop__sum_less_15_2000_2045 = read.xlsx2(paste(output_data_path,'naive_pop__sum_less_15_2000_2045.xlsx', sep = ""),1)
Rou_demo_2000_2045 = read.xlsx2(paste(output_data_path, "rou_pop__sum_less_15_2000_2045.xlsx", sep = ""),1)
Cam_demo_2000_2045 = read.xlsx2(paste(output_data_path, "cam_pop__sum_less_15_2000_2045.xlsx", sep = ""), 1)
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
  int year; //number of col of pop_gen_case_less_15: year
  vector[year] pop_gen_case_less_15; //total population that < 15 years old with diff scenario
  vector[year] naive_pop_less_15; //total population that < 15 years old
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
  vector[year] c_less_15; //the generated cases that <15 years old from 2000 to 2045
  vector[year] ir_less_15; //the generated incidence rate that <15 years old from 2000 to 2045
  vector[year] mo_less_15; //the mortality that <15 years old from 2000 to 2045
  vector[year] dis_less_15; //the disability proportion that <15 years old from 2000 to 2045
  vector[year] DALY_less_15; //the DALY that <15 years old from 2000 to 2045 = YLL + YLDacute + YLDchronic

  c_less_15 = rho*(1 - e()^(-lambda))*pop_gen_case_less_15;
  ir_less_15 = c_less_15 ./ naive_pop_less_15*100000;
  mo_less_15 = c_less_15*uniform_rng(0.1,0.3);
  dis_less_15 = (c_less_15 - mo_less_15)*uniform_rng(0.3, 0.5);
  DALY_less_15 = mo_less_15*72 + (c_less_15 - mo_less_15)*uniform_rng(0.012,0.024) + dis_less_15*uniform_rng(1.43,44.5);

}
"

stan_sero_model = stan_model(model_code = ll_case_je)

#loading country data
data_inci = read.xlsx(paste(data_inci_path, "India", "_JE.xlsx", sep = ""), 1)
subnation = unlist(unique(data_inci$subnation))
data_inci_subnation = data_inci[data_inci$subnation == subnation[2],]

#loading generate data from diff
#Calculate for each country with different scenario:
stan_model_for_each_scenario <- function(source_pop_data, naive_pop_data ,scenario){
  #Calculate the age_l and age_u:
  age_group_split = strsplit(as.character(data_inci_subnation$Age_group),"-")
  age_l = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[1])}))
  age_u = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[2])}))
  
  #Get demo data from country:
  choosen_country = as.character(unique(data_inci$ISO))
  pop_gen_case_less_15 = source_pop_data[source_pop_data$country_ls == choosen_country, 2:ncol(source_pop_data)]
  naive_pop_less_15 = naive_pop__sum_less_15_2000_2045[naive_pop__sum_less_15_2000_2045$country_ls == choosen_country, 2:ncol(naive_pop__sum_less_15_2000_2045)]
  #MCMC:
  data_for_HMC = list(N = nrow(data_inci_subnation), year = length(pop_gen_case_less_15),
                      # year_study = year_study_taken, 
                      age_l = age_l, age_u = age_u, 
                      case_age = data_inci_subnation$Case_Sero, 
                      pop_age = data_inci_subnation$Pop_all_age, 
                      pop_gen_case_less_15 = as.numeric(as.character(unlist(pop_gen_case_less_15))),
                      naive_pop_less_15 = as.numeric(as.character(unlist(naive_pop_less_15))),
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

source_pop_data = naive_pop__sum_less_15_2000_2045
naive_scen = stan_model_for_each_scenario(source_pop_data, naive_pop, "Naive")

source_pop_data = Cam_demo_2000_2045
Cam_scen = stan_model_for_each_scenario(source_pop_data, naive_pop, "Campaign")

source_pop_data = Rou_demo_2000_2045
Rou_scen = stan_model_for_each_scenario(source_pop_data, naive_pop, "Routine")

case_proj_total = rbind(naive_scen[[1]], Cam_scen[[1]], Rou_scen[[1]])
ggplot(case_proj_total, aes(x = year, y = mean, group = Scenario)) +
  geom_boxplot(aes(ymax = X97.5., ymin = X2.5., 
                   lower = X25., upper = X75., middle = X50., fill = Scenario),
               stat = "identity")

inci_rate_proj_total = rbind(naive_scen[[2]], Cam_scen[[2]], Rou_scen[[2]])
ggplot(inci_rate_proj_total, aes(x = year, y = mean, group = Scenario)) +
  geom_boxplot(aes(ymax = X97.5., ymin = X2.5., 
                   lower = X25., upper = X75., middle = X50., fill = Scenario),
               stat = "identity")

morta_proj_total = rbind(naive_scen[[3]], Cam_scen[[3]], Rou_scen[[3]])
ggplot(morta_proj_total, aes(x = year, y = mean, group = Scenario)) +
  geom_boxplot(aes(ymax = X97.5., ymin = X2.5., 
                   lower = X25., upper = X75., middle = X50., fill = Scenario),
               stat = "identity")

dis_proj_total = rbind(naive_scen[[4]], Cam_scen[[4]], Rou_scen[[4]])
ggplot(dis_proj_total, aes(x = year, y = mean, group = Scenario)) +
  geom_boxplot(aes(ymax = X97.5., ymin = X2.5., 
                   lower = X25., upper = X75., middle = X50., fill = Scenario),
               stat = "identity")

DALY_proj_total = rbind(naive_scen[[5]], Cam_scen[[5]], Rou_scen[[5]])
ggplot(DALY_proj_total, aes(x = year, y = mean, group = Scenario)) +
  geom_boxplot(aes(ymax = X97.5., ymin = X2.5., 
                   lower = X25., upper = X75., middle = X50., fill = Scenario),
               stat = "identity")
