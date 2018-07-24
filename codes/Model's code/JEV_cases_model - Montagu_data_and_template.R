rm(list = ls())
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
#Pop data load:#Montagu template load:
output_data_path = 'D:/OUCRU/Hannah/JEV/JEV_model/data_JE_clean/Montagu_data/test3/'
Det_No_vac_template = read.csv(paste(output_data_path,' Det_No_vac_template.csv', sep = ""))
Det_Cam_vac_template = read.csv(paste(output_data_path,' Det_Cam_vac_template.csv', sep = ""))
Det_Rou_vac_template = read.csv(paste(output_data_path,' Det_Rou_vac_template.csv', sep = ""))
Sto_No_vac_template = read.csv(paste(output_data_path,' Sto_No_vac_template.csv', sep = ""))
Sto_Cam_vac_template = read.csv(paste(output_data_path,' Sto_Cam_vac_template.csv', sep = ""))
Sto_Rou_vac_template = read.csv(paste(output_data_path,' Sto_Rou_vac_template.csv', sep = ""))
Sto_para_template = read.csv(paste(output_data_path,"stochastic_template_params.csv", sep = ""))
Sto_para_template = data.frame(run_id = Sto_para_template[,1])
#Pop data load:
pop_data_path = 'D:/OUCRU/Hannah/JEV/JEV_model/data_JE_clean/Montagu_data/test1/'
naive_pop_demo = read.xlsx2(paste(pop_data_path,'naive_pop_2000_2044.xlsx', sep = ""), 1, colClasses=NA)
#cases data load:
data_inci_path = "D:/OUCRU/Hannah/JEV/JEV_model/data_JE_raw/cases_sero_data/"
#save country result path:
coun_result_path = "D:/OUCRU/Hannah/JEV/JEV_model/result_model/first_model/"

#List of country need:
ISO_country_list = as.character(unique(Det_No_vac_template$country))

#year select:
year_sel = 2000:2044
#age group select:
age_group_sel = 0:14
#run times of stochastic model:
run_times = length(unique(Sto_No_vac_template$run_id))

ll_case_je="

data {
  
  int N;  //number of age groups
  int l_age_group; //length of selected age group
  int age_seq[l_age_group]; // age group sequence
  vector[N] age_l; //lower bound of the age group
  vector[N] age_u; //upper bound of the age group
  vector[N] case_age; //number of cases in each age group.
  vector[N] pop_age; //population demo of the study year.
  int year; //number of col of pop_gen_case_less_15: year
  vector[year*l_age_group] pop_gen_case_less_15; //total population that < 15 years old with diff scenario
  vector[year*l_age_group] naive_pop_less_15; //total population that < 15 years old
  real t_case; //total cases in every groups
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
  vector[year*l_age_group] c_less_15; //the generated cases that <15 years old from 2000 to 2044
  vector[year*l_age_group] ir_less_15; //the generated incidence rate that <15 years old from 2000 to 2044
  vector[year] c_total_less_15; //the generated total cases that <15 years old from 2000 to 2044
  vector[year] ir_total_less_15; //the generated total incidence rate that <15 years old from 2000 to 2044
  vector[year*l_age_group] mo_less_15; //the mortality that <15 years old from 2000 to 2044
  vector[year*l_age_group] dis_less_15; //the disability proportion that <15 years old from 2000 to 2044
  vector[year*l_age_group] DALY_less_15; //the DALY that <15 years old from 2000 to 2044 = YLL + YLDacute + YLDchronic
  
  c_less_15 = uniform_rng(1.0/500, 1.0/250)*(1 - e()^(-lambda))*pop_gen_case_less_15 .* to_vector(rep_matrix(exp((to_vector(age_seq)' - 1)*(-lambda)), year));
  //c_less_15 = rho*(1 - e()^(-lambda))*pop_gen_case_less_15 .* to_vector(rep_matrix(exp((to_vector(age_seq)' - 1)*(-lambda)), year));
  c_total_less_15 = rep_vector(0.0, year);
  for(year_sum in 1:year){
    for(age_sum in 1:l_age_group){
      c_total_less_15[year_sum] = c_total_less_15[year_sum] + c_less_15[(age_sum - 1)*year + year_sum];
    }
  }

  ir_less_15 = c_less_15 ./ naive_pop_less_15*100000;
  ir_total_less_15 = rep_vector(0.0, year);
  for(year_sum in 1:year){
    for(age_sum in 1:l_age_group){
      ir_total_less_15[year_sum] = ir_total_less_15[year_sum] + ir_less_15[(age_sum - 1)*year + year_sum];
    }
  }

  mo_less_15 = c_less_15*uniform_rng(0.1,0.3);
  dis_less_15 = (c_less_15 - mo_less_15)*uniform_rng(0.3, 0.5);
  DALY_less_15 = mo_less_15*72 + (c_less_15 - mo_less_15)*uniform_rng(0.012,0.024) + dis_less_15*uniform_rng(1.43,44.5);

}
"

stan_sero_model = stan_model(model_code = ll_case_je)

#loading generate data from diff
#Calculate for each country with different scenario:
stan_model_for_each_scenario <- function(data_inci_subnation, pop_gen_case_less_15, naive_pop_less_15, scenario, combine){
  #Calculate the age_l and age_u:
  age_group_split = strsplit(as.character(data_inci_subnation$Age_group),"-")
  age_l = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[1])}))
  age_u = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[2])}))
  
  #Calculate the time taken to conduct the study and when in average the study conducted:
  temp_date_data = strsplit(as.character(data_inci_subnation$Year),"_")
  date_study_taken = lapply(temp_date_data, FUN = f <- function(x){
    as.Date(paste("01/",x[2], sep = ""), "%d/%m/%Y") - 
      as.Date(paste("01/",x[1], sep = ""), "%d/%m/%Y")})
  year_study_taken = unlist(date_study_taken)/365.25
  
  #MCMC:
  data_for_HMC = list(N = nrow(data_inci_subnation), year = length(year_sel),
                      # year_study = year_study_taken, 
                      age_l = age_l, age_u = age_u, l_age_group = length(age_group_sel), age_seq = 1:length(age_group_sel),
                      case_age = data_inci_subnation$Case_Sero/year_study_taken, 
                      pop_age = data_inci_subnation$Pop_all_age, 
                      pop_gen_case_less_15 = pop_gen_case_less_15,
                      naive_pop_less_15 = naive_pop_less_15,
                      t_case = sum(data_inci_subnation$Case_Sero/year_study_taken))
  
  stan_FOI_fit = sampling(object = stan_sero_model, data = data_for_HMC, 
                          chains = 4, iter = 8000)
  #print(stan_FOI_fit)
  sum_data <- data.frame(summary(stan_FOI_fit)$summary)
  
  
  #1st and 2nd column of output data:
  year = rep(year_sel, length(age_group_sel))
  age =  rep(age_group_sel, each = length(year_sel))
  #Cases estimated:
  case_proj_index = grep("c_less_15" ,rownames(sum_data))
  case_proj = sum_data[case_proj_index,]
  case_proj = cbind(year, age, case_proj)
  case_proj$Scenario = scenario
  
  #total Cases estimated:
  total_case_proj_index = grep("c_total_less_15" ,rownames(sum_data))
  total_case_proj = sum_data[total_case_proj_index,]
  total_case_proj = cbind(year_sel, total_case_proj)
  total_case_proj$Scenario = scenario
  
  #incidence rate estimated:
  ir_proj_index = grep("ir_less_15" ,rownames(sum_data))
  ir_proj = sum_data[ir_proj_index,]
  ir_proj = cbind(year, age, ir_proj)
  ir_proj$Scenario = scenario
  
  #total incidence rate estimated:
  total_ir_proj_index = grep("ir_total_less_15" ,rownames(sum_data))
  total_ir_proj = sum_data[total_ir_proj_index,]
  total_ir_proj = cbind(year_sel, total_ir_proj)
  total_ir_proj$Scenario = scenario
  
  #Mortality estimated:
  mo_proj_index = grep("mo_less_15" ,rownames(sum_data))
  mo_proj = sum_data[mo_proj_index,]
  mo_proj = cbind(year, age, mo_proj)
  mo_proj$Scenario = scenario
  
  #Disability number estimated:
  dis_proj_index = grep("dis_less_15" ,rownames(sum_data))
  dis_proj = sum_data[dis_proj_index,]
  dis_proj = cbind(year, age, dis_proj)
  dis_proj$Scenario = scenario
  
  #DALY estimated:
  DALY_proj_index = grep("DALY_less_15" ,rownames(sum_data))
  DALY_proj = sum_data[DALY_proj_index,]
  DALY_proj = cbind(year, age, DALY_proj)
  DALY_proj$Scenario = scenario
  
  #samples lambda:
  lambda_samples = sample(extract(stan_FOI_fit, pars = "lambda")$lambda,run_times)
  
  #Posterior distribution of parameter:
  mcmc_chain = as.array(stan_FOI_fit)
  lambda_plot = mcmc_areas(mcmc_chain, "lambda")
  rho_plot = mcmc_areas(mcmc_chain, "rho")
  i_age_plot = mcmc_areas(mcmc_chain, regex_pars = "i_age") 
  e_age_plot = mcmc_areas(mcmc_chain, regex_pars = "e_age")
  com_post_plot = gridExtra::grid.arrange(lambda_plot, rho_plot, i_age_plot, e_age_plot, nrow = 2, ncol = 2)
  
  #If we only use regions of a country to generate cases, the code below will extract all output from stan object => use to combine different regions into a whole country
  if(combine){
    cases_ext = extract(stan_FOI_fit, pars = "c_less_15")$c_less_15
    death_ext = extract(stan_FOI_fit, pars = "mo_less_15")$mo_less_15
    dalys_ext = extract(stan_FOI_fit, pars = "DALY_less_15")$DALY_less_15
    return(list(cases = case_proj, ir = ir_proj, deaths = mo_proj, disabilities = dis_proj, dalys = DALY_proj, 
                summary = sum_data, samples = list(lambda = lambda_samples), 
                total = list(cases = total_case_proj, ir = total_ir_proj), 
                extract = list(cases = cases_ext, deaths = death_ext, dalys = dalys_ext), 
                posterior_plot = com_post_plot))
  }
  
  return(list(cases = case_proj, ir = ir_proj, deaths = mo_proj, disabilities = dis_proj, dalys = DALY_proj, 
              summary = sum_data, samples = list(lambda = lambda_samples), 
              total = list(cases = total_case_proj, ir = total_ir_proj),
              posterior_plot = com_post_plot))
}

######################
##India:

##Population data of India overtime:
naive_pop_India = naive_pop_demo[naive_pop_demo$country == "IND",]
#Getting province population overtime:
JE_India_inc = read.xlsx("D:/OUCRU/Hannah/JEV/JEV_model/data_JE_raw/cases_sero_data/India_JE_case_timeseries.xlsx",1)
  #data from: http://www.census2011.co.in/states.php
province_pop_data = data.frame(matrix(cbind(c(as.character(JE_India_inc$Affected.States[-nrow(JE_India_inc)]),"Pondicherry"),
                                             c(84580777,1383727,31205576, 104099452, 16787941, 1458545,25351462,
                                               32988134, 61095297, 33406061, 112374333, 2855794, 2966889, 1978502,
                                               41974218, 27743338, 72147030, 35286757, 3673917, 99812341, 10086292, 
                                               91276115, 241773)),
                                       nrow = 23, ncol = 2))
  #function to generate subnation population over time by the national population data and subnation population of a certain time:
sub_pop_gen <-function(sub_pop_year, year, naive_pop_country){
  naive_pop = naive_pop_country[,4:ncol(naive_pop_country)]
  naive_pop_prop_year = colSums(naive_pop)/colSums(naive_pop)[year]
  naive_pop_sub_year = naive_pop_prop_year*sub_pop_year
  naive_pop_prop_age = sweep(naive_pop, 2, colSums(naive_pop), FUN = "/")
  naive_pop_sub_age_year = sweep(naive_pop_prop_age, 2, naive_pop_sub_year, FUN = "*")
  return(naive_pop_sub_age_year[1:15,])
}

##Lowest incidence population:
lowest_prov = c("Haryana","Punjab")
pop_lowest_prov = sum(as.numeric(as.character(province_pop_data$X2[province_pop_data$X1 %in% lowest_prov])))
pop_lowest_gen = sub_pop_gen(pop_lowest_prov, "X2011", naive_pop_India)
pop_lowest_gen_naive = as.vector(t(pop_lowest_gen))
#get the function of generating campaign and routine population:
source("D:/OUCRU/Hannah/JEV/JEV_model/code_JEV_model/JEV_Montagu_data - Scen_pop_generate_function.R")
pop_lowest_gen_cam = as.vector(t(cam_scen_pop_gen("IND",pop_lowest_gen)))
pop_lowest_gen_rou = as.vector(t(rou_scen_pop_gen("IND",pop_lowest_gen)))
#incidence data to generate Lowest FOI:
data_inci_lowest = read.xlsx(paste(data_inci_path, "Japan", "_JE.xlsx", sep = ""), 1)
subnation_lowest = unlist(unique(data_inci_lowest$subnation))
data_inci_subnation_lowest = data_inci_lowest[data_inci_lowest$subnation == subnation_lowest[1],]
#generate lowest inci in India:
naive_pop = pop_lowest_gen_naive
Ind_lowest_naive_scen = stan_model_for_each_scenario(data_inci_subnation_lowest, pop_lowest_gen_naive, naive_pop, "Naive", T)
Ind_lowest_cam_scen = stan_model_for_each_scenario(data_inci_subnation_lowest, pop_lowest_gen_cam, naive_pop, "Campaign", T)
Ind_lowest_rou_scen = stan_model_for_each_scenario(data_inci_subnation_lowest, pop_lowest_gen_rou, naive_pop, "Routine", T)
#save posterior plot:
ggsave("posterior_plot_l.PNG", plot = Ind_lowest_naive_scen$posterior_plot, path = paste(coun_result_path, "India/", sep = ""), 
       width = 1600/100, height = 800/100)

##Medium to High incidence population:
m_high_prov = c("Andhra Pradesh", "Goa", "Kerala", "Karnataka", "Maharashtra", "Pondicherry", "Tamil Nadu")
pop_m_high_prov = sum(as.numeric(as.character(province_pop_data$X2[province_pop_data$X1 %in% m_high_prov])))
pop_m_high_gen = sub_pop_gen(pop_m_high_prov, "X2011", naive_pop_India)
pop_m_high_gen_naive = as.vector(t(pop_m_high_gen))
#get the function of generating campaign and routine population:
source("D:/OUCRU/Hannah/JEV/JEV_model/code_JEV_model/JEV_Montagu_data - Scen_pop_generate_function.R")
pop_m_high_gen_cam = as.vector(t(cam_scen_pop_gen("IND",pop_m_high_gen)))
pop_m_high_gen_rou = as.vector(t(rou_scen_pop_gen("IND",pop_m_high_gen)))
#incidence data to generate m_high FOI:
data_inci_m_high = read.xlsx(paste(data_inci_path, "Nepal", "_JE.xlsx", sep = ""), 1)
subnation_m_high = unlist(unique(data_inci_m_high$subnation))
data_inci_subnation_m_high = data_inci_m_high[data_inci_m_high$subnation == subnation_m_high[4],]
#data_inci_subnation_m_high = data_inci_subnation_m_high[1:4,]
#generate m_high inci in India:
naive_pop = pop_m_high_gen_naive
Ind_m_high_naive_scen = stan_model_for_each_scenario(data_inci_subnation_m_high, pop_m_high_gen_naive, naive_pop, "Naive", T)
Ind_m_high_cam_scen = stan_model_for_each_scenario(data_inci_subnation_m_high, pop_m_high_gen_cam, naive_pop, "Campaign", T)
Ind_m_high_rou_scen = stan_model_for_each_scenario(data_inci_subnation_m_high, pop_m_high_gen_rou, naive_pop, "Routine", T)
ggsave("posterior_plot_m_high.PNG", plot = Ind_m_high_naive_scen$posterior_plot, path = paste(coun_result_path, "India/", sep = ""), 
       width = 1600/100, height = 800/100)

##High incidence population:
high_prov = c("Assam", "West Bengal", "Bihar", "Manipur", "Uttar Pradesh")
pop_high_prov = sum(as.numeric(as.character(province_pop_data$X2[province_pop_data$X1 %in% high_prov])))
pop_high_gen = sub_pop_gen(pop_high_prov, "X2011", naive_pop_India)
pop_high_gen_naive = as.vector(t(pop_high_gen))
#get the function of generating campaign and routine population:
source("D:/OUCRU/Hannah/JEV/JEV_model/code_JEV_model/JEV_Montagu_data - Scen_pop_generate_function.R")
pop_high_gen_cam = as.vector(t(cam_scen_pop_gen("IND",pop_high_gen)))
pop_high_gen_rou = as.vector(t(rou_scen_pop_gen("IND",pop_high_gen)))
#incidence data to generate high FOI:
data_inci_high = read.xlsx(paste(data_inci_path, "Nepal", "_JE.xlsx", sep = ""), 1)
subnation_high = unlist(unique(data_inci_high$subnation))
data_inci_subnation_high = data_inci_high[data_inci_high$subnation == subnation_high[4],]
#data_inci_subnation_high = data_inci_subnation_high[1:3,]
#generate high inci in India:
naive_pop = pop_high_gen_naive
Ind_high_naive_scen = stan_model_for_each_scenario(data_inci_subnation_high, pop_high_gen_naive, naive_pop, "Naive", T)
Ind_high_cam_scen = stan_model_for_each_scenario(data_inci_subnation_high, pop_high_gen_cam, naive_pop, "Campaign", T)
Ind_high_rou_scen = stan_model_for_each_scenario(data_inci_subnation_high, pop_high_gen_rou, naive_pop, "Routine", T)
ggsave("posterior_plot_high.PNG", plot = Ind_high_naive_scen$posterior_plot, path = paste(coun_result_path, "India/", sep = ""), 
       width = 1600/100, height = 800/100)

######################
##Parkistan:

##Population data of Parkistan overtime:
naive_pop_Parkistan = naive_pop_demo[naive_pop_demo$country == "PAK",]
sindh_prov_gen = sub_pop_gen(29991161, "X2000", naive_pop_Parkistan)
sindh_prov_gen_naive = as.vector(t(sindh_prov_gen))
#get the function of generating campaign and routine population:
source("D:/OUCRU/Hannah/JEV/JEV_model/code_JEV_model/JEV_Montagu_data - Scen_pop_generate_function.R")
pop_sindh_gen_cam = as.vector(t(rou_scen_pop_gen("PAK",sindh_prov_gen)))
pop_sindh_gen_rou = as.vector(t(cam_scen_pop_gen("PAK",sindh_prov_gen)))
#incidence data to generate Lowest FOI:
data_inci_lowest = read.xlsx(paste(data_inci_path, "Japan", "_JE.xlsx", sep = ""), 1)
subnation_lowest = unlist(unique(data_inci_lowest$subnation))
data_inci_subnation_lowest = data_inci_lowest[data_inci_lowest$subnation == subnation_lowest[1],]
#generate lowest inci in India:###fix
naive_pop = sindh_prov_gen_naive
Sind_naive_scen = stan_model_for_each_scenario(data_inci_subnation_lowest, sindh_prov_gen_naive, naive_pop, "Naive", F)
Sind_cam_scen = stan_model_for_each_scenario(data_inci_subnation_lowest, pop_sindh_gen_cam, naive_pop, "Campaign", F)
Sind_rou_scen = stan_model_for_each_scenario(data_inci_subnation_lowest, pop_sindh_gen_rou, naive_pop, "Routine", F)
ggsave("posterior_plot.PNG", plot = Sind_naive_scen$posterior_plot, path = paste(coun_result_path, "Parkistan/", sep = ""), 
       width = 1600/100, height = 800/100)
###################################################

##############PLOTING

###################################################
plotting_proj = function(naive_scen, Cam_scen, Rou_scen){
  case_proj_total = rbind(naive_scen$cases, Cam_scen$cases, Rou_scen$cases)
  cases_plot = ggplot(case_proj_total, aes(x = year, y = mean, group = age)) +
    geom_line(aes(color = age))+
    geom_ribbon(aes(ymax = X97.5., ymin = X2.5., fill = age), alpha = 0.15,
                stat = "identity")+
    facet_grid(.~Scenario) +
    ggtitle("Cases")
  
  inci_rate_proj_total = rbind(naive_scen$ir, Cam_scen$ir, Rou_scen$ir)
  inci_rate_plot = ggplot(inci_rate_proj_total, aes(x = year, y = mean, group = age)) +
    geom_line(aes(color = age))+
    geom_ribbon(aes(ymax = X97.5., ymin = X2.5., fill = age), alpha = 0.15,
                stat = "identity")+
    facet_grid(.~Scenario) +
    ggtitle("Incidence rate")
  
  case_proj_total_sum = rbind(naive_scen$total$cases, Cam_scen$total$cases, Rou_scen$total$cases)
  cases_sum_plot = ggplot(case_proj_total_sum, aes(x = year_sel, y = mean)) +
    geom_line()+
    geom_ribbon(aes(ymax = X97.5., ymin = X2.5.), alpha = 0.15,
                stat = "identity")+
    facet_grid(.~Scenario) +
    ggtitle("Total Cases")
  
  inci_rate_proj_total_sum = rbind(naive_scen$total$ir, Cam_scen$total$ir, Rou_scen$total$ir)
  inci_rate_sum_plot = ggplot(inci_rate_proj_total_sum, aes(x = year_sel, y = mean)) +
    geom_line()+
    geom_ribbon(aes(ymax = X97.5., ymin = X2.5.), alpha = 0.15,
                stat = "identity")+
    facet_grid(.~Scenario) +
    ggtitle("Total Incidence rate")
  
  morta_proj_total = rbind(naive_scen$deaths, Cam_scen$deaths, Rou_scen$deaths)
  morta_plot = ggplot(morta_proj_total, aes(x = year, y = mean, group = age)) +
    geom_line(aes(color = age))+
    geom_ribbon(aes(ymax = X97.5., ymin = X2.5., fill = age), alpha = 0.15,
                stat = "identity")+
    facet_grid(.~Scenario)+
    ggtitle("Mortality")
  
  dis_proj_total = rbind(naive_scen$disabilities, Cam_scen$disabilities, Rou_scen$disabilities)
  dis_plot = ggplot(dis_proj_total, aes(x = year, y = mean, group = age)) +
    geom_line(aes(color = age))+
    geom_ribbon(aes(ymax = X97.5., ymin = X2.5., fill = age), alpha = 0.15,
                stat = "identity")+
    facet_grid(.~Scenario)+
    ggtitle("Disability")
  
  DALY_proj_total = rbind(naive_scen$dalys, Cam_scen$dalys, Rou_scen$dalys)
  dalys_plot = ggplot(DALY_proj_total, aes(x = year, y = mean, group = age)) +
    geom_line(aes(color = age))+
    geom_ribbon(aes(ymax = X97.5., ymin = X2.5., fill = age), alpha = 0.15,
                stat = "identity")+
    facet_grid(.~Scenario)+
    ggtitle("DALYS")
  return(list(cases_plot = cases_plot, incidence_rate_plot = inci_rate_plot, 
              deaths_plot = morta_plot, disability_plot = dis_plot, DALYS_plot = dalys_plot, 
              Sum_cases_plot = cases_sum_plot, Sum_incidence_rate_plot = inci_rate_sum_plot))
}

#India lowest plot:
India_lowest_plot = plotting_proj(Ind_lowest_naive_scen, Ind_lowest_cam_scen, Ind_lowest_rou_scen)

#India mid high plot:
India_m_high_plot = plotting_proj(Ind_m_high_naive_scen, Ind_m_high_cam_scen, Ind_m_high_rou_scen)

#India high plot:
India_high_plot = plotting_proj(Ind_high_naive_scen, Ind_high_cam_scen, Ind_high_rou_scen)

#Parkistan plot:
Parkistan_plot = plotting_proj(Sind_naive_scen, Sind_cam_scen, Sind_rou_scen)

##################################################

#################Combine data of regions into a whole countries

##################################################
#function to combine data:
combine_data <- function(list_data, scenario){
  #list_data = list(...)
  #cases:
  list_cases = lapply(list_data, FUN = f <- function(x)x$extract$cases)
  sum_cases = Reduce("+", list_cases)
  #colnames(sum_cases) = paste("c_less_15", "[",1:ncol(sum_cases),"]", sep = "")
  
  #deaths:
  list_deaths = lapply(list_data, FUN = f <- function(x){x$extract$deaths})
  sum_deaths = Reduce("+", list_deaths)
  #colnames(sum_deaths) = paste("mo_less_15", "[",1:ncol(sum_deaths),"]", sep = "")
  
  #dalys:
  list_dalys = lapply(list_data, FUN = f <- function(x){x$extract$dalys})
  sum_dalys = Reduce("+", list_dalys)
  #colnames(sum_dalys) = paste("DALY_less_15", "[",1:ncol(sum_dalys),"]", sep = "")
  
  summary_data_function = function(data) {data.frame(mean = colMeans(data), 
                            X2.5. = apply(data, 2, f <- function(x)quantile(x, 0.025)),
                            X97.5. = apply(data, 2, f <- function(x)quantile(x, 0.975)))}
  
  #1st and 2nd column of output data:
  year = rep(year_sel, length(age_group_sel))
  age =  rep(age_group_sel, each = length(year_sel))
  #Cases estimated:
  case_proj = summary_data_function(sum_cases)
  case_proj = cbind(year, age, case_proj)
  case_proj$Scenario = scenario
  
  #Mortality estimated:
  mo_proj = summary_data_function(sum_deaths)
  mo_proj = cbind(year, age, mo_proj)
  mo_proj$Scenario = scenario
  
  #DALY estimated:
  DALY_proj = summary_data_function(sum_dalys)
  DALY_proj = cbind(year, age, DALY_proj)
  DALY_proj$Scenario = scenario
  
  return(list(cases = case_proj, deaths = mo_proj, dalys = DALY_proj))
}

#India:
com_India_naive_scen = combine_data(list(Ind_lowest_naive_scen, Ind_m_high_naive_scen, Ind_high_naive_scen), "Naive")
com_India_cam_scen = combine_data(list(Ind_lowest_cam_scen, Ind_m_high_cam_scen, Ind_high_cam_scen), "Campaign")
com_India_rou_scen = combine_data(list(Ind_lowest_rou_scen, Ind_m_high_rou_scen, Ind_high_rou_scen), "Routine")

com_India_plot = plotting_proj(com_India_naive_scen, com_India_cam_scen, com_India_rou_scen)
com_India_plot$cases_plot
com_India_plot$deaths_plot
com_India_plot$DALYS_plot

##################################################

#################Generating stochastic model run_id

##################################################
#set random variable equal for every country, age group, 
rho_samples = runif(run_times, 1/500, 1/250)
death_prop_samples = runif(run_times, 0.1, 0.3)
disability_prop_samples = runif(run_times, 0.3, 0.5)
dalys1_prop_samples = runif(run_times, 0.012, 0.024)
dalys2_prop_samples = runif(run_times, 1.43, 44.5)

sto_model_gen <-function(model_scen, pop_gen_case_less_15){
  lambda_samples = model_scen$samples$lambda
  l_data_output = length(pop_gen_case_less_15)
  cases_gen = rep(rho_samples*(1 - exp(-lambda_samples)), l_data_output)*rep(pop_gen_case_less_15, each = run_times)*
    exp(rep(age_group_sel, each = length(year_sel)*run_times)*(-rep(lambda_samples, l_data_output)))
  mort_gen = cases_gen*rep(death_prop_samples, l_data_output)
  disa_gen = (cases_gen - mort_gen)*rep(disability_prop_samples, l_data_output)
  DALYS_gen = mort_gen*72 + (cases_gen - mort_gen)*rep(dalys1_prop_samples, l_data_output) + 
    disa_gen*rep(dalys2_prop_samples, l_data_output)
  return(list(lambda = lambda_samples, rho = rho_samples, death_prop = death_prop_samples, 
              disability_prop = disability_prop_samples, dalys1_prop = dalys1_prop_samples, dalys2_prop = dalys2_prop_samples,
              cases = cases_gen, deaths = mort_gen, dalys = DALYS_gen))
}

com_sto_model_gen <- function(list_model, list_pop_gen_case_less_15){
  com_list = list()
  for(list in 1:length(list_model)){
    com_list = c(com_list, list(sto_model_gen(list_model[[list]], list_pop_gen_case_less_15[[list]])))
  }
  cases_sum = Reduce("+",lapply(com_list, FUN = f<-function(x)x$cases))
  deaths_sum = Reduce("+",lapply(com_list, FUN = f<-function(x)x$deaths))
  dalys_sum = Reduce("+",lapply(com_list, FUN = f<-function(x)x$dalys))
  lambda_com = Reduce(data.frame,lapply(com_list, FUN = f<-function(x)x$lambda))
  return(list(cases = cases_sum, deaths = deaths_sum, dalys = dalys_sum, 
              lambda = lambda_com))
}
#India:
Ind_naive_sto_gen = com_sto_model_gen(list(Ind_lowest_naive_scen, Ind_m_high_naive_scen, Ind_high_naive_scen), 
                                      list(pop_lowest_gen_naive, pop_m_high_gen_naive, pop_high_gen_naive))
Ind_cam_sto_gen = com_sto_model_gen(list(Ind_lowest_naive_scen, Ind_m_high_naive_scen, Ind_high_naive_scen), 
                                      list(pop_lowest_gen_cam, pop_m_high_gen_cam, pop_high_gen_cam))
Ind_rou_sto_gen = com_sto_model_gen(list(Ind_lowest_naive_scen, Ind_m_high_naive_scen, Ind_high_naive_scen), 
                                      list(pop_lowest_gen_rou, pop_m_high_gen_rou, pop_high_gen_rou))
#Pakistan:
Pak_naive_sto_gen = sto_model_gen(Sind_naive_scen, sindh_prov_gen_naive)
Pak_cam_sto_gen = sto_model_gen(Sind_naive_scen, pop_sindh_gen_cam)
Pak_rou_sto_gen = sto_model_gen(Sind_naive_scen, pop_sindh_gen_rou)

##################################################

################Fill in template:

##################################################
fill_in_Det_temp <- function(Det_template, list_Det_scen_data){
  deaths_vector = c()
  dalys_vector = c()
  cases_vector = c()
  for(i in 1:length(list_Det_scen_data)){
    deaths_vector = c(deaths_vector, list_Det_scen_data[[i]]$deaths$mean)
    dalys_vector = c(dalys_vector, list_Det_scen_data[[i]]$dalys$mean)
    cases_vector = c(cases_vector, list_Det_scen_data[[i]]$cases$mean)
  }
  Det_template$deaths = deaths_vector
  Det_template$dalys = dalys_vector
  Det_template$cases = cases_vector
  return(Det_template)
}
fill_in_Sto_temp <- function(Sto_template, list_Sto_scen_data){
  deaths_vector = c()
  dalys_vector = c()
  cases_vector = c()
  for(ii in 1:length(list_Sto_scen_data)){
    deaths_vector = c(deaths_vector, list_Sto_scen_data[[ii]]$deaths)
    dalys_vector = c(dalys_vector, list_Sto_scen_data[[ii]]$dalys)
    cases_vector = c(cases_vector, list_Sto_scen_data[[ii]]$cases)
  }
  Sto_template$deaths = deaths_vector
  Sto_template$dalys = dalys_vector
  Sto_template$cases = cases_vector
  return(Sto_template)
}
fill_in_Sto_par <- function(Sto_para_template, list_Sto_para, v_model_name){
  data_Sto_lambda = Reduce(data.frame, lapply(list_Sto_para, FUN = f <-function(x)x$lambda))
  filled_Sto_para_template = data.frame(Sto_para_template, data_Sto_lambda, rho_samples, death_prop_samples, 
                                        disability_prop_samples, dalys1_prop_samples, dalys2_prop_samples)
  colnames(filled_Sto_para_template) = c("run_id", paste("lambda", v_model_name, sep = "_"), "rho", "death_prop",
                                         "disability_prop", "dalys1_prop", "dalys2_prop")                                     
  return(filled_Sto_para_template)
}

countries_order = as.character(unique(Det_No_vac_template$country))

Det_No_vac = fill_in_Det_temp(Det_No_vac_template, list(Sind_naive_scen, com_India_naive_scen))
Det_Cam_vac = fill_in_Det_temp(Det_Cam_vac_template, list(Sind_cam_scen, com_India_cam_scen))
Det_Rou_vac = fill_in_Det_temp(Det_Rou_vac_template, list(Sind_rou_scen, com_India_rou_scen))

Sto_No_vac = fill_in_Sto_temp(Sto_No_vac_template, list(Pak_naive_sto_gen, Ind_naive_sto_gen))
Sto_Cam_vac = fill_in_Sto_temp(Sto_Cam_vac_template, list(Pak_cam_sto_gen, Ind_cam_sto_gen))
Sto_Rou_vac = fill_in_Sto_temp(Sto_Rou_vac_template, list(Pak_rou_sto_gen, Ind_rou_sto_gen))

list_Sto_para = list(Pak_naive_sto_gen, Ind_naive_sto_gen)
v_model_name = c("PAK", "IND_lowest", "IND_medium_high", "IND_high")
Sto_para = fill_in_Sto_par(Sto_para_template, list_Sto_para, v_model_name)

#save filled in template:

Det_template_save_path = "D:/OUCRU/Hannah/JEV/JEV_model/result_model/first_model/Det_template/Det"
write.csv(Det_No_vac, file = paste(Det_template_save_path, "_No_vac.csv", sep = ""), row.names = F)
write.csv(Det_Cam_vac, file = paste(Det_template_save_path, "_Cam_vac.csv", sep = ""), row.names = F)
write.csv(Det_Rou_vac, file = paste(Det_template_save_path, "_Rou_vac.csv", sep = ""), row.names = F)

Sto_template_save_path = "D:/OUCRU/Hannah/JEV/JEV_model/result_model/first_model/Sto_template/stochastic_burden_template_JE-OUCRU-Clapham"
write.csv(Sto_No_vac, file = paste(Sto_template_save_path, "-No-vaccination.csv", sep = ""), row.names = F)
write.csv(Sto_Cam_vac, file = paste(Sto_template_save_path, "-Campaign.csv", sep = ""), row.names = F)
write.csv(Sto_Rou_vac, file = paste(Sto_template_save_path, "-Routine.csv", sep = ""), row.names = F)
write.csv(Sto_para, file = paste(Sto_template_save_path, "-parameters.csv", sep = ""), row.names = F)

