rm(list = ls())
#JEV model: Multinomial + Poisson to estimate case reported, based on Hannah's code
#library loading:
library(rstan)
library(bayesplot)
library(xlsx)
library(ggplot2)
library(dplyr)
{
#####################################################################################################
ll_case_je="
data {

  int N;  //number of age groups
  vector[N] age_l; //lower bound of the age group
  vector[N] age_u; //upper bound of the age group
  vector[N] case_age; //number of cases in each age group.
  vector[N] pop_age; //population demo
  int t_case; //total cases in every groups

}

parameters {

  real<lower = 0> lambda; //constant FOI
  real<lower = 0> rho;
}

transformed parameters{
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
"

stan_sero_model_2 = stan_model(model_code = ll_case_je)

model_run <- function(data){
  data_for_HMC = list(N = nrow(data), age_l = data$age_l, age_u = data$age_u, 
                      case_age = data$cases, pop_age = data$pop, 
                      t_case = sum(data$cases))
  
  stan_FOI_fit = sampling(object = stan_sero_model_2, data = data_for_HMC, 
                          chains = 4, iter = 4000)
  return(stan_FOI_fit)
}

data_inci_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_raw/cases_sero_data/"


#Model for "A hospital-based surveillance for Japanese encephalitis in Bali, Indonesia"
Data_je_indo = data.frame(c("00-01","01-04", "05-09", "10-15"),
                          c(0, 1, 5, 10),
                          c(0, 4, 9, 15),
                          c(15, 43, 20, 4), 
                          c(15145699, 44796619, 68044091, 104012947))
colnames(Data_je_indo) = c("age_group", "age_l", "age_u", "cases", "pop")

stan_FOI_fit <- model_run(Data_je_indo)
print(stan_FOI_fit)

mcmc_areas(as.array(stan_FOI_fit),regex_pars = "i_age")
mcmc_areas(as.array(stan_FOI_fit),regex_pars = "e_age")
mcmc_areas(as.array(stan_FOI_fit),regex_pars = "lambda")

mcmc_trace(as.array(stan_FOI_fit))

mcmc_neff(neff_ratio(stan_FOI_fit))

mcmc_acf(as.array(stan_FOI_fit))

mcmc_pairs(as.array(stan_FOI_fit),regex_pars = "i_age")

sum_data = data.frame(summary(stan_FOI_fit)$summary)
est_cases = sum_data[grep("e_age",rownames(sum_data)),]
est_cases$age_group = Data_je_indo$age_group
  
ggplot(data = Data_je_indo, aes(x = age_group))+
  geom_point(aes(y = cases))+
  geom_path(aes(y = est_cases$mean), group = 1)+
  geom_ribbon(data = est_cases, aes(x = 1:3, ymin = X2.5., ymax = X97.5.), alpha = 0.1)



#Touch paper:
Data_je_touch = read.xlsx("D:/OUCRU/Hannah/JEV/JE_shared_documents/Recreating WHO 2011 estimates/Papers included in their analysis/Touch_paper.xlsx", 1)
colnames(Data_je_touch) = c("age", "cases", "pop")
Data_je_touch$age_group = paste(formatC(Data_je_touch$age, digits = 2, width = 2, flag = "0"),"-",
                          formatC(Data_je_touch$age + 1, digits = 2, width = 2, flag = "0"),sep = "")
data = Data_je_touch
data_for_HMC = list(N = nrow(data), age_l = data$age, age_u = data$age + 1, 
                    case_age = data$cases, pop_age = data$pop, 
                    t_case = sum(data$cases))

stan_FOI_fit = sampling(object = stan_sero_model_2, data = data_for_HMC, 
                        chains = 4, iter = 4000)

print(stan_FOI_fit)

sum_data = data.frame(summary(stan_FOI_fit)$summary)
est_cases = sum_data[grep("e_age",rownames(sum_data)),]
est_cases$age_group = Data_je_touch$age_group

ggplot(data = Data_je_touch, aes(x = age_group))+
  geom_point(aes(y = cases))+
  geom_path(aes(y = est_cases$mean), group = 1)+
  geom_ribbon(data = est_cases, aes(x = 1:nrow(data), ymin = X2.5., ymax = X97.5.), alpha = 0.1)


#Wierzba paper:
Data_je_wierzba = read.xlsx("D:/OUCRU/Hannah/JEV/JE_shared_documents/Recreating WHO 2011 estimates/Wierzba_et_al_data.xlsx", 1)
Data_je_wierzba$numbervaccine = c(0,0,331400.53,74666.67,44800.00,0,0,0,151244.80,0,0,0)
Data_je_wierzba$age_l = c(0, 1, 5, 15, 20, 35)
Data_je_wierzba$age_u = c(1, 5, 15, 20, 35, 70)
Data_je_wierzba$age_group = paste(formatC(Data_je_wierzba$age_l, digits = 2, width = 2, flag = "0"),"-",
                                formatC(Data_je_wierzba$age_u, digits = 2, width = 2, flag = "0"),sep = "")
#Western Terai:
W_Terai_data = Data_je_wierzba[1:6,]

data = W_Terai_data
data_for_HMC = list(N = nrow(data), age_l = data$age_l, age_u = data$age_u, 
                    case_age = data$Cases, pop_age = data$population - data$numbervaccine, 
                    t_case = sum(data$Cases))

stan_FOI_fit = sampling(object = stan_sero_model_2, data = data_for_HMC, 
                        chains = 4, iter = 4000)

summary(stan_FOI_fit)$summary

sum_data = data.frame(summary(stan_FOI_fit)$summary)
est_cases = sum_data[grep("e_age",rownames(sum_data)),]
est_cases$age_group = W_Terai_data$age_group

ggplot(data = W_Terai_data, aes(x = age_group))+
  geom_point(aes(y = Cases))+
  geom_path(aes(y = est_cases$mean), group = 1)+
  geom_ribbon(data = est_cases, aes(x = 1:nrow(data), ymin = X2.5., ymax = X97.5.), alpha = 0.1)

#Non-Western Terai:
Non_W_Terai_data = Data_je_wierzba[7:12,]

data = Non_W_Terai_data
data_for_HMC = list(N = nrow(data), age_l = data$age_l, age_u = data$age_u, 
                    case_age = data$Cases, pop_age = data$population - data$numbervaccine, 
                    t_case = sum(data$Cases))

stan_FOI_fit = sampling(object = stan_sero_model_2, data = data_for_HMC, 
                        chains = 4, iter = 4000)

summary(stan_FOI_fit)$summary

sum_data = data.frame(summary(stan_FOI_fit)$summary)
est_cases = sum_data[grep("e_age",rownames(sum_data)),]
est_cases$age_group = Non_W_Terai_data$age_group

mcmc_areas(as.array(stan_FOI_fit),regex_pars = "i_age")
mcmc_areas(as.array(stan_FOI_fit),regex_pars = "e_age")
mcmc_areas(as.array(stan_FOI_fit),regex_pars = "lambda")

mcmc_trace(as.array(stan_FOI_fit))

mcmc_neff(neff_ratio(stan_FOI_fit))

mcmc_acf(as.array(stan_FOI_fit))

mcmc_pairs(as.array(stan_FOI_fit),regex_pars = "i_age")


ggplot(data = Non_W_Terai_data, aes(x = age_group))+
  geom_point(aes(y = Cases))+
  geom_path(aes(y = est_cases$mean), group = 1)+
  geom_ribbon(data = est_cases, aes(x = 1:nrow(data), ymin = X2.5., ymax = X97.5.), alpha = 0.1)

}

#######################################################################################################

ll_case_je_2="

data {

int N;  //number of age groups
vector[N] age_l; //lower bound of the age group
vector[N] age_u; //upper bound of the age group
vector[N] case_age; //number of cases in each age group.
vector[N] pop_age; //population demo of the study year.
real t_case; //total cases in every groups

}

parameters {

real<lower = 0> lambda;
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
"
stan_case_model = stan_model(model_code = ll_case_je_2)

###################################################
#Model for cohort data:
ll_case_cohort_same_rho_je_2="
data {
  int N;//number of age groups
  int N_cohort;// number of studies in cohort
  int start_cohort[N_cohort];//starting index of the cohort
  int stop_cohort[N_cohort];//stoping index of the cohort
  vector[N] age_l; //lower bound of the age group
  vector[N] age_u; //upper bound of the age group
  vector[N] case_age; //number of cases in each age group.
  vector[N] pop_age; //population demo of the study year.
  real t_case; //total cases in every groups
}

parameters {
  real<lower = 0> lambda[N_cohort];
  real<lower = 0> rho;
}

transformed parameters {
  vector[N] i_age; //the expected proportion of infections in each age group
  vector[N] e_age; //The expected number of cases in each age group
  for(year in 1:N_cohort){
    i_age[(start_cohort[year]):(stop_cohort[year])] = 
    exp(-lambda[year]*age_l[(start_cohort[year]):(stop_cohort[year])]) - 
    exp(-lambda[year]*(age_u[(start_cohort[year]):(stop_cohort[year])] + 1));
  }
  e_age = (pop_age .* i_age)*rho;
}

model {
  real l_MN;
  //prior distribution
  lambda[N_cohort] ~ normal(0, 1000);
  rho ~ uniform(0, 1);
  //MN likelihood function:
  l_MN = lgamma(t_case + 1) - sum(lgamma(case_age + 1)) + sum(case_age .* log(e_age/sum(e_age))); 
  //likelihood function, included poisson for total cases across all age group:
  target += l_MN + t_case*log(sum(e_age)) - sum(e_age) - lgamma(t_case + 1);
}
"
stan_case_cohort_same_rho_model = stan_model(model_code = ll_case_cohort_same_rho_je_2)

###################################################
#Model for cohort data:
ll_case_cohort_same_FOI_je_2="
data {
  int N;//number of age groups
  int N_cohort;// number of studies in cohort
  int start_cohort[N_cohort];//starting index of the cohort
  int stop_cohort[N_cohort];//stoping index of the cohort
  vector[N] age_l; //lower bound of the age group
  vector[N] age_u; //upper bound of the age group
  vector[N] case_age; //number of cases in each age group.
  vector[N] pop_age; //population demo of the study year.
  real t_case; //total cases in every groups
}

parameters {
  real<lower = 0> lambda;
  real<lower = 0> rho[N_cohort];
}

transformed parameters {
  vector[N] i_age; //the expected proportion of infections in each age group
  vector[N] e_age; //The expected number of cases in each age group
  i_age = exp(-lambda*age_l) - exp(-lambda*(age_u + 1));
  for(year in 1:N_cohort){
    e_age[(start_cohort[year]):(stop_cohort[year])] = (pop_age[(start_cohort[year]):(stop_cohort[year])] .* i_age[(start_cohort[year]):(stop_cohort[year])])*rho[year];
  }
  
}

model {
  real l_MN;
  //prior distribution
  lambda ~ normal(0, 1000);
  rho[N_cohort] ~ uniform(0, 1);
  //MN likelihood function:
  l_MN = lgamma(t_case + 1) - sum(lgamma(case_age + 1)) + sum(case_age .* log(e_age/sum(e_age))); 
  //likelihood function, included poisson for total cases across all age group:
  target += l_MN + t_case*log(sum(e_age)) - sum(e_age) - lgamma(t_case + 1);
}
"
stan_case_cohort_same_FOI_model = stan_model(model_code = ll_case_cohort_same_FOI_je_2)

model_run_case <- function(data, model){
  #Calculate the age_l and age_u:
  age_group_split = strsplit(as.character(data$Age_group),"-")
  age_l = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[1])}))
  age_u = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[2])}))
  
  #Calculate the time taken to conduct the study and when in average the study conducted:
  temp_date_data = strsplit(as.character(data$Year),"_")
  date_study_taken = lapply(temp_date_data, FUN = f <- function(x){
    as.Date(paste("01/",x[2], sep = ""), "%d/%m/%Y") - 
      as.Date(paste("01/",x[1], sep = ""), "%d/%m/%Y")})
  year_study_taken = unlist(date_study_taken)/365.25
  
  data_for_HMC = list(N = nrow(data), age_l = age_l, age_u = age_u, 
                      case_age = data$Case_Sero/year_study_taken, pop_age = data$Pop_all_age, 
                      t_case = sum(data$Case_Sero/year_study_taken))
  
  stan_FOI_fit = sampling(object = model, data = data_for_HMC, 
                          chains = 4, iter = 8000)
  
  #Posterior distribution of parameter:
  mcmc_chain = as.array(stan_FOI_fit)
  lambda_plot = mcmc_areas(mcmc_chain, "lambda")
  rho_plot = mcmc_areas(mcmc_chain, regex_pars = "rho")
  i_age_plot = mcmc_areas(mcmc_chain, regex_pars = "i_age") 
  e_age_plot = mcmc_areas(mcmc_chain, regex_pars = "e_age")
  com_post_plot = gridExtra::grid.arrange(lambda_plot, rho_plot, i_age_plot, e_age_plot, nrow = 2, ncol = 2)
  
  return(list(summary = summary(stan_FOI_fit)$summary, plot = com_post_plot))
}
model_run_case_year_sum <- function(data, model){
  #Calculate the age_l and age_u:
  age_group_split = strsplit(as.character(data$Age_group),"-")
  age_l = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[1])}))
  age_u = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[2])}))
  
  data_for_HMC = list(N = nrow(data), age_l = age_l, age_u = age_u, 
                      case_age = data$Case_Sero, pop_age = data$Pop_all_age_year_sum, 
                      t_case = sum(data$Case_Sero))
  
  stan_FOI_fit = sampling(object = model, data = data_for_HMC, 
                          chains = 4, iter = 8000)
  
  #Posterior distribution of parameter:
  mcmc_chain = as.array(stan_FOI_fit)
  lambda_plot = mcmc_areas(mcmc_chain, "lambda")
  rho_plot = mcmc_areas(mcmc_chain, regex_pars = "rho")
  i_age_plot = mcmc_areas(mcmc_chain, regex_pars = "i_age") 
  e_age_plot = mcmc_areas(mcmc_chain, regex_pars = "e_age")
  com_post_plot = gridExtra::grid.arrange(lambda_plot, rho_plot, i_age_plot, e_age_plot, nrow = 2, ncol = 2)
  
  #Summary data:
  sum_data = data.frame(summary(stan_FOI_fit)$summary)
  extra_data = data[,c("Index", "Year", "subnation", "ISO")] %>% unique
  
  return(list(summary = sum_data, model_info = extra_data, plot = com_post_plot))
}
model_run_case_year_sum_cohort_same_rho <- function(data, model){
  #Calculate the age_l and age_u:
  age_group_split = strsplit(as.character(data$Age_group),"-")
  age_l = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[1])}))
  age_u = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[2])}))
  
  l_cohort = unname(table(as.vector(data$subnation)))
  start_cohort = cumsum(c(1,head(l_cohort, -1)))
  stop_cohort = cumsum(l_cohort)
  
  data_for_HMC = list(N = nrow(data), age_l = age_l, age_u = age_u, 
                      case_age = data$Case_Sero, pop_age = data$Pop_all_age_year_sum, 
                      start_cohort = start_cohort, stop_cohort = stop_cohort, N_cohort = length(start_cohort),
                      t_case = sum(data$Case_Sero))
  
  stan_FOI_fit = sampling(object = model, data = data_for_HMC, 
                          chains = 4, iter = 12000)
  
  #Posterior distribution of parameter:
  mcmc_chain = as.array(stan_FOI_fit)
  lambda_plot = mcmc_areas(mcmc_chain, regex_pars = "lambda")
  rho_plot = mcmc_areas(mcmc_chain, regex_pars = "rho")
  i_age_plot = mcmc_areas(mcmc_chain, regex_pars = "i_age") 
  e_age_plot = mcmc_areas(mcmc_chain, regex_pars = "e_age")
  com_post_plot = gridExtra::grid.arrange(lambda_plot, rho_plot, i_age_plot, e_age_plot, nrow = 2, ncol = 2)
  
  #Summary data:
  sum_data = summary(stan_FOI_fit)$summary %>% data.frame
  extra_data = data[,c("Index", "Year", "subnation", "ISO")] %>% unique
  
  return(list(summary = sum_data, model_info = extra_data, plot = com_post_plot))
}

data_inci_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_raw/cases_sero_data/"

##############
###Model fitting for available countries:
##############
##Japan:
Japan_inci_data = read.xlsx(paste(data_inci_path, "Japan_JE.xlsx", sep = ""), 1)
JPN_model_01 <- model_run_case_year_sum(Japan_inci_data[1:10,], stan_case_model) #FOI: 0.000014 =>not fit well
JPN_model_02 <- model_run_case_year_sum(Japan_inci_data[11:15,], stan_case_model) #FOI: 0.003 =>not fit well
JPN_model_03 <- model_run_case_year_sum(Japan_inci_data[16:25,], stan_case_model) #FOI: 0.000015 =>not fit well
JPN_model_04 <- model_run_case_year_sum(Japan_inci_data[26:35,], stan_case_model) #FOI: 0.0001 =>not fit well

#Laos:
#A Prospective Assessment of the Accuracy of Commercial IgM ELISAs in Diagnosis of Japanese Encephalitis Virus Infections in Patients with Suspected Central Nervous System Infections in Laos
Laos_inci_data = read.xlsx(paste(data_inci_path, "Laos_JE.xlsx", sep = ""), 1)
LAO_model <- model_run_case_year_sum(Laos_inci_data, stan_case_model) #FOI: 0.11|FOI: 0.074

#Cam:
Cam_inci_data = read.xlsx(paste(data_inci_path, "Cambodia_JE.xlsx", sep = ""), 1)
KHM_model_01 <- model_run_case_year_sum(Cam_inci_data[1:16,], stan_case_model) #FOI: 0.06 => not fit well => groupping?|FOI: 0.063 => better fit|FOI: 0.064
KHM_model_02 <- model_run_case_year_sum(Cam_inci_data[17:19,], stan_case_model) #FOI: 0.72|FOI: 0.69|FOI: 0.69
KHM_model_03 <- model_run_case_year_sum(Cam_inci_data[20:23,], stan_case_model) #FOI: 0.12 => not fit quite well|FOI:0.099 => fit well|0.102
KHM_model_04 <- model_run_case_year_sum(Cam_inci_data[24:28,], stan_case_model) #FOI: 0.084 => not fit quite well|FOI: 0.086 => fit well|0.088

#South Korea:cohort data
Kor_inci_data = read.xlsx(paste(data_inci_path, "South_Korea_JE.xlsx", sep = ""), 1)
KOR_model_01_01 <- model_run_case_year_sum(Kor_inci_data[1:7,], stan_case_model) #FOI: 0.011|0.016 #rho: 0.04
KOR_model_01_02 <- model_run_case_year_sum(Kor_inci_data[8:14,], stan_case_model) #FOI: 0.033|0.041 #rho: 0.005
KOR_model_01_03 <- model_run_case_year_sum(Kor_inci_data[15:21,], stan_case_model) #FOI: 0|0
KOR_model_01_04 <- model_run_case_year_sum(Kor_inci_data[22:28,], stan_case_model) #FOI: 0|0
KOR_model_01_05 <- model_run_case_year_sum(Kor_inci_data[29:35,], stan_case_model) #FOI: 0.020|0.026 #rho: 0.022
KOR_model_01_06 <- model_run_case_year_sum(Kor_inci_data[36:42,], stan_case_model) #FOI: 0|0
KOR_model_01_07 <- model_run_case_year_sum(Kor_inci_data[43:49,], stan_case_model) #FOI: 0.005|0.015 #rho: 0.027
KOR_model_01_08 <- model_run_case_year_sum(Kor_inci_data[50:56,], stan_case_model) #FOI: 0.005|0.014 #rho: 0.04
KOR_model_01_09 <- model_run_case_year_sum(Kor_inci_data[57:63,], stan_case_model) #FOI: 0.021|0.043 #rho: 0.003
KOR_model_01_10 <- model_run_case_year_sum(Kor_inci_data[64:70,], stan_case_model) #FOI: 0.008|0.025 #rho: 0.003
KOR_model_01_11 <- model_run_case_year_sum(Kor_inci_data[71:77,], stan_case_model) #FOI: 0.014|0.040 #rho: 0.014
KOR_model_01_12 <- model_run_case_year_sum(Kor_inci_data[78:84,], stan_case_model) #FOI: 0.012|0.032 => the best result but not totally fit  #rho: 0.0001
KOR_model_01_13 <- model_run_case_year_sum(Kor_inci_data[85:91,], stan_case_model) #FOI: 0.002|0.006 #rho: 0.050
KOR_model_01_14 <- model_run_case_year_sum(Kor_inci_data[92:98,], stan_case_model) #FOI: 0.003|0.019 #rho: 0.010
KOR_model_01_15 <- model_run_case_year_sum(Kor_inci_data[99:105,], stan_case_model) #all year combined #FOI: |0.029 => fit well  #rho: 0.00005
KOR_model_01_sameFOIrho <- model_run_case_year_sum(Kor_inci_data[1:98,], stan_case_model) #all year have the same rho, same FOI: 0.0205
KOR_model_01_samerho <- model_run_case_year_sum_cohort_same_rho(Kor_inci_data[1:98,], stan_case_cohort_same_rho_model) #all year have the same rho: not fit well
KOR_model_01_samerho_e <- model_run_case_year_sum_cohort_same_rho(Kor_inci_data[c(1:14, 29:35, 43:98),], stan_case_cohort_same_rho_model) #all year have the same rho, eliminate years that have 0 case: fit better, FOI range from 0.001 to 0.06 (most likely 0.02)
KOR_model_01_sameFOI <- model_run_case_year_sum_cohort_same_rho(Kor_inci_data[1:98,], stan_case_cohort_same_FOI_model)# FOI: 0.028 
#All of result not fit well => how to include the information of vaccination
#All better fit but not fit well after including vaccination information.

KOR_model_02 <- model_run_case_year_sum(Kor_inci_data[106:122,], stan_case_model) #FOI: 0.0063 => not fit well
KOR_model_03 <- model_run_case_year_sum(Kor_inci_data[123:130,], stan_case_model) #FOI: 0.0098 => not fit quite well

#Indonesia:
Indo_inci_data = read.xlsx(paste(data_inci_path, "Indonesia_JE.xlsx", sep = ""), 1)
IDN_model_01 <- model_run_case_year_sum(Indo_inci_data[1:3,], stan_case_model) #FOI: 0.126 => not fit well|FOI: 0.19 => fit well
IDN_model_02 <- model_run_case_year_sum(Indo_inci_data[4:7,], stan_case_model) #FOI: 0.39|FOI: 0.35|FOI: 0.37

#India:
India_inci_data = read.xlsx(paste(data_inci_path, "India_JE.xlsx", sep = ""), 1)

IND_model_01_01 <- model_run_case_year_sum(India_inci_data[1:4,], stan_case_model) #FOI: 0.001 => not fit well, after vaccine
IND_model_01_02 <- model_run_case_year_sum(India_inci_data[5:9,], stan_case_model) #FOI: 0.06 => before vaccine
IND_model_01_sameFOIrho <- model_run_case_year_sum(India_inci_data[1:9,], stan_case_model) #FOI: 0.016 => fit better
IND_model_01_samerho <- model_run_case_year_sum_cohort_same_rho(India_inci_data[1:9,], stan_case_cohort_same_rho_model) #FOI 1: 0.019, FOI 2: 0.14 => fit well
IND_model_01_sameFOI <- model_run_case_year_sum_cohort_same_rho(India_inci_data[1:9,], stan_case_cohort_same_FOI_model) #FOI: 0.09

IND_model_02 <- model_run_case_year_sum(India_inci_data[10:13,], stan_case_model) #FOI: 0.002 => not fit well| #FOI: 0.002 not fit well

IND_model_03_01 <- model_run_case_year_sum(India_inci_data[14:18,], stan_case_model) #FOI: 0.51
IND_model_03_02 <- model_run_case_year_sum(India_inci_data[19:23,], stan_case_model) #FOI: 0.207
IND_model_03_03 <- model_run_case_year_sum(India_inci_data[24:28,], stan_case_model) #FOI: 0.187
IND_model_03_all <- model_run_case_year_sum(India_inci_data[29:33,], stan_case_model) #FOI: 0.21
IND_model_03_sameFOIrho <- model_run_case_year_sum(India_inci_data[14:28,], stan_case_model) #FOI: 0.21
IND_model_03_samerho <- model_run_case_year_sum_cohort_same_rho(India_inci_data[14:28,], stan_case_cohort_same_rho_model)#FOI: 0.5, 0.2, 0.18
IND_model_03_sameFOI <- model_run_case_year_sum_cohort_same_rho(India_inci_data[14:28,], stan_case_cohort_same_FOI_model) #FOI: 0.21

IND_model_04_01 <- model_run_case_year_sum(India_inci_data[34:39,], stan_case_model) #FOI: 0.101
IND_model_04_02 <- model_run_case_year_sum(India_inci_data[40:45,], stan_case_model) #FOI: 0.062
IND_model_04_03 <- model_run_case_year_sum(India_inci_data[46:51,], stan_case_model) #FOI: 0.0.051
IND_model_04_all <- model_run_case_year_sum(India_inci_data[52:57,], stan_case_model) #FOI: 0.068|FOI:0.063|FOI: 0.063
IND_model_04_sameFOIrho <- model_run_case_year_sum(India_inci_data[34:51,], stan_case_model) #FOI: 0.063
IND_model_04_samerho <- model_run_case_year_sum_cohort_same_rho(India_inci_data[34:51,], stan_case_cohort_same_rho_model)#FOI: 0.101 0.063 0.048 => not fit well
IND_model_04_sameFOI <- model_run_case_year_sum_cohort_same_rho(India_inci_data[34:51,], stan_case_cohort_same_FOI_model) #FOI: 0.064

IND_model_05 <- model_run_case_year_sum(India_inci_data[58:60,], stan_case_model) #FOI: 0.014 => not fit well|#FOI: 0.014 => not fit well
IND_model_06 <- model_run_case_year_sum(India_inci_data[61:65,], stan_case_model) #FOI: 0.0009 => not fit well (included vaccination program)

IND_model_07_01 <- model_run_case_year_sum(India_inci_data[66:68,], stan_case_model) #FOI: 0.161
IND_model_07_02 <- model_run_case_year_sum(India_inci_data[69:71,], stan_case_model) #FOI: 0.155
IND_model_07_03 <- model_run_case_year_sum(India_inci_data[72:74,], stan_case_model) #FOI: 0.153
IND_model_07_all <- model_run_case_year_sum(India_inci_data[75:77,], stan_case_model) #FOI: 0.157
IND_model_07_sameFOIrho <- model_run_case_year_sum(India_inci_data[66:74,], stan_case_model) #FOI: 0.157
IND_model_07_samerho <- model_run_case_year_sum_cohort_same_rho(India_inci_data[66:74,], stan_case_cohort_same_rho_model)#FOI: 0.140, 0.197, 0.155
IND_model_07_sameFOI <- model_run_case_year_sum_cohort_same_rho(India_inci_data[66:74,], stan_case_cohort_same_FOI_model)#FOI: 0.16

IND_model_08 <- model_run_case_year_sum(India_inci_data[78:92,], stan_case_model) #FOI: 0.002 => not fit well (included vaccination program)
IND_model_09 <- model_run_case_year_sum(India_inci_data[93:99,], stan_case_model) #FOI: 0.170
IND_model_10 <- model_run_case_year_sum(India_inci_data[100:104,], stan_case_model) #FOI: 0.148
IND_model_11 <- model_run_case_year_sum(India_inci_data[105:112,], stan_case_model) #FOI: 0.054 => not fit well (included vaccination program)
IND_model_12 <- model_run_case_year_sum(India_inci_data[113:128,], stan_case_model) #FOI: 0.034 => not fit well
IND_model_13 <- model_run_case_year_sum(India_inci_data[129:132,], stan_case_model) #FOI: 0.179
IND_model_14 <- model_run_case_year_sum(India_inci_data[133:149,], stan_case_model) #FOI: 0.276
IND_model_15 <- model_run_case_year_sum(India_inci_data[150:154,], stan_case_model) #FOI: 0.051 => not fit quite well
IND_model_16 <- model_run_case_year_sum(India_inci_data[155:158,], stan_case_model) #FOI: 0.137 => not fit quite well

#Nepal:
Nepal_inci_data = read.xlsx(paste(data_inci_path, "Nepal_JE.xlsx", sep = ""), 1)
NPL_model_01 <- model_run_case_year_sum(Nepal_inci_data[1:3,], stan_case_model) #FOI: 0.104|FOI: 0.084|FOI: 0.084 => Kathmandu valley is supposed to be higher than non-Kathmandu.
NPL_model_02 <- model_run_case_year_sum(Nepal_inci_data[4:6,], stan_case_model) #FOI: NA|FOI: 0.074|FOI: 0.074
NPL_model_12_samerho <- model_run_case_year_sum_cohort_same_rho(Nepal_inci_data[1:6,], stan_case_cohort_same_rho_model) #FOI: 0.15 0.59

NPL_model_03 <- model_run_case_year_sum(Nepal_inci_data[7:12,], stan_case_model) #FOI: 0.062|FOI: 0.0675|FOI: 0.074 => Western Terai is supposed to be higher than Non-west.
NPL_model_04 <- model_run_case_year_sum(Nepal_inci_data[13:18,], stan_case_model) #FOI: 0.072|FOI: 0.072|FOI: 0.070
NPL_model_34_samerho <- model_run_case_year_sum_cohort_same_rho(Nepal_inci_data[7:18,], stan_case_cohort_same_rho_model)#FOI: 0.05 0.98

NPL_model_05 <- model_run_case_year_sum(Nepal_inci_data[19:24,], stan_case_model) #FOI: 0.017 => during outbreak|FOI:0.023|FOI: 0.024
NPL_model_06 <- model_run_case_year_sum(Nepal_inci_data[25:41,], stan_case_model) #FOI: 0.016 => not fit well|FOI:0.00079 => not fit well|FOI: 0.0008 not fit well
NPL_model_07 <- model_run_case_year_sum(Nepal_inci_data[42:48,], stan_case_model) #FOI: 0.114|FOI: 0.12|FOI: 0.12
NPL_model_08 <- model_run_case_year_sum(Nepal_inci_data[49:51,], stan_case_model) #FOI: 0.075 => combined 1st and 2nd model|FOI: 0.075
NPL_model_09 <- model_run_case_year_sum(Nepal_inci_data[52:57,], stan_case_model) #FOI: 0.058 => combined 3rd and 4th model|FOI: 0.058
NPL_model_10 <- model_run_case_year_sum(Nepal_inci_data[58:60,], stan_case_model) #FOI: 0.010 => not fit well

NPL_model_11_01 <- model_run_case_year_sum(Nepal_inci_data[61:64,], stan_case_model) #FOI: 0.11
NPL_model_11_02 <- model_run_case_year_sum(Nepal_inci_data[65:68,], stan_case_model) #FOI: 0.12
NPL_model_11_03 <- model_run_case_year_sum(Nepal_inci_data[69:72,], stan_case_model) #FOI: 0.13
NPL_model_11_04 <- model_run_case_year_sum(Nepal_inci_data[73:76,], stan_case_model) #FOI: 0.18
NPL_model_11_all <- model_run_case_year_sum(Nepal_inci_data[77:80,], stan_case_model) #FOI: 0.136
NPL_model_sameFOIrho <- model_run_case_year_sum(Nepal_inci_data[61:76,], stan_case_model) #FOI: 0.136
NPL_model_samerho <- model_run_case_year_sum_cohort_same_rho(Nepal_inci_data[61:76,], stan_case_cohort_same_rho_model)#FOI: 0.13 0.11 0.12 0.16
NPL_model_sameFOI <- model_run_case_year_sum_cohort_same_rho(Nepal_inci_data[61:76,], stan_case_cohort_same_FOI_model)#FOI: 0.136

NPL_model_12 <- model_run_case_year_sum(Nepal_inci_data[81:108,], stan_case_model) #FOI: 0.014 => not fit quite well
NPL_model_13 <- model_run_case_year_sum(Nepal_inci_data[109:111,], stan_case_model) #FOI: 0.106

#Taiwan:cohort data
Taiwan_inci_data = read.xlsx(paste(data_inci_path, "Taiwan_JE.xlsx", sep = ""), 1)
TWN_model_01_01 <- model_run_case_year_sum(Taiwan_inci_data[1:7,], stan_case_model) #FOI: 0.0007|0.084|0.084
TWN_model_01_02 <- model_run_case_year_sum(Taiwan_inci_data[8:14,], stan_case_model) #FOI: 0.0014|0.092|0.094
TWN_model_01_03 <- model_run_case_year_sum(Taiwan_inci_data[15:21,], stan_case_model) #FOI: 0.0003|0.065|0.061
TWN_model_01_04 <- model_run_case_year_sum(Taiwan_inci_data[22:28,], stan_case_model) #FOI: 0.0003|0.062|0.057
TWN_model_01_05 <- model_run_case_year_sum(Taiwan_inci_data[29:35,], stan_case_model) #FOI: 0.0003|0.0567|0.052
TWN_model_01_06 <- model_run_case_year_sum(Taiwan_inci_data[36:42,], stan_case_model) #FOI: 0.0003|0.060|0.058
TWN_model_01_07 <- model_run_case_year_sum(Taiwan_inci_data[43:49,], stan_case_model) #FOI: 0.0022|0.060|0.081
TWN_model_01_08 <- model_run_case_year_sum(Taiwan_inci_data[50:56,], stan_case_model) #FOI: 0.0003|0.055|0.053
TWN_model_01_09 <- model_run_case_year_sum(Taiwan_inci_data[57:63,], stan_case_model) #FOI: 0.0003|0.062|0.062
TWN_model_01_10 <- model_run_case_year_sum(Taiwan_inci_data[64:70,], stan_case_model) #FOI: 0.0004|0.063|0.0625
TWN_model_01_11 <- model_run_case_year_sum(Taiwan_inci_data[71:77,], stan_case_model) #FOI: 0.0003|0.062|0.062
TWN_model_01_all <- model_run_case_year_sum(Taiwan_inci_data[78:84,], stan_case_model) #FOI: 0.0612
TWN_model_01_sameFOIrho <- model_run_case_year_sum(Taiwan_inci_data[1:77,], stan_case_model) #FOI: 0.061
TWN_model_01_samerho <- model_run_case_year_sum_cohort_same_rho(Taiwan_inci_data[1:77,], stan_case_cohort_same_rho_model)#FOI: 0.080 0.078 0.06 0.055 0.058 0.053 0.08 0.067 0.056 0.062 0.055
TWN_model_01_sameFOI <- model_run_case_year_sum_cohort_same_rho(Taiwan_inci_data[1:77,], stan_case_cohort_same_FOI_model) #FOI: 0.063

TWN_model_02_01 <- model_run_case_year_sum(Taiwan_inci_data[85:88 + 0*4,], stan_case_model) #FOI: 0.062 => not fit quite well
TWN_model_02_02 <- model_run_case_year_sum(Taiwan_inci_data[85:88 + 1*4,], stan_case_model) #FOI: 0.076
TWN_model_02_03 <- model_run_case_year_sum(Taiwan_inci_data[85:88 + 2*4,], stan_case_model) #FOI: 0.070
TWN_model_02_04 <- model_run_case_year_sum(Taiwan_inci_data[85:88 + 3*4,], stan_case_model) #FOI: 0.085
TWN_model_02_05 <- model_run_case_year_sum(Taiwan_inci_data[85:88 + 4*4,], stan_case_model) #FOI: 0.063
TWN_model_02_06 <- model_run_case_year_sum(Taiwan_inci_data[85:88 + 5*4,], stan_case_model) #FOI: 0.060
TWN_model_02_07 <- model_run_case_year_sum(Taiwan_inci_data[85:88 + 6*4,], stan_case_model) #FOI: 0.054
TWN_model_02_08 <- model_run_case_year_sum(Taiwan_inci_data[85:88 + 7*4,], stan_case_model) #FOI: 0.056
TWN_model_02_09 <- model_run_case_year_sum(Taiwan_inci_data[85:88 + 8*4,], stan_case_model) #FOI: 0.042
TWN_model_02_10 <- model_run_case_year_sum(Taiwan_inci_data[85:88 + 9*4,], stan_case_model) #FOI: 0.042
TWN_model_02_11 <- model_run_case_year_sum(Taiwan_inci_data[85:88 + 10*4,], stan_case_model) #FOI: 0.057
TWN_model_02_12 <- model_run_case_year_sum(Taiwan_inci_data[85:88 + 11*4,], stan_case_model) #FOI: 0.061
TWN_model_02_13 <- model_run_case_year_sum(Taiwan_inci_data[85:88 + 12*4,], stan_case_model) #FOI: 0.053
TWN_model_02_14 <- model_run_case_year_sum(Taiwan_inci_data[85:88 + 13*4,], stan_case_model) #FOI: 0.067
TWN_model_02_15 <- model_run_case_year_sum(Taiwan_inci_data[85:88 + 14*4,], stan_case_model) #FOI: 0.045
TWN_model_02_all <- model_run_case_year_sum(Taiwan_inci_data[145:151,], stan_case_model) #FOI: 0.068
TWN_model_02_sameFOIrho <- model_run_case_year_sum(Taiwan_inci_data[85:144,], stan_case_model) #FOI: 0.058
TWN_model_02_samerho <- model_run_case_year_sum_cohort_same_rho(Taiwan_inci_data[85:144,], stan_case_cohort_same_rho_model)#FOI: some are not well fit
TWN_model_02_sameFOI <- model_run_case_year_sum_cohort_same_rho(Taiwan_inci_data[85:144,], stan_case_cohort_same_FOI_model) #FOI: 0.061

TWN_model_03_01 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 0*5,], stan_case_model) #FOI: 0.178
TWN_model_03_02 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 1*5,], stan_case_model) #FOI: 0.182
TWN_model_03_03 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 2*5,], stan_case_model) #FOI: 0.280
TWN_model_03_04 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 3*5,], stan_case_model) #FOI: 0.240
TWN_model_03_05 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 4*5,], stan_case_model) #FOI: 0.255
TWN_model_03_06 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 5*5,], stan_case_model) #FOI: 0.333
TWN_model_03_07 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 6*5,], stan_case_model) #FOI: 0.322
TWN_model_03_08 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 7*5,], stan_case_model) #FOI: 0.321
TWN_model_03_09 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 8*5,], stan_case_model) #FOI: 0.433
TWN_model_03_10 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 9*5,], stan_case_model) #FOI: 0.333
TWN_model_03_11 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 10*5,], stan_case_model) #FOI: 0.333
TWN_model_03_12 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 11*5,], stan_case_model) #FOI: 0.306 => not fit quite well
TWN_model_03_13 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 12*5,], stan_case_model) #FOI: 0.277
TWN_model_03_14 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 13*5,], stan_case_model) #FOI: 0.275
TWN_model_03_15 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 14*5,], stan_case_model) #FOI: 0.293
TWN_model_03_16 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 15*5,], stan_case_model) #FOI: 0.279
TWN_model_03_17 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 16*5,], stan_case_model) #FOI: 0.331
TWN_model_03_18 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 17*5,], stan_case_model) #FOI: 0.370
TWN_model_03_19 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 18*5,], stan_case_model) #FOI: 0.358
TWN_model_03_20 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 19*5,], stan_case_model) #FOI: 0.311
TWN_model_03_21 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 20*5,], stan_case_model) #FOI: 0.312
TWN_model_03_22 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 21*5,], stan_case_model) #FOI: 0.287
TWN_model_03_23 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 22*5,], stan_case_model) #FOI: 0.240
TWN_model_03_24 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 23*5,], stan_case_model) #FOI: 0.316
TWN_model_03_25 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 24*5,], stan_case_model) #FOI: 0.203
TWN_model_03_26 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 25*5,], stan_case_model) #FOI: 0.196 => not fit quite well
TWN_model_03_27 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 26*5,], stan_case_model) #FOI: 0.246 => not fit quite well
TWN_model_03_28 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 27*5,], stan_case_model) #FOI: not fit
TWN_model_03_29 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 28*5,], stan_case_model) #FOI: not fit
TWN_model_03_30 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 29*5,], stan_case_model) #FOI: 0.113 => not fit quite well
TWN_model_03_31 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 30*5,], stan_case_model) #FOI: 0.152 => not fit quite well
TWN_model_03_32 <- model_run_case_year_sum(Taiwan_inci_data[152:156 + 31*5,], stan_case_model) #FOI: 0.114 => not fit quite well
TWN_model_03_all <- model_run_case_year_sum(Taiwan_inci_data[312:316,], stan_case_model) #FOI: 0.199
TWN_model_03_sameFOIrho <- model_run_case_year_sum(Taiwan_inci_data[152:311,], stan_case_model) #FOI: 0.201
TWN_model_03_samerho <- model_run_case_year_sum_cohort_same_rho(Taiwan_inci_data[152:311,], stan_case_cohort_same_rho_model)#FOI: some are not well fit
TWN_model_03_samerho_e <- model_run_case_year_sum_cohort_same_rho(Taiwan_inci_data[c(152:286,297:311),], stan_case_cohort_same_rho_model)#FOI: eliminate some year that not fit
TWN_model_03_sameFOI <- model_run_case_year_sum_cohort_same_rho(Taiwan_inci_data[152:311,], stan_case_cohort_same_FOI_model) #FOI: 0.249
#All of result not fit well => how to include the information of vaccination
#already include the vaccination programe

#Vietnam:
Vietnam_inci_data = read.xlsx(paste(data_inci_path, "Vietnam_JE.xlsx", sep = ""), 1)
VNM_model_01 <- model_run_case_year_sum(Vietnam_inci_data[1:7,], stan_case_model) #FOI: 0.26|FOI: 0.23
VNM_model_02 <- model_run_case_year_sum(Vietnam_inci_data[8:12,], stan_case_model) #FOI: 0.23:FOI: 0.20

#Thailand:
Thailand_inci_data = read.xlsx(paste(data_inci_path, "Thailand_JE.xlsx", sep = ""), 1)
THL_model_01 <- model_run_case_year_sum(Thailand_inci_data[1:13,], stan_case_model) #FOI: 0.054|FOI: 0.056|FOI: 0.056

#Bangladesh:
Bangladesh_inci_data = read.xlsx(paste(data_inci_path, "Bangladesh_JE.xlsx", sep = ""), 1)
BGD_model_01 <- model_run_case_year_sum(Bangladesh_inci_data[1:6,], stan_case_model) #FOI: 0.064|FOI: 0.062|FOI: 0.062

#China:
China_inci_data = read.xlsx(paste(data_inci_path, "China_JE.xlsx", sep = ""), 1)

CHN_model_01_01 <- model_run_case_year_sum(China_inci_data[1:3,], stan_case_model) #FOI: 0.404
CHN_model_01_02 <- model_run_case_year_sum(China_inci_data[4:6,], stan_case_model) #FOI: 0.441
CHN_model_01_03 <- model_run_case_year_sum(China_inci_data[7:9,], stan_case_model) #FOI: 0.453
CHN_model_01_04 <- model_run_case_year_sum(China_inci_data[10:12,], stan_case_model) #FOI: 0.470
CHN_model_01_05 <- model_run_case_year_sum(China_inci_data[13:15,], stan_case_model) #FOI: 1.25
CHN_model_01_06 <- model_run_case_year_sum(China_inci_data[16:18,], stan_case_model) #FOI: 0.85
CHN_model_01_all <- model_run_case_year_sum(China_inci_data[19:21,], stan_case_model) #FOI: 0.40
CHN_model_01_sameFOIrho <- model_run_case_year_sum(China_inci_data[1:18,], stan_case_model) #FOI: 0.40
CHN_model_01_samerho <- model_run_case_year_sum_cohort_same_rho(China_inci_data[1:18,], stan_case_cohort_same_rho_model) #FOI: 0.42 0.46 0.45 0.36 0.42 0.35
CHN_model_01_sameFOI <- model_run_case_year_sum_cohort_same_rho(China_inci_data[1:18,], stan_case_cohort_same_FOI_model) #FOI: 0.64

CHN_model_02 <- model_run_case_year_sum(China_inci_data[22:39,], stan_case_model) #FOI: 0.288

CHN_model_03_01 <- model_run_case_year_sum(China_inci_data[40:45,], stan_case_model) #FOI: 0.063
CHN_model_03_02 <- model_run_case_year_sum(China_inci_data[46:51,], stan_case_model) #FOI: 0.112
CHN_model_03_03 <- model_run_case_year_sum(China_inci_data[52:57,], stan_case_model) #FOI: 0.117
CHN_model_03_04 <- model_run_case_year_sum(China_inci_data[58:63,], stan_case_model) #FOI: 0.49
CHN_model_03_samerho <- model_run_case_year_sum_cohort_same_rho(China_inci_data[40:63,], stan_case_cohort_same_rho_model) #FOI: 0.072 0.176 0.31 0.30

CHN_model_04_01 <- model_run_case_year_sum(China_inci_data[64:67,], stan_case_model) #FOI: 0.158
CHN_model_04_02 <- model_run_case_year_sum(China_inci_data[68:71,], stan_case_model) #FOI: 0.182
CHN_model_04_03 <- model_run_case_year_sum(China_inci_data[72:75,], stan_case_model) #FOI: 0.184
CHN_model_04_04 <- model_run_case_year_sum(China_inci_data[76:79,], stan_case_model) #FOI: 0.174
CHN_model_04_05 <- model_run_case_year_sum(China_inci_data[80:83,], stan_case_model) #FOI: 0.189
CHN_model_04_06 <- model_run_case_year_sum(China_inci_data[84:87,], stan_case_model) #FOI: 0.191
CHN_model_04_07 <- model_run_case_year_sum(China_inci_data[90:93,], stan_case_model) #FOI: 0.181
CHN_model_04_08 <- model_run_case_year_sum(China_inci_data[96:99,], stan_case_model) #FOI: 0.184
CHN_model_04_09 <- model_run_case_year_sum(China_inci_data[100:103,], stan_case_model) #FOI: 0.183
CHN_model_04_all <- model_run_case_year_sum(China_inci_data[104:107,], stan_case_model) #FOI: 0.181
CHN_model_04_sameFOIrho <- model_run_case_year_sum(China_inci_data[64:103,], stan_case_model) #FOI: 0.181
CHN_model_04_samerho <- model_run_case_year_sum_cohort_same_rho(China_inci_data[64:103,], stan_case_cohort_same_rho_model) #FOI: 0.16 0.18 0.17 0.18 0.18 0.19 0.18 0.19 0.2 0.23 0.20
CHN_model_04_sameFOI <- model_run_case_year_sum_cohort_same_rho(China_inci_data[64:103,], stan_case_cohort_same_FOI_model) #FOI: 0.181

CHN_model_05_01 <- model_run_case_year_sum(China_inci_data[108:115,], stan_case_model) #FOI: 0.126
CHN_model_05_02 <- model_run_case_year_sum(China_inci_data[116:123,], stan_case_model) #FOI: 0.036
CHN_model_05_03 <- model_run_case_year_sum(China_inci_data[124:131,], stan_case_model) #FOI: 0.040
CHN_model_05_04 <- model_run_case_year_sum(China_inci_data[132:139,], stan_case_model) #FOI: 0.024 => fit not well
CHN_model_05_05 <- model_run_case_year_sum(China_inci_data[140:147,], stan_case_model) #FOI: 0.068
CHN_model_05_06 <- model_run_case_year_sum(China_inci_data[148:155,], stan_case_model) #FOI: 0.040 => fit not well
CHN_model_05_07 <- model_run_case_year_sum(China_inci_data[156:163,], stan_case_model) #FOI: 0.054 => fit not well
CHN_model_05_all <- model_run_case_year_sum(China_inci_data[164:171,], stan_case_model) #FOI: 0.046
CHN_model_05_sameFOIrho <- model_run_case_year_sum(China_inci_data[108:163,], stan_case_model) #FOI: 0.046
CHN_model_05_samerho <- model_run_case_year_sum_cohort_same_rho(China_inci_data[108:163,], stan_case_cohort_same_rho_model) #FOI: 0.13 0.04 0.04 0.03 0.07 0.02 0.02
CHN_model_05_sameFOI <- model_run_case_year_sum_cohort_same_rho(China_inci_data[108:163,], stan_case_cohort_same_FOI_model) #FOI: 0.044

CHN_model_06 <- model_run_case_year_sum(China_inci_data[172:180,], stan_case_model) #FOI: 0.64
CHN_model_07 <- model_run_case_year_sum(China_inci_data[181:183,], stan_case_model) #FOI: 0.082
CHN_model_08 <- model_run_case_year_sum(China_inci_data[184:191,], stan_case_model) #FOI: 0.747
CHN_model_09 <- model_run_case_year_sum(China_inci_data[192:194,], stan_case_model) #FOI: 0.040
CHN_model_10 <- model_run_case_year_sum(China_inci_data[195:204,], stan_case_model) #FOI: 0.943

CHN_model_11_first_half <- model_run_case_year_sum(China_inci_data[205:207,], stan_case_model) #FOI: 0.143
CHN_model_11_second_half <- model_run_case_year_sum(China_inci_data[208:210,], stan_case_model) #FOI: 0.121
CHN_model_11_all <- model_run_case_year_sum(China_inci_data[211:213,], stan_case_model) #FOI: 0.135
CHN_model_11_01 <- model_run_case_year_sum(China_inci_data[214:216,], stan_case_model) #FOI: 0.142
CHN_model_11_02 <- model_run_case_year_sum(China_inci_data[217:219,], stan_case_model) #FOI: 0.153
CHN_model_11_03 <- model_run_case_year_sum(China_inci_data[220:222,], stan_case_model) #FOI: 0.128
CHN_model_11_04 <- model_run_case_year_sum(China_inci_data[223:225,], stan_case_model) #FOI: 0.137
CHN_model_11_05 <- model_run_case_year_sum(China_inci_data[226:228,], stan_case_model) #FOI: 0.146
CHN_model_11_06 <- model_run_case_year_sum(China_inci_data[229:231,], stan_case_model) #FOI: 0.135
CHN_model_11_07 <- model_run_case_year_sum(China_inci_data[232:234,], stan_case_model) #FOI: 0.114
CHN_model_11_08 <- model_run_case_year_sum(China_inci_data[235:237,], stan_case_model) #FOI: 0.155
CHN_model_11_09 <- model_run_case_year_sum(China_inci_data[238:240,], stan_case_model) #FOI: 0.127
CHN_model_11_10 <- model_run_case_year_sum(China_inci_data[241:243,], stan_case_model) #FOI: 0.063
CHN_model_11_11 <- model_run_case_year_sum(China_inci_data[244:246,], stan_case_model) #FOI: 0.104
CHN_model_11_sameFOIrho <- model_run_case_year_sum(China_inci_data[214:246,], stan_case_model) #FOI: 0.129
CHN_model_11_samerho <- model_run_case_year_sum_cohort_same_rho(China_inci_data[214:246,], stan_case_cohort_same_rho_model) #FOI: 0.13 0.14 0.11 0.13 0.15 0.13 0.13 0.17 0.14 0.08 0.15
CHN_model_11_sameFOI <- model_run_case_year_sum_cohort_same_rho(China_inci_data[214:246,], stan_case_cohort_same_FOI_model) #FOI: 0.126

#Philippines:
Philippines_inci_data = read.xlsx(paste(data_inci_path, "Philippines_JE.xlsx", sep = ""), 1)
PHL_model_01 <- model_run_case_year_sum(Philippines_inci_data[1:6,], stan_case_model) #FOI: 0.17
PHL_model_02 <- model_run_case_year_sum(Philippines_inci_data[7:10,], stan_case_model) #FOI: 0.098 => pretty old study
PHL_model_03 <- model_run_case_year_sum(Philippines_inci_data[11:18,], stan_case_model) #FOI: 0.043 => this study have some results that are cross-reactive with DENV
 
#Get the estimated FOI from different models:
sel_countries_list = c("JPN", "LAO", "KHM", "KOR", "IND", "IDN", "NPL", "TWN", "VNM", "THL", "BGD", "CHN", "PHL")
FOI_summary_table = c()
Rho_summary_table = c()
for(coun in sel_countries_list){
  coun_model = mget(apropos(paste(coun, "_model", sep = "")))
  summary_coun_model = lapply(coun_model, function(x){
    y <- x$summary %>% data.frame %>% select(mean, n_eff, Rhat)
    y <- list(lambda = y[agrep("lambda",rownames(y)),], rho = y[agrep("rho",rownames(y)),])
    if(nrow(x$model_info) == 1){
      edited_model_info <- mutate(x$model_info, Index = as.factor(Index))
      y <- Reduce(rbind,y) %>% cbind(edited_model_info)
    } else {
      model_info_head_tail_comp <- data.frame(Index = do.call(paste, c(as.list(unique(x$model_info$Index)), sep = "_")), 
                                              Year = paste(substring(x$model_info$Year[1],1,7), substring(tail(x$model_info$Year,1), first = 9), sep = "_"),
                                              subnation = paste(x$model_info$subnation[1], tail(x$model_info$subnation, 1), sep = "_"), 
                                              ISO = do.call(paste, c(as.list(unique(x$model_info$ISO)), sep = "_")))
      for(par in 1:length(y)){
        if(nrow(y[[par]]) == 1) y[[par]] = cbind(y[[par]], model_info_head_tail_comp)
        else y[[par]] = cbind(y[[par]], x$model_info)
      }
      y <- Reduce(rbind,y)
    }
    return(y)
  })
  for(i in 1:length(coun_model)){
    summary_coun_model[[i]]$country_model = names(coun_model)[i]
  }
  summary_coun_model = Reduce(rbind, summary_coun_model)
  para_names = row.names(summary_coun_model)
  summary_coun_model = mutate(summary_coun_model, N_eff_ratio = n_eff/24000)
  summary_coun_model$para_names = para_names
  
  sel_col = c("country_model", "Index", "Year", "subnation", "ISO", "mean", "Rhat", "N_eff_ratio")
  store_lambda_coun = summary_coun_model[agrep("lambda",summary_coun_model$para_names), sel_col]
  rownames(store_lambda_coun) <- NULL
  store_rho_coun = summary_coun_model[agrep("rho",summary_coun_model$para_names), sel_col]
  rownames(store_rho_coun) <- NULL
  FOI_summary_table = c(FOI_summary_table, list(store_lambda_coun))
  Rho_summary_table = c(Rho_summary_table, list(store_rho_coun))
}
FOI_summary_table = Reduce(rbind,FOI_summary_table)
Rho_summary_table = Reduce(rbind,Rho_summary_table)

save_FOI_table_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/result_model/first_model/"
write.csv(FOI_summary_table, paste(save_FOI_table_path, "FOI_estimation_all_models.csv"), row.names = F)
write.csv(Rho_summary_table, paste(save_FOI_table_path, "Rho_estimation_all_models.csv"), row.names = F)
#######################
######model with different admission rate for children and adult.
#######################
ll_case_je_2="

data {

int S; //number of study
int S_begin_index[S]; //the begin index of the study
int S_end_index[S]; //the end index of the study
int total_row; //total number of row of the data frame

vector[total_row] age_l; //lower bound of the age group
vector[total_row] age_u; //upper bound of the age group
vector[total_row] case_age; //number of cases in each age group.
vector[total_row] pop_age; //population demo of the study year.
real t_case[S]; //total cases in every groups in each study

}

parameters {

real<lower = 0> lambda[S];
real<lower = 0> rho_report[S];
real<lower = 0> rho_young;
real<lower = 0> rho_old;
real<lower = 0> age_thres;

}

transformed parameters {

vector[total_row] i_age_dept; //the expected proportion of infections in each age group in each study
vector[total_row] e_age; //The expected number of cases in each age group in each study

  for(study in 1:S){
    for(i in (S_begin_index[study]):(S_end_index[study])){
      i_age_dept[i] = rho_young*(exp(-lambda[study]*fmin(age_l[i], age_thres)) - exp(-lambda[study]*(fmin(age_u[i], age_thres) + 1))) + rho_old*(exp(-lambda[study]*fmax(age_l[i], age_thres)) - exp(-lambda[study]*(fmax(age_u[i], age_thres) + 1)));
      e_age[i] = pop_age[i]*i_age_dept[i]*rho_report[study];
    }
  }

}

model {

real l_MN[S];
//prior distribution
lambda[S] ~ normal(0, 1000);
rho_report[S] ~ uniform(0, 1);
rho_young ~ uniform(0, 1);
rho_old ~ uniform(0, 1);
age_thres ~ uniform(0, max(age_u));

  for(study in 1:S){
    //MN likelihood function:
    l_MN[study] = lgamma(t_case[study] + 1) - sum(lgamma(case_age[(S_begin_index[study]):(S_end_index[study])] + 1)) + sum(case_age[(S_begin_index[study]):(S_end_index[study])] .* log(e_age[(S_begin_index[study]):(S_end_index[study])]/sum(e_age[(S_begin_index[study]):(S_end_index[study])]))); 
    //likelihood function, included poisson for total cases across all age group:
    target += l_MN[study] + t_case[study]*log(sum(e_age[(S_begin_index[study]):(S_end_index[study])])) - sum(e_age[(S_begin_index[study]):(S_end_index[study])]) - lgamma(t_case[study] + 1);
  }
}
"

stan_case_model = stan_model(model_code = ll_case_je_2)

data_inci_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_raw/cases_sero_data/"

#Gather all data:
sel_subnation_data <-function(l_country_inci_data, l_study_index){
  comb_study_data = data.frame(Year = character(), subnation = character(), ISO = character(), Age_group = character(),
                               Case_Sero = numeric(), Pop_all_age_year_sum = numeric())
  for(study in 1:length(l_country_inci_data)){
    country_data = l_country_inci_data[[study]]
    country_study_index = l_study_index[[study]]
    sel_study_index = unlist(unique(country_data$subnation))[country_study_index]
    sel_coutry_data = country_data[country_data$subnation %in% sel_study_index,]
    comb_study_data = rbind(comb_study_data, sel_coutry_data[,c("Year","subnation", "ISO", "Age_group", "Case_Sero", "Pop_all_age_year_sum")])
  }
  return(comb_study_data)
}

Laos_inci_data = read.xlsx(paste(data_inci_path, "Laos_JE.xlsx", sep = ""), 1) #1st study
Indo_inci_data = read.xlsx(paste(data_inci_path, "Indonesia_JE.xlsx", sep = ""), 1) # 1st, 2nd study
Cam_inci_data = read.xlsx(paste(data_inci_path, "Cambodia_JE.xlsx", sep = ""), 1) #2nd, 3rd, 4th study
Nepal_inci_data = read.xlsx(paste(data_inci_path, "Nepal_JE.xlsx", sep = ""), 1) #1st, 3rd, 4th, 5th, 6th, 7th
Vietnam_inci_data = read.xlsx(paste(data_inci_path, "Vietnam_JE.xlsx", sep = ""), 1) #1st, 2nd study
Thailand_inci_data = read.xlsx(paste(data_inci_path, "Thailand_JE.xlsx", sep = ""), 1) #1st study
Bangladesh_inci_data = read.xlsx(paste(data_inci_path, "Bangladesh_JE.xlsx", sep = ""), 1) #1st study
China_inci_data = read.xlsx(paste(data_inci_path, "China_JE.xlsx", sep = ""), 1) #6th -> 16th study

All_sel_coun_name = c("Laos", "Indo", "Cam", "Nepal", "Vietnam", "Thailand", "Bangladesh", "China")
All_sel_coun_name = c("China")
All_sel_countries_data = mget(paste(All_sel_coun_name, "_inci_data", sep = ""))
All_sel_coun_study_index = list(1, 1:2, 2:4, c(1, 3:5, 7), c(1,2), 1, 1, 6:15)
All_sel_coun_study_index = list(6:15)

All_sel_coun_study = sel_subnation_data(l_country_inci_data = All_sel_countries_data, l_study_index = All_sel_coun_study_index)

#Calculate the age_l and age_u:
age_group_split = strsplit(as.character(All_sel_coun_study$Age_group),"-")
age_l = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[1])}))
age_u = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[2])}))

all_study_length_count = table(as.character(All_sel_coun_study$subnation))
all_study_length_count = as.vector(all_study_length_count)[match(as.character(unique(All_sel_coun_study$subnation)), names(all_study_length_count))]
S_end_index = cumsum(all_study_length_count)
S_begin_index = c(1, head(S_end_index + 1, length(S_end_index) - 1))
t_case = c()
for(i in 1:length(S_begin_index)){
  t_case = c(t_case, sum(All_sel_coun_study$Case_Sero[S_begin_index[i]:S_end_index[i]]))
}

All_sel_coun_study_for_HMC = list(S = length(unique(All_sel_coun_study$subnation)), S_begin_index = S_begin_index, S_end_index = S_end_index, total_row = nrow(All_sel_coun_study), 
                                  age_l = age_l, age_u = age_u, case_age = All_sel_coun_study$Case_Sero, pop_age = All_sel_coun_study$Pop_all_age_year_sum, 
                                  t_case = t_case)

stan_FOI_fit = sampling(object = stan_case_model, data = All_sel_coun_study_for_HMC, 
                        chains = 4, iter = 8000)

#Posterior distribution of parameter:
mcmc_chain = as.array(stan_FOI_fit)
lambda_plot = mcmc_areas(mcmc_chain, regex_pars ="lambda")
rho_plot = mcmc_areas(mcmc_chain, regex_pars = "rho")
i_age_plot = mcmc_areas(mcmc_chain, regex_pars = "i_age") 
e_age_plot = mcmc_areas(mcmc_chain, regex_pars = "e_age")
com_post_plot = gridExtra::grid.arrange(lambda_plot, rho_plot, i_age_plot, e_age_plot, nrow = 2, ncol = 2)

