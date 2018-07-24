rm(list = ls())
#Test 3 - 16 countries, single model run for det model
#library loading:
library(rstan)
library(xlsx)
library(dplyr)

###data loading:
#Pop data load:#Montagu template load:
output_data_path = 'C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_clean/Montagu_data/test3/'
Det_template = read.csv(paste(output_data_path,'central_burden_template_JE-OUCRU-Clapham.csv', sep = ""))
Sto_template = read.csv(paste(output_data_path,'stochastic_burden_template_JE-OUCRU-Clapham.csv', sep = ""))
#cases data load:
data_inci_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_raw/cases_sero_data/"
#List of country need:
ISO_country_list = as.character(unique(Det_template$country))
#year select:
year_sel = 2000:2100
#age group select:
age_group_sel = 0:99
#Pop data load:
naive_pop_demo = read.xlsx2(paste(output_data_path,'naive_pop_1950_2100.xlsx', sep = ""), 1, colClasses=NA)
naive_pop_demo = naive_pop_demo[, c("country", "age_from", "age_to", paste("X", year_sel, sep = ""))]
#Expected life loss for each age in each country:
raw_expected_YLL_by_age = read.csv('C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_raw/Montagu-data/201710gavi-2_dds-201710_life_ex_both.csv')
raw_expected_YLL_by_age = raw_expected_YLL_by_age[raw_expected_YLL_by_age$year %in% year_sel, c("country_code", "age_from", "age_to", "year","value")]

###Reduced model which only generate lambda samples:
ll_case_je_reduced="

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
stan_sero_model_reduced = stan_model(model_code = ll_case_je_reduced)

#############################################

##########MODEL RUN FOR EACH COUNTRY

#############################################

#function to generate subnation population over time by the national population data and subnation population of a certain time:
#sub_pop_year, year: population of subnational region in year X....
#naive_pop_country: dataframe of population of nation over selected time and age group
sub_pop_gen <-function(sub_pop_year, year, naive_pop_country){
  naive_pop = naive_pop_country[,paste("X", year_sel, sep = "")]
  naive_pop_prop_year = colSums(naive_pop)/colSums(naive_pop)[year]
  naive_pop_sub_year = naive_pop_prop_year*sub_pop_year
  naive_pop_prop_age = sweep(naive_pop, 2, colSums(naive_pop), FUN = "/")
  naive_pop_sub_age_year = sweep(naive_pop_prop_age, 2, naive_pop_sub_year, FUN = "*")
  naive_pop_sub_age_year = cbind(naive_pop_country[,c("country", "age_from", "age_to")], naive_pop_sub_age_year)
  return(naive_pop_sub_age_year[naive_pop_country$age_from == age_group_sel,])
}

#get the function of generating campaign and routine population:
source("C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/code_JEV_model/JEV_Montagu_data - Scen_pop_generate_function.R")

#function to run model for each country => output is sampled of lambda:
model_run_each_country <- function(sel_country){
  print(sel_country)
  data_county_model = filter(All_country_model, country == sel_country)
  l_country_of_incidence_data = unlist(data_county_model$countries_inci_data)
  l_sel_subnation_index = unlist(data_county_model$index_sel_inci_data)
  
  
  #Get incidence data (literature review)
  l_data_inci_subnation = c()
  for(n_inci_data in 1:length(l_country_of_incidence_data)){
    data_inci = read.xlsx(paste(data_inci_path, l_country_of_incidence_data[n_inci_data], "_JE.xlsx", sep = ""), 1)
    subnation_inci = unlist(unique(data_inci$subnation))
    l_data_inci_subnation = c(l_data_inci_subnation,list(data_inci[data_inci$subnation == subnation_inci[l_sel_subnation_index[n_inci_data]],]))
  }
  
  #run model, get lambda for each data set:
  l_lambda_sampled = c()
  l_rho_sampled = c()
  l_data_sum = c()
  for(n_inci_data in 1:length(l_data_inci_subnation)){
    data = l_data_inci_subnation[[n_inci_data]]
    #Calculate the age_l and age_u:
    age_group_split = strsplit(as.character(data$Age_group),"-")
    age_l = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[1])}))
    age_u = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[2])}))
    
    data_for_HMC = list(N = nrow(data), age_l = age_l, age_u = age_u, 
                        case_age = data$Case_Sero, pop_age = data$Pop_all_age_year_sum, 
                        t_case = sum(data$Case_Sero))
    
    stan_FOI_fit = sampling(object = stan_sero_model_reduced, data = data_for_HMC, 
                            chains = 4, iter = 8000)
    lambda_all = extract(stan_FOI_fit, "lambda")$lambda
    rho_all = extract(stan_FOI_fit, "rho")$rho
    l_lambda_sampled = c(l_lambda_sampled, list(lambda_all))
    l_rho_sampled = c(l_rho_sampled, list(rho_all))
    l_data_sum = c(l_data_sum, list(summary(stan_FOI_fit)$summary))
  }
  return(list(lambda = l_lambda_sampled, rho = l_rho_sampled, sum_data = l_data_sum))
}

#function to generate cases, deaths, dalys: parameter = list(lambdas, sym_rate, death_rate, dis_rate)
#for Sto model, run can range from 1 to run_times, 
#for Det model, run = 1, lambdas is the mean, sym_rate = (1/500 + 1/250)/2.0, death_rate = 0.2, dis_rate = 0.4:
model_gen_each_country <- function(parameters, sel_country, run){
  print(sel_country)
  data_county_model = filter(All_country_model, country == sel_country)
  pop_gen <- filter(naive_pop_demo, country == sel_country)
  subnational_y_n = data_county_model$subnation
  sub_info_pop = unlist(data_county_model$subnation_pop_info_pop)
  sub_info_year = unlist(data_county_model$subnation_pop_info_year)
  ##data for generate quantities:
  #remained life expectancy:
  raw_expected_remain_life_by_country = raw_expected_YLL_by_age[raw_expected_YLL_by_age$country == sel_country,]
  infer_remain_life_1 = c() #infer the same value for each age group
  for(i in 1:nrow(raw_expected_remain_life_by_country)){
    if(i %% 22 == 1){
      infer_remain_life_1 = c(infer_remain_life_1, raw_expected_remain_life_by_country$value[i])
    } else if(i %% 22 == 2){
      infer_remain_life_1 = c(infer_remain_life_1, rep(raw_expected_remain_life_by_country$value[i], 4))
    } else if(i %% 22 != 0){
      infer_remain_life_1 = c(infer_remain_life_1, rep(raw_expected_remain_life_by_country$value[i], 5))
    }
  }
  year_rep_1 = rep(seq(year_sel[1],year_sel[length(year_sel)], 5), each = 100)
  infer_remain_life_2 = c() #infer the same value for each between years
  for(year in unique(raw_expected_remain_life_by_country$year)){
    year_data = infer_remain_life_1[year_rep_1 == year]
    infer_remain_life_2 = c(infer_remain_life_2, rep(year_data, 5))
  }
  infer_remain_life_2 = c(infer_remain_life_2, infer_remain_life_2[9901:10000]) #year 2100
  
  #infer subnational data from vaccination population:
  l_pop_gen_naive = list()
  l_pop_gen_cam = list()
  l_pop_gen_rou = list()
  #if vaccination campaign is only on subnational scale: get the local pop first, then vaccinate
  if(subnational_y_n){
    subnational_gen = sub_pop_gen(sub_info_pop, sub_info_year, pop_gen)
    
    l_pop_gen_naive = list(unlist(subnational_gen[,paste("X", year_sel, sep = "")]))
    l_pop_gen_cam = list(unlist(cam_scen_pop_gen(sel_country, subnational_gen, year_sel)[,paste("X", year_sel, sep = "")]))
    l_pop_gen_rou = list(unlist(rou_scen_pop_gen(sel_country, subnational_gen, year_sel)[,paste("X", year_sel, sep = "")]))
  } else { #if vaccination campaign is on national scale: get the pop vaccinate first, then get the local pop if there are many subregions
    
    pop_gen_naive = pop_gen
    pop_gen_cam = cam_scen_pop_gen(sel_country,pop_gen, year_sel)
    pop_gen_rou = rou_scen_pop_gen(sel_country,pop_gen, year_sel)
    
    if(length(sub_info_pop) != 1){
      for(sub in 1:length(sub_info_pop)){
        l_pop_gen_naive[[sub]] = unlist(pop_gen_naive[,paste("X", year_sel, sep = "")])*(sub_info_pop[sub]/sum(pop_gen[,sub_info_year[sub]]))
        l_pop_gen_cam[[sub]] = unlist(pop_gen_cam[,paste("X", year_sel, sep = "")])*(sub_info_pop[sub]/sum(pop_gen[,sub_info_year[sub]]))
        l_pop_gen_rou[[sub]] = unlist(pop_gen_rou[,paste("X", year_sel, sep = "")])*(sub_info_pop[sub]/sum(pop_gen[,sub_info_year[sub]]))
      }
    } else {
      l_pop_gen_naive = list(unlist(pop_gen_naive[,paste("X", year_sel, sep = "")]))
      l_pop_gen_cam = list(unlist(pop_gen_cam[,paste("X", year_sel, sep = "")]))
      l_pop_gen_rou = list(unlist(pop_gen_rou[,paste("X", year_sel, sep = "")]))
    }
  }
  l_pop_gen_scen <- list(naive = l_pop_gen_naive, cam = l_pop_gen_cam, rou = l_pop_gen_rou)
  
  #generate cases, deaths, DALYS from parameters:
  lambdas = parameters$lambdas #list of lambda
  sym_rate = parameters$sym_rate
  death_rate = parameters$death_rate
  dis_rate = parameters$dis_rate
  
  whole_country_output <- c()
  for(scen in 1:length(l_pop_gen_scen)){
    pop_scen <- l_pop_gen_scen[[scen]]
    l_data_output <- length(pop_scen[[1]])
    cases_scen_list <- list() #combine if there are multiple data
    
    cases_all_runs <- c()
    deaths_all_runs <- c()
    dalys_all_runs <- c()
    
    for(n_inci_data in 1:length(lambdas)){
      lambda_sampled = lambdas[[n_inci_data]][run]
      cases_scen_list[[n_inci_data]] = rep(sym_rate[run]*(1 - exp(-lambda_sampled)), l_data_output)*pop_scen[[n_inci_data]]* 
        exp(rep(age_group_sel, length(year_sel))*(-rep(lambda_sampled, l_data_output)))
    }
    cases_gen <- Reduce("+", cases_scen_list) %>% unname
    mort_gen = cases_gen*death_rate[run]
    disa_gen = (cases_gen - mort_gen)*dis_rate[run]
    DALYS_gen = mort_gen*infer_remain_life_2 + cases_gen*(0.133*2.5/52.0) + 
      disa_gen*0.542*infer_remain_life_2
    
    whole_country_output[[scen]] <- list(cases = cases_gen, deaths = mort_gen, dalys = DALYS_gen) #gen data for each scenario:
  }
  names(whole_country_output) <- c("naive", "cam", "rou")
  
  return(whole_country_output)
}

#####################all subnational population getting from Campbell 2011
####################
All_country_model = data.frame(country = c(ISO_country_list), 
                               subnation = c(F, T, F, T, F, F, F, F, F, F, T, F, F, F, F, F), 
                               countries_inci_data = I(list("Bangladesh", "Nepal", "Cambodia", "China", c(NULL, "Thailand", "Nepal"),
                                                       c("Thailand", "Indonesia"), "Thailand", "Laos", "Thailand", "Nepal", NA, "Thailand",
                                                       "Philippines", "Thailand", "Indonesia" ,"Laos")), 
                               index_sel_inci_data = I(list(1, 9, 3, 23, c(NULL, 1, 8), c(1, 1), 1, 1, 1, 9, NA, 1, 1, 1, 1, 1)),
                               subnation_pop_info_pop = I(list(NA, 400000, NA, 1302300000, c(53100000, 366400000, 428600000), c(187200000, 50400000),
                                                          NA, NA, NA, NA, 30400000, NA, NA, NA, NA, NA)),
                               subnation_pop_info_year = I(list(NA, "X2005", NA, "X2010", c("X2008","X2008","X2008"), c("X2005", "X2010"),
                                                           NA, NA, NA, NA, "X2000", NA, NA, NA, NA, NA)))

####################

####Det model

####################
#Run all model for every countries at the same time:
sel_country_for_model_run = ISO_country_list

for(coun in sel_country_for_model_run){
  if(coun == "PAK"){
    coun_model_run <- list(0.01)
  } else if(coun == "IND"){
    coun_model_run <- model_run_each_country(coun)$lambda %>% lapply(mean)
    coun_model_run <- c(list(0.01),coun_model_run)
  } else {
    coun_model_run <- model_run_each_country(coun)$lambda %>% lapply(mean)
  }
  parameters = list(lambdas = coun_model_run, sym_rate = (1/500 + 1/250)/2.0, death_rate = 0.2, dis_rate = 0.4)
  coun_model_gen <- model_gen_each_country(parameters, coun, 1)
  assign(paste(coun, "_model", sep = ""), coun_model_gen, envir = .GlobalEnv)
}

# l_coun_model_run = c()
# temp_dist = rlnorm(16000, log(0.01), 1)
# temp_dist_data = data.frame(rowname = "lambda_make", mean = mean(temp_dist), se_mean = sd(temp_dist)/sqrt(length(temp_dist)), sd = sd(temp_dist)
#                             , matrix(quantile(temp_dist, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)), nrow = 1, 
#                                      dimnames = list(1, c("2.5%", "25%", "50%", "75%", "97.5%"))))
# for(coun in sel_country_for_model_run){
#   if(coun == "PAK"){
#     coun_model_run <- cbind(country = coun, temp_dist_data)
#   } else if(coun == "IND"){
#     coun_model_run <- model_run_each_country(coun)$sum_data %>% lapply(function(x) data.frame(x[rownames(x) == "lambda" | rownames(x) == "rho",-c(9,10)])) %>% lapply(rownames_to_column)
#     coun_model_run <- c(lambda_make = list(temp_dist_data),coun_model_run)
#     coun_model_run <- data.frame(country = coun, Reduce(rbind, coun_model_run))
#   } else {
#     coun_model_run <- model_run_each_country(coun)$sum_data %>% lapply(function(x) data.frame(x[rownames(x) == "lambda" | rownames(x) == "rho",-c(9,10)])) %>% lapply(rownames_to_column)
#     coun_model_run <- data.frame(country = coun, Reduce(rbind, coun_model_run))
#   }
#   l_coun_model_run = c(l_coun_model_run, list(coun_model_run))
# }
# data_coun_run = Reduce(rbind, l_coun_model_run)
# write.csv(data_coun_run, "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/result_model/first_model/FOI_rho_sum.csv", row.names = F)
  
#Fill in template:
fill_in_Det_temp <- function(Det_template_input, l_country_code, Det_scen_data, scenario){
  for(country in 1:length(l_country_code)){
    Det_scen_data_coun = Det_scen_data[[country]]
    #Country index in the template:
    index_country = Det_template_input$country %in% l_country_code[country]
    
    Det_template_input[index_country, "deaths"] = Det_scen_data_coun[[scenario]]$deaths
    Det_template_input[index_country, "dalys"] = Det_scen_data_coun[[scenario]]$dalys
    Det_template_input[index_country, "cases"] = Det_scen_data_coun[[scenario]]$cases
    #cohort column - always the original population
    sel_naive_pop_demo = naive_pop_demo[(naive_pop_demo$age_from %in% age_group_sel) & (naive_pop_demo$country == l_country_code[country]),paste("X",year_sel, sep = "")]
    Det_template_input[index_country, "cohort_size"] = unlist(sel_naive_pop_demo)
  }
  return(Det_template_input)
}

sel_fill_model = ISO_country_list
All_sel_country_model = mget(paste(sel_fill_model, "model", sep = "_"))
###Deteministic output:
Det_No_vac = fill_in_Det_temp(Det_template_input = Det_template, l_country_code = sel_fill_model, Det_scen_data = All_sel_country_model, scenario = "naive")
Det_Cam_vac = fill_in_Det_temp(Det_template_input = Det_template, l_country_code = sel_fill_model, Det_scen_data = All_sel_country_model, scenario = "cam")
Det_Rou_vac = fill_in_Det_temp(Det_template_input = Det_template, l_country_code = sel_fill_model, Det_scen_data = All_sel_country_model, scenario = "rou")

#Add more countries or modify the output to the template:
# Det_template_save_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/result_model/first_model/Det_template/"
# Det_No_vac = read.csv(paste(Det_template_save_path, "test3_Det_No_vac.csv", sep = ""))
# Det_Cam_vac = read.csv(paste(Det_template_save_path, "test3_Det_Cam_vac.csv", sep = ""))
# Det_Rou_vac = read.csv(paste(Det_template_save_path, "test3_Det_Rou_vac.csv", sep = ""))
#China:
# Det_No_vac = fill_in_Det_temp(Det_template_input = Det_No_vac, l_country_code = "CHN", Det_scen_data = list(CHN_model), scenario = "naive")
# Det_Cam_vac = fill_in_Det_temp(Det_template_input = Det_Cam_vac, l_country_code = "CHN", Det_scen_data = list(CHN_model), scenario = "cam")
# Det_Rou_vac = fill_in_Det_temp(Det_template_input = Det_Rou_vac, l_country_code = "CHN", Det_scen_data = list(CHN_model), scenario = "rou")
#save filled in template:

Det_template_save_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/result_model/first_model/Det_template/"
write.csv(Det_No_vac, file = paste(Det_template_save_path, "test3_singlerun_Det_No_vac.csv", sep = ""), row.names = F)
write.csv(Det_Cam_vac, file = paste(Det_template_save_path, "test3_singlerun_Det_Cam_vac.csv", sep = ""), row.names = F)
write.csv(Det_Rou_vac, file = paste(Det_template_save_path, "test3_singlerun_Det_Rou_vac.csv", sep = ""), row.names = F)

################

######Generating stochastic model run_id

###############
run_times = 200
#sel set of countries
sel_country_for_model_run = ISO_country_list
#Sample random variables equal for every country, age group: asymptomatic rate, mortality rate, disability rate from # of run_id:
sample_sym_rate <- runif(run_times, 1/500, 1/250)
sample_death_rate <- runif(run_times, 0.1, 0.3)
sample_disa_rate <- runif(run_times, 0.3, 0.5)

###Save parameters file:
Sto_para_template = data.frame(run_id = 1:run_times) #Save para results
Sto_para_template[, c("<symptomatic_rate>", "<deaths_rate>", "<disability_rate>")] <- cbind(sample_sym_rate, sample_death_rate, sample_disa_rate)
for(coun in sel_country_for_model_run){
  
  if(coun == "PAK"){
    coun_model_run <- list(rlnorm(run_times, log(0.01), 1))
  } else if(coun == "IND"){
    coun_model_run <- model_run_each_country(coun)
    coun_model_run <- c(list(rlnorm(run_times, log(0.01), 1)),coun_model_run)
  } else {
    coun_model_run <- model_run_each_country(coun)
  }
  
  #save parameters:
  if(length(Sto_coun_model_run) == 1){
    Sto_para_template[,paste("<lambda>:<",coun, ">",sep ="")] <- Sto_coun_model_run[[1]]
  } else {
    for(i in 1:length(Sto_coun_model_run)){
      Sto_para_template[,paste("<lambda_",i,">:<",coun, ">",sep ="")] <- Sto_coun_model_run[[i]]
    }
  }
}
Sto_template_save_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/result_model/first_model/Sto_template/test3_16coun_200runid_every_countries/"
template_file_name = "stochastic_burden_template_JE-OUCRU-Clapham_"
write.csv(Sto_para_template, file = paste(Sto_template_save_path, template_file_name, "-parameters.csv", sep = ""), row.names = F)

###generate quantities from estimated FOI:
#read parameters from file:
Sto_template_save_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/result_model/first_model/Sto_template/test3_16coun_200runid_every_countries/"
template_file_name = "stochastic_burden_template_JE-OUCRU-Clapham_"
Sto_para_template <-read.csv(paste(Sto_template_save_path, template_file_name, "-parameters.csv", sep = ""), check.names = F)

#Save generated file:
for(i in 1:run_times){
  print(i)
  Sto_No_vac = c()
  Sto_Cam_vac = c()
  Sto_Rou_vac = c()
  for(coun in sel_country_for_model_run){
    Sto_coun_template <- filter(Sto_template, country == coun)
    Sto_coun_template$run_id <- rep(i, length(year_sel)*length(age_group_sel))
    Sto_coun_template$cohort_size = filter(naive_pop_demo, country == coun) %>% select(starts_with("X")) %>% unlist %>% unname
    
    sel_lambda = select(Sto_para_template, contains(paste("<",coun,">",sep = ""))) %>% lapply(c)
    parameters = list(lambdas = sel_lambda, sym_rate = Sto_para_template$`<symptomatic_rate>`, 
                      death_rate = Sto_para_template$`<deaths_rate>`, dis_rate = Sto_para_template$`<disability_rate>`)
    Sto_coun_model_gen <- model_gen_each_country(parameters, coun, i)
    
    sel_col_Sto_temp <- c("cases", "deaths", "dalys")
    
    #fill in template:
    Sto_No_vac_coun = Sto_coun_template
    Sto_No_vac_coun[,sel_col_Sto_temp] <- Sto_coun_model_gen$naive[sel_col_Sto_temp]
    Sto_Cam_vac_coun = Sto_coun_template
    Sto_Cam_vac_coun[,sel_col_Sto_temp] <- Sto_coun_model_gen$cam[sel_col_Sto_temp]
    Sto_Rou_vac_coun = Sto_coun_template
    Sto_Rou_vac_coun[,sel_col_Sto_temp] <- Sto_coun_model_gen$rou[sel_col_Sto_temp]
    
    #Combine all countries
    Sto_No_vac = c(Sto_No_vac, list(Sto_No_vac_coun))
    Sto_Cam_vac = c(Sto_Cam_vac, list(Sto_Cam_vac_coun))
    Sto_Rou_vac = c(Sto_Rou_vac, list(Sto_Rou_vac_coun))
  }
  
  Sto_No_vac = Reduce(rbind, Sto_No_vac)
  Sto_Cam_vac = Reduce(rbind, Sto_Cam_vac)
  Sto_Rou_vac = Reduce(rbind, Sto_Rou_vac)
  #save sto file:
  write.csv(Sto_No_vac, file = paste(Sto_template_save_path, template_file_name, "je-routine-no-vaccination_",i,".csv", sep = ""), row.names = F)
  write.csv(Sto_Cam_vac, file = paste(Sto_template_save_path, template_file_name, "je-campaign-gavi_",i,".csv", sep = ""), row.names = F)
  write.csv(Sto_Rou_vac, file = paste(Sto_template_save_path, template_file_name, "je-routine-gavi_",i,".csv", sep = ""), row.names = F)
}

###############

###Generate vaccinated population (for testing)

###############
#Function to generate no vac, cam, rou population overtime
model_gen_pop_scenario_each_country <- function(sel_country){
  print(sel_country)
  data_county_model = filter(All_country_model, country == sel_country)
  pop_gen <- filter(naive_pop_demo, country == sel_country)
  subnational_y_n = data_county_model$subnation
  sub_info_pop = unlist(data_county_model$subnation_pop_info_pop)
  sub_info_year = unlist(data_county_model$subnation_pop_info_year)
  
  #infer subnational data from vaccination population:
  l_pop_gen_naive = list()
  l_pop_gen_cam = list()
  l_pop_gen_rou = list()
  #if vaccination campaign is only on subnational scale: get the local pop first, then vaccinate
  if(subnational_y_n){
    subnational_gen = sub_pop_gen(sub_info_pop, sub_info_year, pop_gen)
    
    l_pop_gen_naive = unlist(subnational_gen[,paste("X", year_sel, sep = "")])
    l_pop_gen_cam = unlist(cam_scen_pop_gen(sel_country, subnational_gen, year_sel)[,paste("X", year_sel, sep = "")])
    l_pop_gen_rou = unlist(rou_scen_pop_gen(sel_country, subnational_gen, year_sel)[,paste("X", year_sel, sep = "")])
  } else { #if vaccination campaign is on national scale: get the pop vaccinate first, then get the local pop if there are many subregions
    
    pop_gen_naive = pop_gen
    pop_gen_cam = cam_scen_pop_gen(sel_country,pop_gen, year_sel)
    pop_gen_rou = rou_scen_pop_gen(sel_country,pop_gen, year_sel)
    
    if(length(sub_info_pop) != 1){
      for(sub in 1:length(sub_info_pop)){
        l_pop_gen_naive[[sub]] = unlist(pop_gen_naive[,paste("X", year_sel, sep = "")])*(sub_info_pop[sub]/sum(pop_gen[,sub_info_year[sub]]))
        l_pop_gen_cam[[sub]] = unlist(pop_gen_cam[,paste("X", year_sel, sep = "")])*(sub_info_pop[sub]/sum(pop_gen[,sub_info_year[sub]]))
        l_pop_gen_rou[[sub]] = unlist(pop_gen_rou[,paste("X", year_sel, sep = "")])*(sub_info_pop[sub]/sum(pop_gen[,sub_info_year[sub]]))
      }
      l_pop_gen_naive = Reduce("+", l_pop_gen_naive)
      l_pop_gen_cam = Reduce("+", l_pop_gen_cam)
      l_pop_gen_rou = Reduce("+", l_pop_gen_rou)
    } else {
      l_pop_gen_naive = unlist(pop_gen_naive[,paste("X", year_sel, sep = "")])
      l_pop_gen_cam = unlist(pop_gen_cam[,paste("X", year_sel, sep = "")])
      l_pop_gen_rou = unlist(pop_gen_rou[,paste("X", year_sel, sep = "")])
    }
  }
  l_pop_gen_scen <- list(naive = l_pop_gen_naive, cam = l_pop_gen_cam, rou = l_pop_gen_rou)
  return(l_pop_gen_scen)
}

Det_pop_each_scenario = data.frame(Det_template[,c("year", "age", "country", "country_name")], No_vac = NA, Cam = NA, Rou = NA)
for(coun in ISO_country_list){
  l_pop_scenario = model_gen_pop_scenario_each_country(coun)
  Det_pop_each_scenario[Det_pop_each_scenario$country == coun,c("No_vac", "Cam", "Rou")] = l_pop_scenario[c("naive", "cam", "rou")]
}
Det_pop_save_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/result_model/first_model/"
write.csv(Det_pop_each_scenario, file = paste(Det_pop_save_path, "Pop_scen.csv", sep = ""), row.names = F)
