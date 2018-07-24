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
library(dplyr)
library(data.table)

#data loading:
#Pop data load:#Montagu template load:
output_data_path = 'C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_clean/Montagu_data/test3/'
Det_template = read.csv(paste(output_data_path,'central_burden_template_JE-OUCRU-Clapham.csv', sep = ""))
Sto_template = read.csv(paste(output_data_path,'stochastic_burden_template_JE-OUCRU-Clapham.csv', sep = ""))
#Sto_para_template = read.csv(paste(output_data_path,"stochastic_template_params.csv", sep = ""))
#Sto_para_template = data.frame(run_id = Sto_para_template[,1])

#cases data load:
data_inci_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_raw/cases_sero_data/"
#save country result path:
coun_result_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/result_model/test3/"

#List of country need:
ISO_country_list = as.character(unique(Det_template$country))

#year select:
year_sel = 2000:2100
#age group select:
age_group_sel = 0:99
#Pop data load:
naive_pop_demo = read.xlsx2(paste(output_data_path,'naive_pop_1950_2100.xlsx', sep = ""), 1, colClasses=NA)
naive_pop_demo = naive_pop_demo[, c("country", "age_from", "age_to", paste("X", year_sel, sep = ""))]
temp_naive_pop_demo = read.xlsx2('C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_clean/pop_data/IDB_0_100_2000_2045_24.xlsx', 1, colClasses=NA)
#Expected life loss for each age in each country:
raw_expected_YLL_by_age = read.csv('C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/data_JE_raw/Montagu-data/201710gavi-2_dds-201710_life_ex_both.csv')
raw_expected_YLL_by_age = raw_expected_YLL_by_age[raw_expected_YLL_by_age$year %in% year_sel, c("country_code", "age_from", "age_to", "year","value")]

#run times of stochastic model:
run_times = 200

ll_case_je="

data {
  
  int N;  //number of age groups
  int l_age_group; //length of selected age group
  int age_seq[l_age_group]; // age group sequence
  vector[N] age_l; //lower bound of the age group
  vector[N] age_u; //upper bound of the age group
  vector[N] case_age; //number of cases in each age group.
  vector[N] pop_age; //population demo of the study year.
  int year; //number of col of pop_gen_case: year
  vector[year*l_age_group] pop_gen_case; //total population  with diff scenario
  vector[year*l_age_group] remained_life; // expected remained life for each age, from Montagu data.
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

generated quantities {
  vector[year*l_age_group] cases; //the generated cases from 2000 to 2044
  vector[year*l_age_group] mortality; //the mortality from 2000 to 2044
  vector[year*l_age_group] DALYs; //the DALYs from 2000 to 2044 = YLL + YLDacute + YLDchronic
  
  cases = uniform_rng(1.0/500, 1.0/250)*(1 - e()^(-lambda))*pop_gen_case .* to_vector(rep_matrix(exp((to_vector(age_seq)' - 1)*(-lambda)), year)');
  mortality = cases*uniform_rng(0.1,0.3);
  DALYs = mortality .* remained_life + cases*0.133*2.5/52.0 + (cases - mortality)*uniform_rng(0.3, 0.5)*0.542 .* remained_life;

}
"

stan_sero_model = stan_model(model_code = ll_case_je)

#loading generate data from diff
#Calculate for each country with different scenario:
stan_model_for_each_scenario <- function(l_data_inci_subnation, country_code, l_pop_gen_case, scenario){
  #remained life expectancy:
  raw_expected_remain_life_by_country = raw_expected_YLL_by_age[raw_expected_YLL_by_age$country_code == country_code,]
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
  
  #for each subnation data, we got a model output:
  model_output = list()
  for(sub_country in 1:length(l_data_inci_subnation)){
    #Calculate the age_l and age_u:
    age_group_split = strsplit(as.character(l_data_inci_subnation[[sub_country]]$Age_group),"-")
    age_l = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[1])}))
    age_u = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[2])}))
    #MCMC:
    data_for_HMC = list(N = nrow(l_data_inci_subnation[[sub_country]]), year = length(year_sel),
                        age_l = age_l, age_u = age_u, l_age_group = length(age_group_sel), age_seq = 1:length(age_group_sel),
                        case_age = l_data_inci_subnation[[sub_country]]$Case_Sero, 
                        pop_age = l_data_inci_subnation[[sub_country]]$Pop_all_age, 
                        pop_gen_case = l_pop_gen_case[[sub_country]],
                        t_case = sum(l_data_inci_subnation[[sub_country]]$Case_Sero),
                        remained_life = infer_remain_life_2)
    
    stan_FOI_fit = sampling(object = stan_sero_model, data = data_for_HMC, 
                            chains = 4, iter = 8000, thin = 4)
    #print(stan_FOI_fit)
    sum_data <- data.frame(summary(stan_FOI_fit)$summary)
    
    
    #1st and 2nd column of output data:
    year = rep(year_sel, length(age_group_sel))
    age =  rep(age_group_sel, each = length(year_sel))
    #Cases estimated:
    case_proj_index = grep("cases" ,rownames(sum_data))
    case_proj = sum_data[case_proj_index,]
    case_proj = cbind(year, age, case_proj)
    case_proj$Scenario = scenario
    
    #Mortality estimated:
    mo_proj_index = grep("mortality" ,rownames(sum_data))
    mo_proj = sum_data[mo_proj_index,]
    mo_proj = cbind(year, age, mo_proj)
    mo_proj$Scenario = scenario
    
    #DALY estimated:
    DALY_proj_index = grep("DALYs" ,rownames(sum_data))
    DALY_proj = sum_data[DALY_proj_index,]
    DALY_proj = cbind(year, age, DALY_proj)
    DALY_proj$Scenario = scenario
    
    #samples lambda:
    lambda_extracted = extract(stan_FOI_fit, pars = "lambda")$lambda
    lambda_samples = sample(lambda_extracted, run_times)
    
    if(length(l_data_inci_subnation) != 1){
      cases_ext = extract(stan_FOI_fit, pars = "cases")$cases
      death_ext = extract(stan_FOI_fit, pars = "mortality")$mortality
      dalys_ext = extract(stan_FOI_fit, pars = "DALYs")$DALYs
      model_output[[sub_country]] = list(cases = case_proj, deaths = mo_proj, dalys = DALY_proj, 
           summary = sum_data, samples = list(lambda = lambda_samples), 
           extract = list(cases = cases_ext, deaths = death_ext, dalys = dalys_ext))
    } else {
      model_output = list(cases = case_proj, deaths = mo_proj, dalys = DALY_proj, 
                          summary = sum_data, samples = list(lambda = lambda_samples))
    }
  }
  
  
  return(model_output)
}

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

#function to combine data:
#list_data: list of subnational data
#scenario: "Naive", "Campaign". or "Routine"
combine_data <- function(list_data, scenario){
  #cases:
  list_cases = lapply(list_data, FUN = f <- function(x)x$extract$cases)
  sum_cases = Reduce("+", list_cases)
  
  #deaths:
  list_deaths = lapply(list_data, FUN = f <- function(x){x$extract$deaths})
  sum_deaths = Reduce("+", list_deaths)
  
  #dalys:
  list_dalys = lapply(list_data, FUN = f <- function(x){x$extract$dalys})
  sum_dalys = Reduce("+", list_dalys)
  
  summary_data_function = function(data) {data.frame(mean = colMeans(data), 
                                                     X2.5. = apply(data, 2, f <- function(x)quantile(x, 0.025)),
                                                     X97.5. = apply(data, 2, f <- function(x)quantile(x, 0.975)))}
  
  #lambda:
  list_lambda = lapply(list_data, FUN = f <- function(x)x$samples)
  
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
  
  return(list(cases = case_proj, deaths = mo_proj, dalys = DALY_proj, samples = list(lambda = list_lambda)))
}

#get the function of generating campaign and routine population:
source("C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/code_JEV_model/JEV_Montagu_data - Scen_pop_generate_function.R")

#function to generate model results with 3 different scenario:
#pop_gen: data frame of population over selected time by each selected group age 
#country_code: 3 letter code the country
#subnational: does the campaign vaccination is on subnational scale or national scale
#l_country_of_incidence_data: name of the file to generate FOI
#l_sel_subnation_index: the index of selected incidence data from l_country_of_incidence_data
#out_put_name: name of the output file
model_result_gen <- function(pop_gen, subnational, l_country_of_incidence_data, l_sel_subnation_index, l_subnation_info){
  country_code = as.character(pop_gen$country[1])
  
  #incidence data to generate FOI:
  l_data_inci_subnation = c()
  for(n_inci_data in 1:length(l_country_of_incidence_data)){
    data_inci = read.xlsx(paste(data_inci_path, l_country_of_incidence_data[n_inci_data], "_JE.xlsx", sep = ""), 1)
    subnation_inci = unlist(unique(data_inci$subnation))
    l_data_inci_subnation = c(l_data_inci_subnation,list(data_inci[data_inci$subnation == subnation_inci[l_sel_subnation_index[n_inci_data]],]))
  }
  
  #infer subnational data from vaccination population:
  l_pop_gen_naive = list()
  l_pop_gen_cam = list()
  l_pop_gen_rou = list()
  #if vaccination campaign is only on subnational scale: get the local pop first, then vaccinate
  if(subnational){
    subnational_gen = sub_pop_gen(l_subnation_info$pop, l_subnation_info$year, pop_gen)

    l_pop_gen_naive = list(unlist(subnational_gen[,paste("X", year_sel, sep = "")]))
    l_pop_gen_cam = list(unlist(cam_scen_pop_gen(country_code, subnational_gen, year_sel)[,paste("X", year_sel, sep = "")]))
    l_pop_gen_rou = list(unlist(rou_scen_pop_gen(country_code, subnational_gen, year_sel)[,paste("X", year_sel, sep = "")]))
  } else { #if vaccination campaign is on national scale: get the pop vaccinate first, then get the local pop if there are many subregions

    pop_gen_naive = pop_gen
    pop_gen_cam = cam_scen_pop_gen(country_code,pop_gen, year_sel )
    pop_gen_rou = rou_scen_pop_gen(country_code,pop_gen, year_sel)
    
    if(length(l_data_inci_subnation) != 1){
      for(sub in 1:length(l_subnation_info)){
        l_pop_gen_naive[[sub]] = unlist(pop_gen_naive[,paste("X", year_sel, sep = "")])*(l_subnation_info$pop[sub]/sum(pop_gen[,l_subnation_info$year[sub]]))
        l_pop_gen_cam[[sub]] = unlist(pop_gen_cam[,paste("X", year_sel, sep = "")])*(l_subnation_info$pop[sub]/sum(pop_gen[,l_subnation_info$year[sub]]))
        l_pop_gen_rou[[sub]] = unlist(pop_gen_rou[,paste("X", year_sel, sep = "")])*(l_subnation_info$pop[sub]/sum(pop_gen[,l_subnation_info$year[sub]]))
      }
    } else {
      l_pop_gen_naive = list(unlist(pop_gen_naive[,paste("X", year_sel, sep = "")]))
      l_pop_gen_cam = list(unlist(pop_gen_cam[,paste("X", year_sel, sep = "")]))
      l_pop_gen_rou = list(unlist(pop_gen_rou[,paste("X", year_sel, sep = "")]))
    }
  }
  
  #model run:
  naive_scen = stan_model_for_each_scenario(l_data_inci_subnation, country_code, l_pop_gen_case = l_pop_gen_naive, scenario = "Naive")
  cam_scen = stan_model_for_each_scenario(l_data_inci_subnation, country_code, l_pop_gen_cam, "Campaign")
  rou_scen = stan_model_for_each_scenario(l_data_inci_subnation, country_code, l_pop_gen_rou, "Routine")
  
  #combine data if applicable:
  if(length(l_data_inci_subnation) != 1){
    comb_naive = combine_data(naive_scen, "Naive")
    comb_cam = combine_data(cam_scen, "Campaign")
    comb_rou = combine_data(rou_scen, "Routine")
    #generate output:
    return(list(naive = comb_naive, cam = comb_cam, rou = comb_rou))
  } else {
    #generate output:
    return(list(naive = naive_scen, cam = cam_scen, rou = rou_scen))
  }
}

#Function to generate cases, mortality, DALYS from a value of lambda:
det_model_from_lambda <-function(lambda, pop_gen, subnational, sub_info){
    country_code = as.character(pop_gen$country[1])
    #remained life expectancy:
    raw_expected_remain_life_by_country = raw_expected_YLL_by_age[raw_expected_YLL_by_age$country_code == country_code,]
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

    if(subnational){
      #get local pop first:
      pop_gen = sub_pop_gen(sub_info$pop, sub_info$year, pop_gen)
      #vaccinate later:
      pop_gen_naive = unlist(pop_gen[,paste("X", year_sel, sep = "")])
      pop_gen_cam = unlist(cam_scen_pop_gen(country_code, pop_gen, year_sel )[,paste("X", year_sel, sep = "")])
      pop_gen_rou = unlist(rou_scen_pop_gen(country_code, pop_gen, year_sel)[,paste("X", year_sel, sep = "")])
    } else {
      #get the vaccination first:
      pop_gen_naive = pop_gen
      pop_gen_cam = cam_scen_pop_gen(country_code, pop_gen, year_sel )
      pop_gen_rou = rou_scen_pop_gen(country_code, pop_gen, year_sel)
      #get the proportional local pop after vaccinated:
      pop_gen_naive = unlist(pop_gen_naive[,paste("X", year_sel, sep = "")])*sub_info$pop/sum(pop_gen[,sub_info$year])
      pop_gen_cam = unlist(pop_gen_cam[,paste("X", year_sel, sep = "")])*sub_info$pop/sum(pop_gen[,sub_info$year])
      pop_gen_rou = unlist(pop_gen_rou[,paste("X", year_sel, sep = "")])*sub_info$pop/sum(pop_gen[,sub_info$year])
    }
    
    list_pop_gen = list(pop_gen_naive, pop_gen_cam, pop_gen_rou)
    scenarios = c("naive", "cam", "rou")
    scen_list = list()
    for(i in 1:length(scenarios)){
      pop_gen_case_less_15 = list_pop_gen[[i]]
      l_data_output = length(pop_gen_case_less_15)
      cases_gen = rep((1/500 + 1/250)/2.0*(1 - exp(-lambda)), l_data_output)*pop_gen_case_less_15* #central estimation of symptomatic rate
        exp(rep(age_group_sel, each = length(year_sel))*(-rep(lambda, l_data_output)))
      mort_gen = cases_gen*rep(0.2, l_data_output) #central estimation of death rate
      disa_gen = (cases_gen - mort_gen)*rep(0.4, l_data_output)
      DALYS_gen = mort_gen*infer_remain_life_2 + cases_gen*rep(0.133*2.5/52.0, l_data_output) + 
        disa_gen*rep(0.542, l_data_output)*infer_remain_life_2
      scen_list[[i]] = list(cases = list(mean = cases_gen), 
                            deaths = list(mean = mort_gen), 
                            dalys = list(mean = DALYS_gen))
    }
    names(scen_list) = scenarios
    return(scen_list)
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

#Run all model for every countries at the same time:
model_Rdata_save_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/result_model/first_model/Model_Rdata_storage/"
sel_country_for_model_run = setdiff(ISO_country_list,c("IND", "PAK"))

for(country_ID in sel_country_for_model_run){
  country_model_data = filter(All_country_model, country == country_ID)
  naive_pop_country = filter(naive_pop_demo, country == country_ID) #population of the country overtime
  country_sunation_info = list(pop = unlist(country_model_data$subnation_pop_info_pop), year = unlist(country_model_data$subnation_pop_info_year))
  country_model = model_result_gen(pop_gen = naive_pop_country, subnational = country_model_data$subnation, 
                                   l_country_of_incidence_data = unlist(country_model_data$countries_inci_data), 
                                   l_sel_subnation_index = unlist(country_model_data$index_sel_inci_data), 
                                   l_subnation_info = country_sunation_info)
  saveRDS(country_model, file = paste(model_Rdata_save_path,country_ID,"_model.RData", sep = ""))
  #a = readRDS(paste(model_Rdata_save_path,country_ID,"_model.RData", sep = ""))
}

##############each country run:
{##India:

##Population data of India overtime:
naive_pop_India = naive_pop_demo[naive_pop_demo$country == "IND",]

##Lowest incidence population: "Haryana","Punjab"
IND_lowest_info = list(pop = 53100000, year = "X2008")
#Run model:
IND_lowest_model = det_model_from_lambda(lambda = 0.01, pop_gen = naive_pop_India, subnational = F, sub_info = IND_lowest_info)

India_sub_info = list(pop = c(366400000, 428600000), year = c("X2008", "X2008"))
#Run model:
IND_model = model_result_gen(pop_gen = naive_pop_India, subnational = F, 
                 l_country_of_incidence_data = c("Thailand", "Nepal"), l_sel_subnation_index = c(1, 8), l_subnation_info = India_sub_info)

######################
##Parkistan:

##Population data of Parkistan overtime:
naive_pop_Parkistan = naive_pop_demo[naive_pop_demo$country == "PAK",]
PAK_sub_info = list(pop = 30400000, year = "X2000")
#Run model:
PAK_model = det_model_from_lambda(lambda = 0.01, pop_gen = naive_pop_Parkistan, subnational = T, sub_info = PAK_sub_info)

######################
##Cambodia:
##Population data of Cambodia overtime:
naive_pop_Cambodia = naive_pop_demo[naive_pop_demo$country == "KHM",]
#Run model:
KHM_model = model_result_gen(pop_gen = naive_pop_Cambodia, subnational = F, 
                 l_country_of_incidence_data = "Cambodia", l_sel_subnation_index = 3)

######################
##Indonesia: => stratify to high and med
##Population data of Indonesia overtime:
naive_pop_Indonesia = naive_pop_demo[naive_pop_demo$country == "IDN",]
##low incidence population:all provinces of sumatra, Java and Bapua, plus Kepulauan, Bangka Belitung and Riau
##high incidence population: Bali, Nusa Tenggara, all prov in Borneo and Sulawesi, and the Moluccas
#Run model:
IDN_sub_info = list(pop = c(187200000, 50400000), year = c("X2005", "X2010"))
IDN_model = model_result_gen(pop_gen = naive_pop_Indonesia, subnational = F, 
                             l_country_of_incidence_data = c("Thailand", "Indonesia"), l_sel_subnation_index = c(1, 1), l_subnation_info = IDN_sub_info)
######################
##Laos:
##Population data of Laos overtime:
naive_pop_Laos = naive_pop_demo[naive_pop_demo$country == "LAO",]
#Run model:
LAO_model = model_result_gen(pop_gen = naive_pop_Laos, subnational = F, 
                              l_country_of_incidence_data = "Laos", l_sel_subnation_index = 1)

######################
##Vietnam:
##Population data of Vietnam overtime:
naive_pop_Vietnam = naive_pop_demo[naive_pop_demo$country == "VNM",]
#Run model:
VNM_model = model_result_gen(pop_gen = naive_pop_Vietnam, subnational = F, 
                              l_country_of_incidence_data = "Laos", l_sel_subnation_index = 1)

######################
##Bangladesh:
##Population data of Bangladesh overtime:
naive_pop_Bangladesh = naive_pop_demo[naive_pop_demo$country == "BGD",]
#Run model:
BGD_model = model_result_gen(pop_gen = naive_pop_Bangladesh, subnational = F, 
                              l_country_of_incidence_data = "Bangladesh", l_sel_subnation_index = 1)

######################
##Nepal: =>stratify to high (24 Terai and inner Terai districts) and low (51 mountain and hill districts)
##Population data of Nepal overtime:
naive_pop_Nepal = naive_pop_demo[naive_pop_demo$country == "NPL",]
#Run model:
NPL_model = model_result_gen(pop_gen = naive_pop_Nepal, subnational = F, 
                                  l_country_of_incidence_data = "Nepal", l_sel_subnation_index = 9)

######################
##Butan: => data from Bangladesh, but effect only local: southern foothills
##Population data of Butan overtime:
naive_pop_Butan = naive_pop_demo[naive_pop_demo$country == "BTN",]
BTN_sub_info = list(pop = 400000, year = "X2005")
#Run model:
BTN_model = model_result_gen(pop_gen = naive_pop_Butan, subnational = T, 
                             l_country_of_incidence_data = "Nepal", l_sel_subnation_index = 9, l_subnation_info = BTN_sub_info)

######################
##North Korea:
naive_pop_N_Kor = naive_pop_demo[naive_pop_demo$country == "PRK",]
#Run model:
PRK_model = model_result_gen(pop_gen = naive_pop_N_Kor, subnational = F, 
                                 l_country_of_incidence_data = "Thailand", l_sel_subnation_index = 1)

######################
##Burma: => data from Thailand
##Population data of Burma overtime:
naive_pop_Burma = naive_pop_demo[naive_pop_demo$country == "MMR",]
#Run model:
MMR_model = model_result_gen(pop_gen = naive_pop_Burma, subnational = F, 
                             l_country_of_incidence_data = "Thailand", l_sel_subnation_index = 1)

######################
##PNG:
##Population data of Burma overtime:
naive_pop_PNG = naive_pop_demo[naive_pop_demo$country == "PNG",]
#Run model:
PNG_model = model_result_gen(pop_gen = naive_pop_PNG, subnational = F, 
                               l_country_of_incidence_data = "Thailand", l_sel_subnation_index = 1)

######################
##Srilanka: => data from Thailand
##Population data of Srilanka overtime:
naive_pop_Srilanka = naive_pop_demo[naive_pop_demo$country == "LKA",]
#Run model:
LKA_model = model_result_gen(pop_gen = naive_pop_Srilanka, subnational = F, 
                               l_country_of_incidence_data = "Thailand", l_sel_subnation_index = 1)

######################
##Timor-leste: 
##Population data of Timor overtime:
naive_pop_Timor = naive_pop_demo[naive_pop_demo$country == "TLS",]
#Run model:
TLS_model = model_result_gen(pop_gen = naive_pop_Timor, subnational = F, 
                               l_country_of_incidence_data = "Indonesia", l_sel_subnation_index = 1)

######################
##China:
##Population data of China overtime:
naive_pop_China = naive_pop_demo[naive_pop_demo$country == "CHN",]
CHN_sub_info = list(pop = (1025.7 + 276.6)*1000000, year = "X2010")
#Run model:
CHN_model = model_result_gen(pop_gen = naive_pop_China, subnational = T, 
                             l_country_of_incidence_data = "China", l_sel_subnation_index = 23, l_subnation_info = CHN_sub_info)

######################
##Philippines:
##Population data of Philippines overtime:
naive_pop_Philippines = naive_pop_demo[naive_pop_demo$country == "PHL",]
#Run model:
PHL_model = model_result_gen(pop_gen = naive_pop_Philippines, subnational = F, 
                             l_country_of_incidence_data = "Philippines", l_sel_subnation_index = 1)

}

##################################################

################Fill in template:

##################################################
#India and pak are different:
fill_in_Det_temp <- function(Det_template_input, l_country_code, Det_scen_data, scenario){
  for(country in 1:length(l_country_code)){
    Det_scen_data_coun = Det_scen_data[[country]]
    #Country index in the template:
    index_country = Det_template_input$country %in% l_country_code[country]
    
    Det_template_input[index_country, "deaths"] = Det_scen_data_coun[[scenario]]$deaths$mean
    Det_template_input[index_country, "dalys"] = Det_scen_data_coun[[scenario]]$dalys$mean
    Det_template_input[index_country, "cases"] = Det_scen_data_coun[[scenario]]$cases$mean
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

#Add more countries to the template:
Det_template_save_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/result_model/first_model/Det_template/"
Det_No_vac = read.csv(paste(Det_template_save_path, "test3_Det_No_vac.csv", sep = ""))
Det_Cam_vac = read.csv(paste(Det_template_save_path, "test3_Det_Cam_vac.csv", sep = ""))
Det_Rou_vac = read.csv(paste(Det_template_save_path, "test3_Det_Rou_vac.csv", sep = ""))

#India:
comb_IND = list(naive = list(cases = list(mean = IND_model$naive$cases$mean + IND_lowest_model$naive$cases$mean),
                             deaths = list(mean = IND_model$naive$deaths$mean + IND_lowest_model$naive$deaths$mean),
                             dalys = list(mean = IND_model$naive$dalys$mean + IND_lowest_model$naive$dalys$mean)),
                cam = list(cases = list(mean = IND_model$cam$cases$mean + IND_lowest_model$cam$cases$mean),
                           deaths = list(mean = IND_model$cam$deaths$mean + IND_lowest_model$cam$deaths$mean),
                           dalys = list(mean = IND_model$cam$dalys$mean + IND_lowest_model$cam$dalys$mean)),
                rou = list(cases = list(mean = IND_model$rou$cases$mean + IND_lowest_model$rou$cases$mean),
                           deaths = list(mean = IND_model$rou$deaths$mean + IND_lowest_model$rou$deaths$mean),
                           dalys = list(mean = IND_model$rou$dalys$mean + IND_lowest_model$rou$dalys$mean)))
Det_No_vac = fill_in_Det_temp(Det_template_input = Det_No_vac, l_country_code = "IND", Det_scen_data = list(comb_IND), scenario = "naive")
Det_Cam_vac = fill_in_Det_temp(Det_template_input = Det_Cam_vac, l_country_code = "IND", Det_scen_data = list(comb_IND), scenario = "cam")
Det_Rou_vac = fill_in_Det_temp(Det_template_input = Det_Rou_vac, l_country_code = "IND", Det_scen_data = list(comb_IND), scenario = "rou")

#China:
Det_No_vac = fill_in_Det_temp(Det_template_input = Det_No_vac, l_country_code = "CHN", Det_scen_data = list(CHN_model), scenario = "naive")
Det_Cam_vac = fill_in_Det_temp(Det_template_input = Det_Cam_vac, l_country_code = "CHN", Det_scen_data = list(CHN_model), scenario = "cam")
Det_Rou_vac = fill_in_Det_temp(Det_template_input = Det_Rou_vac, l_country_code = "CHN", Det_scen_data = list(CHN_model), scenario = "rou")
#Philippines:
Det_No_vac = fill_in_Det_temp(Det_template_input = Det_No_vac, l_country_code = "PHL", Det_scen_data = list(PHL_model), scenario = "naive")
Det_Cam_vac = fill_in_Det_temp(Det_template_input = Det_Cam_vac, l_country_code = "PHL", Det_scen_data = list(PHL_model), scenario = "cam")
Det_Rou_vac = fill_in_Det_temp(Det_template_input = Det_Rou_vac, l_country_code = "PHL", Det_scen_data = list(PHL_model), scenario = "rou")
#Indonesia:
Det_No_vac = fill_in_Det_temp(Det_template_input = Det_No_vac, l_country_code = "IDN", Det_scen_data = list(IDN_model), scenario = "naive")
Det_Cam_vac = fill_in_Det_temp(Det_template_input = Det_Cam_vac, l_country_code = "IDN", Det_scen_data = list(IDN_model), scenario = "cam")
Det_Rou_vac = fill_in_Det_temp(Det_template_input = Det_Rou_vac, l_country_code = "IDN", Det_scen_data = list(IDN_model), scenario = "rou")

#save filled in template:

Det_template_save_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/result_model/first_model/Det_template/"
write.csv(Det_No_vac, file = paste(Det_template_save_path, "test3_Det_No_vac.csv", sep = ""), row.names = F)
write.csv(Det_Cam_vac, file = paste(Det_template_save_path, "test3_Det_Cam_vac.csv", sep = ""), row.names = F)
write.csv(Det_Rou_vac, file = paste(Det_template_save_path, "test3_Det_Rou_vac.csv", sep = ""), row.names = F)

##################################################

#################Generating stochastic model run_id

##################################################
#Reduced model which only generate lambda samples:
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

#sel set of countries
# sel_country_for_model_run = setdiff(ISO_country_list,c("IND", "PAK"))
# sel_country_for_model_run = c("IND", "PAK", "IDN")
sel_country_for_model_run = ISO_country_list
#sel_country_for_model_run = setdiff(ISO_country_list,c("IND", "PAK", "BGD", "BTN", "KHM", "CHN"))
#Sample random variables equal for every country, age group: asymptomatic rate, mortality rate, disability rate from # of run_id:
sample_sym_rate <- runif(run_times, 1/500, 1/250)
sample_death_rate <- runif(run_times, 0.1, 0.3)
sample_disa_rate <- runif(run_times, 0.3, 0.5)
Sto_para_template = data.frame(run_id = 1:run_times) #Save para results
Sto_para_template[, c("<symptomatic_rate>", "<deaths_rate>", "<disability_rate>")] <- cbind(sample_sym_rate, sample_death_rate, sample_disa_rate)
#Save results:
Sto_template_save_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/result_model/first_model/Sto_template/test3_16coun_200runid_every_countries/"
template_file_name = "stochastic_burden_template_JE-OUCRU-Clapham_"

#read from file:
#Sto_para_template <-read.csv(paste(Sto_template_save_path, template_file_name, "-parameters.csv", sep = ""), check.names = F)
# sample_sym_rate <- Sto_para_template$`<symptomatic_rate>`
# sample_death_rate <- Sto_para_template$`<deaths_rate>`
# sample_disa_rate <- Sto_para_template$`<disability_rate>`

#function to run model for each country => output is sampled of lambda:
Sto_model_run_each_country <- function(sel_country){
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
    lambda_sampled = sample(lambda_all, run_times)
    l_lambda_sampled = c(l_lambda_sampled, list(lambda_sampled))
  }
  return(l_lambda_sampled = l_lambda_sampled)
}

#Save parameters file:
for(coun in sel_country_for_model_run){
  
  if(coun == "PAK"){
    Sto_coun_model_run <- list(rlnorm(run_times, log(0.01), 1))
  } else if(coun == "IND"){
    Sto_coun_model_run <- Sto_model_run_each_country(coun)
    Sto_coun_model_run <- c(list(rlnorm(run_times, log(0.01), 1)),Sto_coun_model_run)
  } else {
    Sto_coun_model_run <- Sto_model_run_each_country(coun)
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
write.csv(Sto_para_template, file = paste(Sto_template_save_path, template_file_name, "-parameters.csv", sep = ""), row.names = F)

###generate quantities from estimated FOI:
Sto_para_template <-read.csv(paste(Sto_template_save_path, template_file_name, "-parameters.csv", sep = ""), check.names = F)
sample_sym_rate <- Sto_para_template$`<symptomatic_rate>`
sample_death_rate <- Sto_para_template$`<deaths_rate>`
sample_disa_rate <- Sto_para_template$`<disability_rate>`
#function to generate cases, deaths, dalys, Sto model:
Sto_model_gen_each_country <- function(l_lambda_sampled, sel_country, run){
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
    l_pop_gen_cam = list(unlist(cam_scen_pop_gen(sel_country, subnational_gen, year_sel )[,paste("X", year_sel, sep = "")]))
    l_pop_gen_rou = list(unlist(rou_scen_pop_gen(sel_country, subnational_gen, year_sel)[,paste("X", year_sel, sep = "")]))
  } else { #if vaccination campaign is on national scale: get the pop vaccinate first, then get the local pop if there are many subregions
    
    pop_gen_naive = pop_gen
    pop_gen_cam = cam_scen_pop_gen(sel_country,pop_gen, year_sel )
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
  whole_country_output <- c()
  for(i in 1:length(l_pop_gen_scen)){
    pop_scen <- l_pop_gen_scen[[i]]
    l_data_output <- length(pop_scen[[1]])
    cases_scen_list <- list() #combine if there are multiple data
    
    cases_all_runs <- c()
    deaths_all_runs <- c()
    dalys_all_runs <- c()
    
      for(n_inci_data in 1:length(l_lambda_sampled)){
        lambda_sampled = l_lambda_sampled[[n_inci_data]][run]
        cases_scen_list[[n_inci_data]] = rep(sample_sym_rate[run]*(1 - exp(-lambda_sampled)), l_data_output)*pop_scen[[n_inci_data]]* 
          exp(rep(age_group_sel,length(year_sel))*(-rep(lambda_sampled, l_data_output)))
      }
      cases_gen <- Reduce("+", cases_scen_list) %>% unname
      mort_gen = cases_gen*rep(sample_death_rate[run], l_data_output) 
      disa_gen = (cases_gen - mort_gen)*rep(sample_disa_rate[run], l_data_output)
      DALYS_gen = mort_gen*infer_remain_life_2 + cases_gen*rep(0.133*2.5/52.0, l_data_output) + 
        disa_gen*rep(0.542, l_data_output)*infer_remain_life_2
      
    whole_country_output[[i]] <- list(cases = cases_gen, deaths = mort_gen, dalys = DALYS_gen) #gen data for each scenario:
  }
  
  names(whole_country_output) <- c("naive", "cam", "rou")
  
  return(whole_country_output)
}

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
    Sto_coun_model_gen <- Sto_model_gen_each_country(sel_lambda, coun, i)
    
    sel_col_Sto_temp <- c("cases", "deaths", "dalys")
    
    #fill in template:
    Sto_No_vac_coun = Sto_coun_template
    Sto_No_vac_coun[,sel_col_Sto_temp] <- lapply(Sto_coun_model_gen$naive[sel_col_Sto_temp], round)
    Sto_Cam_vac_coun = Sto_coun_template
    Sto_Cam_vac_coun[,sel_col_Sto_temp] <- lapply(Sto_coun_model_gen$cam[sel_col_Sto_temp], round)
    Sto_Rou_vac_coun = Sto_coun_template
    Sto_Rou_vac_coun[,sel_col_Sto_temp] <- lapply(Sto_coun_model_gen$rou[sel_col_Sto_temp], round)
    
    #Combine all countries
    Sto_No_vac = c(Sto_No_vac, list(Sto_No_vac_coun))
    Sto_Cam_vac = c(Sto_Cam_vac, list(Sto_Cam_vac_coun))
    Sto_Rou_vac = c(Sto_Rou_vac, list(Sto_Rou_vac_coun))
  }
  
  Sto_No_vac = Reduce(rbind, Sto_No_vac)
  Sto_Cam_vac = Reduce(rbind, Sto_Cam_vac)
  Sto_Rou_vac = Reduce(rbind, Sto_Rou_vac)
  #save sto file:
  fwrite(Sto_No_vac, file = paste(Sto_template_save_path, template_file_name, "je-routine-no-vaccination_",i,".csv", sep = ""))
  fwrite(Sto_Cam_vac, file = paste(Sto_template_save_path, template_file_name, "je-campaign-gavi_",i,".csv", sep = ""))
  fwrite(Sto_Rou_vac, file = paste(Sto_template_save_path, template_file_name, "je-routine-gavi_",i,".csv", sep = ""))
}

################Different Sto model: sample 200 run ids of central estimates distribution, keep the sample symp rate, dis rate, death rate

ll_case_je="

data {

  int N;  //number of age groups
  int l_age_group; //length of selected age group
  int age_seq[l_age_group]; // age group sequence
  vector[N] age_l; //lower bound of the age group
  vector[N] age_u; //upper bound of the age group
  vector[N] case_age; //number of cases in each age group.
  vector[N] pop_age; //population demo of the study year.
  int year; //number of col of pop_gen_case: year
  vector[year*l_age_group] pop_gen_case; //total population  with diff scenario
  vector[year*l_age_group] remained_life; // expected remained life for each age, from Montagu data.
  real t_case; //total cases in every groups

}

parameters {

  real<lower = 0> lambda;
  real<lower = 0> rho;
  real<lower = 0> lambda;
  real<lower = 0> rho;
  real<lower = 0> lambda;

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

  vector[year*l_age_group] cases; //the generated cases from 2000 to 2044
  vector[year*l_age_group] mortality; //the mortality from 2000 to 2044
  vector[year*l_age_group] DALYs; //the DALYs from 2000 to 2044 = YLL + YLDacute + YLDchronic
  
  cases = uniform_rng(1.0/500, 1.0/250)*(1 - e()^(-lambda))*pop_gen_case .* to_vector(rep_matrix(exp((to_vector(age_seq)' - 1)*(-lambda)), year)');
  mortality = cases*uniform_rng(0.1,0.3);
  DALYs = mortality .* remained_life + cases*0.133*2.5/52.0 + (cases - mortality)*uniform_rng(0.3, 0.5)*0.542 .* remained_life;

}
"

stan_sero_model = stan_model(model_code = ll_case_je)

#sample run_times (200) central estimates from the distribution:
stan_model_for_each_scenario <- function(l_data_inci_subnation, country_code, l_pop_gen_case, scenario){
  #remained life expectancy:
  raw_expected_remain_life_by_country = raw_expected_YLL_by_age[raw_expected_YLL_by_age$country_code == country_code,]
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
  
  #for each subnation data, we got a model output:
  model_output = list()
  for(sub_country in 1:length(l_data_inci_subnation)){
    #Calculate the age_l and age_u:
    age_group_split = strsplit(as.character(l_data_inci_subnation[[sub_country]]$Age_group),"-")
    age_l = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[1])}))
    age_u = unlist(lapply(age_group_split, FUN = f <- function(x){as.numeric(x[2])}))
    #MCMC:
    data_for_HMC = list(N = nrow(l_data_inci_subnation[[sub_country]]), year = length(year_sel),
                        age_l = age_l, age_u = age_u, l_age_group = length(age_group_sel), age_seq = 1:length(age_group_sel),
                        case_age = l_data_inci_subnation[[sub_country]]$Case_Sero, 
                        pop_age = l_data_inci_subnation[[sub_country]]$Pop_all_age, 
                        pop_gen_case = l_pop_gen_case[[sub_country]],
                        t_case = sum(l_data_inci_subnation[[sub_country]]$Case_Sero),
                        remained_life = infer_remain_life_2)
    
    stan_FOI_fit = sampling(object = stan_sero_model, data = data_for_HMC, 
                            chains = 4, iter = 8000, thin = 80)
    
    output_extracted = extract(stan_FOI_fit, pars = c("lambda", "cases", "mortality", "DALYs"))
    
    if(length(l_data_inci_subnation) != 1){
      cases_ext = extract(stan_FOI_fit, pars = "cases")$cases
      death_ext = extract(stan_FOI_fit, pars = "mortality")$mortality
      dalys_ext = extract(stan_FOI_fit, pars = "DALYs")$DALYs
      model_output[[sub_country]] = list(cases = case_proj, deaths = mo_proj, dalys = DALY_proj, 
                                         summary = sum_data, samples = list(lambda = lambda_samples), 
                                         extract = list(cases = cases_ext, deaths = death_ext, dalys = dalys_ext))
    } else {
      model_output = list(cases = case_proj, deaths = mo_proj, dalys = DALY_proj, 
                          summary = sum_data, samples = list(lambda = lambda_samples))
    }
  }
  
  
  return(model_output)
}

Sto_template_save_path = "C:/Users/quantm/Dropbox/JEV_GAVI_Gates_model/result_model/first_model/Sto_template/test3_16coun_200runid_every_countries/"
template_file_name = "stochastic_burden_template_JE-OUCRU-Clapham_"

#read from file:
Sto_para_template <-read.csv(paste(Sto_template_save_path, template_file_name, "-parameters.csv", sep = ""), check.names = F)
