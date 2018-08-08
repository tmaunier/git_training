#JEV model: Dirichlet Multinomial + Poisson to estimate case reported , based on Natsuko's work => not work yet
#library loading:
library(rstan)
library(bayesplot)

#Model for "A hospital-based surveillance for Japanese encephalitis in Bali, Indonesia"
Data_je_indo = data.frame(c("0-4", "5-9", "10-11"),
                          c(0, 5, 10),
                          c(4, 9, 11),
                          c(63, 26, 1), 
                          c(275500, 226700, 96920))
colnames(Data_je_indo) = c("age_group", "age_l", "age_u", "cases", "pop")

ll_case_je_indo="
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
  real<lower = 1> psi; //over-dispersion par from DMN
  real<lower = 0> rho; //reporting rate

}

transformed parameters{
  vector[N] I_age; //the incidence of infections in each age group
  vector[N] D_age; //The average annual incidence rate per person in an age group
  vector[N] C_age; //The expected number of cases per year each age group
  vector[N] p_age; //The expected proportion of cases in one age group relative to the total number of cases across all age groups

  I_age = exp(-lambda*age_l) - exp(-lambda*(age_u + 1));
  D_age = rho*(I_age ./ (age_u + 1 - age_l));
  C_age = pop_age .* D_age;
  p_age = C_age/sum(C_age);
}

model {
  real l_DMN;
  //prior distribution
  psi ~ normal(0, 1000);
  lambda ~ normal(0, 1000);
  rho ~ uniform(0,1);
  
  //DMN likelihood function:
  l_DMN = - (lgamma(1.0/psi + t_case) - lgamma(1.0/psi)) + sum(lgamma(p_age/psi + case_age) - lgamma(p_age/psi));
  
  //likelihood function, included poisson for total cases across all age group:
  target += l_DMN + t_case*log(sum(C_age)) - sum(C_age) - lgamma(t_case + 1);
}
"

stan_sero_model_2 = stan_model(model_code = ll_case_je_indo)


data = Data_je_indo
data_for_HMC = list(N = nrow(data), age_l = data$age_l, age_u = data$age_u, 
                    case_age = data$cases, pop_age = data$pop, 
                    t_case = sum(data$cases))

stan_FOI_fit = sampling(object = stan_sero_model_2, data = data_for_HMC, 
                        chains = 4, iter = 4000)

print(stan_FOI_fit)

mcmc_areas(as.array(stan_FOI_fit),regex_pars = "p_age")

mcmc_trace(as.array(stan_FOI_fit))

mcmc_acf(as.array(stan_FOI_fit))

mcmc_pairs(as.array(stan_FOI_fit),regex_pars = "p_age")
