source("0_Library_and_functions.R")

Animal_prev <- read.csv("../Data/Seroprevalence_data_from_lit.csv") %>%
  rowwise() %>%
  mutate(p = n_positive/n,
         prev_low = binconf(n_positive, n, method = "wilson")[2],
         prev_high = binconf(n_positive, n, method = "wilson")[3],
         Host.species = factor(Host.species, levels = c("Goat/sheep", "Cattle", "Poultry", "Dog", "Pig", "Rodent", "Hedgehog", "Weasel", "Hare"))) %>%
  arrange(SurveyID, Host.species)

all.passes <- read.csv("2_model_runs/Results/Passes_1000.csv")  

# paste all sensitivity results together, reuse the results from parameter permutation, but need to estimate R0, since the previous result only has R0i
sens.all.location <- NULL
for(survey.i in 1:9){
  print(survey.i)
  
  prev.dat <- Animal_prev %>%
    filter(SurveyID == survey.i)
  
  n.species <- nrow(prev.dat)
  sens.results <- read_csv(paste("3_Parameter_permutation/3_4_Perturbed_results/Survey_", survey.i, "_Sensitivity_results.csv", sep = ""))
  
  sens.all <- read.csv(paste("3_Parameter_permutation/3_2_Perturbed_parameters/Survey_", survey.i, "_Sensitivity.csv", sep = "")) %>%
    mutate(pass = apply(sens.results, 1, function(x){sum(x <= prev.dat$prev_high & x >= prev.dat$prev_low) == n.species})) %>%
    pivot_longer(cols = !c(SurveyID, ID, mu_T, phi, change.par, change.value, pass),
                 names_to = c(".value", "Host"),
                 names_pattern = "(.*)_(.*)") %>%
    mutate(rho_H = v_H/(1-v_H) * (mu_H + sigma_H),
           kk = 4*beta^2*chi_H/mu_T*gamma_H/(mu_H + gamma_H)/(mu_H + rho_H + sigma_H),
           R0_host = (phi + sqrt(phi^2 + kk))/2) 
  
  sens.all.location <- rbind(sens.all.location, sens.all)
}

# filter out parameters related to common interentions
interv.dat <- sens.all.location %>%   # time-consuming due to the size of sens.all.location
  filter((change.par == "phi" & change.value <= 1) | (grepl("beta", change.par) & change.value <= 1) | (grepl("chi_H", change.par) & change.value <= 1) | (grepl("sigma_H", change.par) & change.value >= 1)) 

interv.dat  <- interv.dat %>%    # estimate R0
  group_by(SurveyID, ID, change.par, change.value) %>%
  mutate(R0 = (phi + sqrt(phi^2 + sum(kk)))/2)

# generate a lookup table for the intervention-related parameters for each survey
Location.Host <- interv.dat %>%
  ungroup() %>%
  select(SurveyID, change.par)%>%
  filter(grepl("beta|chi_H|sigma_H|phi", change.par))  %>%
  unique()

interv.dat.median <- interv.dat %>%     # time-consuming due to the size of sens.all.location
  select(SurveyID, ID, change.par, change.value, R0) %>%
  unique() %>%
  group_by(SurveyID, change.value, change.par) %>%
  summarise(new.R0 = median(R0)) %>%
  rbind(all.passes %>%
          select(SurveyID, R0) %>% 
          unique() %>%
          group_by(SurveyID) %>%
          summarise(new.R0 = median(R0)) %>%
          mutate(change.value = 1) %>%
          left_join(Location.Host)) %>%
  mutate(SurveyID = paste("Survey ", SurveyID, sep = "")) 



# increase mu_T, but not birth rate of ticks, need to reestimate R0 since the results were not run for the simulation
delta <- c(1:10)
all.dat <- NULL
for(survey.i in 1:9)
{
  current.pass <- all.passes %>%
    filter(SurveyID == survey.i) %>%
    select(SurveyID, ID, HostName, mu_T, phi, chi_H, beta, gamma_H, sigma_H, v_H, mu_H) %>%
    pivot_wider(names_from = HostName,
                values_from = c(chi_H, beta, gamma_H, sigma_H, v_H, mu_H))
  
  change.dat <- expand.grid(ID = unique(current.pass$ID),
                            change.value = delta)
  
  current.pass <- current.pass %>%
    left_join(change.dat) %>%
    mutate(new.mu_T = mu_T*change.value) %>%
    pivot_longer(cols = !c(SurveyID, ID, mu_T, new.mu_T, phi, change.value),
                 names_to = c(".value", "Host"),
                 names_pattern = "(.*)_(.*)") %>%
    mutate(rho_H = v_H/(1-v_H) * (mu_H + sigma_H),
           kk = 4*beta^2*chi_H/new.mu_T*gamma_H/(mu_H + gamma_H)/(mu_H + rho_H + sigma_H),
           R0_host = (phi*mu_T/new.mu_T + sqrt((phi*mu_T/new.mu_T)^2 + kk))/2) %>%
    group_by(SurveyID, ID, change.value) %>%
    mutate(R0 = (phi*mu_T/new.mu_T + sqrt((phi*mu_T/new.mu_T)^2 + sum(kk)))/2)
  
  all.dat  <- rbind(all.dat , current.pass)
}


muT.result.median <- all.dat %>% 
  select(SurveyID, ID, change.value, R0) %>%
  unique() %>%
  group_by(SurveyID, change.value) %>%
  summarise(new.R0 = median(R0)) %>%
  mutate(SurveyID = paste("Survey ", SurveyID, sep = ""),
         change.par = "muT")



all.result <- rbind(interv.dat.median,
                    muT.result.median)

write.csv(all.result, "4_Optimal_intervention/4_2_R0_decrease.csv", row.names = FALSE)
