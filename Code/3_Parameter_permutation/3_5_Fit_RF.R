Animal_prev <- read.csv("../Data/Seroprevalence_data_from_lit.csv") %>%
  rowwise() %>%
  mutate(p = n_positive/n,
         prev_low = binconf(n_positive, n, method = "wilson")[2],
         prev_high = binconf(n_positive, n, method = "wilson")[3],
         Host.species = factor(Host.species, levels = c("Goat/sheep", "Cattle", "Poultry", "Dog", "Pig", "Rodent", "Hedgehog", "Weasel", "Hare"))) %>%
  arrange(SurveyID, Host.species)

# paste all sensitivity results together
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




#=== random forest permutation importance ====
sens.all.location.pass <- sens.all.location %>%
  filter(pass)

RF.result <- NULL   # takes more than an hour
for(i in 1:100)
{
  print(paste(i, Sys.time()))
  
  rg.model.pass <- ranger(R0_host ~ ., sens.all.location.pass[sample(1:nrow(sens.all.location.pass), 200000), c("mu_T", "phi", "chi_H", "beta", "gamma_H", "sigma_H", "mu_H", "v_H", "R0_host")], importance = "permutation")
  
  current.result <- enframe(
    rg.model.pass$variable.importance,
    name = "variable",
    value = "importance"
  ) %>%
    mutate(rep = i)
  
  RF.result <- rbind(RF.result, current.result)
}

write.csv(RF.result, "3_Parameter_permutation/3_6_RF_importance.csv", row.names = FALSE)