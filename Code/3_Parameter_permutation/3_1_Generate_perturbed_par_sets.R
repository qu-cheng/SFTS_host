#========= generate perturbed parameter sets ================
all.passes <- read.csv("2_model_runs/Results/Passes_1000.csv")    

delta <- c(seq(0.1, 0.9, 0.1), 2:10)
for(survey.i in 1:9)
{
  current.pass <- all.passes %>%
    filter(SurveyID == survey.i) %>%
    select(SurveyID, ID, HostName, mu_T, phi, chi_H, beta, gamma_H, sigma_H, v_H, mu_H) %>%
    pivot_wider(names_from = HostName,
                values_from = c(chi_H, beta, gamma_H, sigma_H, v_H, mu_H))
  
  change.dat <- expand.grid(ID = unique(current.pass$ID),
                            change.par = colnames(current.pass)[-c(1:2)],
                            change.value = delta)
  
  current.pass <- current.pass %>%
    left_join(change.dat) %>%
    as.data.table()
  
  all.change.par <- colnames(current.pass)[-c(1:2)]
  for(j in 1:length(all.change.par))
  {
    current.var <- all.change.par[j]
    current.pass[change.par == all.change.par[j], (current.var) := get(current.var) * change.value]
  }
  
  write.csv(current.pass, file = paste("3_Parameter_permutation/3_2_Perturbed_parameters/Survey_", survey.i, "_Sensitivity.csv", sep = ""), row.names = FALSE)
}
