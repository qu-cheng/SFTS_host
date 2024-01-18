source("0_Library_and_functions.R")

Animal_par <- read.csv("2_model_runs/Animal_par.csv")

Study.ID <- read.csv("../Data/Seroprevalence_data_from_lit.csv") %>%
  select(SurveyID, References) %>% 
  unique() %>%
  mutate(References = gsub(" ", "_", References))

Animal_prev_all <- read.csv("../Data/Seroprevalence_data_from_lit.csv") %>%
  filter(SurveyID == 9) %>%
  rowwise() %>%
  mutate(p = n_positive/n,
         prev_low = binconf(n_positive, n, method = "wilson")[2],
         prev_high = binconf(n_positive, n, method = "wilson")[3],
         Host.species = factor(Host.species, levels = c("Goat/sheep", "Cattle", "Poultry", "Dog", "Pig", "Rodent", "Hedgehog", "Weasel", "Hare"))) %>%
  arrange(Host.species)

n.species.all = nrow(Animal_prev_all)      # sheep, cattle, dogs, pigs, chickens, rodents

nc = parallel::detectCores()
print(nc)
options('mc.cores' = nc)


for(i in 1:n.species.all)
{
  Animal_prev <- Animal_prev_all %>%
    ungroup() %>%
    slice(-i)
  
  n.species = nrow(Animal_prev)      # sheep, cattle, dogs, pigs, chickens, rodents
  
  survey.i <- i
  
  #============= Round 2 ==================
  ######## !!!!! WARNING: TAKES A LONG TIME ###########
  #========================== Round 2 =========================
  # in the round of analyses, the ranges of betas were cut
  # stopped at R2 if more than 1000 passes were collected, proceeded to Round 3 if more than 50 but less than 1000 passes were collected, repeated R2 if less than 50 passes were collected till we got 50 passes
  R1 <- read.csv(paste("5_Missing_host/Results/Survey_9_missing_", i,  "_R1.csv", sep = ""))
  
  print(nrow(R1))
  
  beta.range <- beta_for_R2(R1, Animal_prev, 0, 2)   # use the function beta_for_R2() to obtain the new ranges for betas
  print(beta.range)
  
  # Rerun the model with a modified sampler
  results <- mclapply(1:nc, function(i){
    simulate.fun(n.species, 21000, i, sample.function2, keep.k = n.species-2)  # set the second par: the number of repetitions to a smaller value if you want to try the code
  }, mc.cores = nc)
  
  result.current<- do.call(rbind, results) %>%
    data.frame()
  
  result.sens.real <- result.current %>%
    apply(2, as.numeric) %>%
    as.data.frame() %>%
    select(-contains("mu_H"))
  
  pass.no <- get_pass_no(result.sens.real)
  
  if(pass.no >= 1000)   # save the result as FINAL if having more than 1000 passes
  {
    write.csv(result.sens.real, paste("5_Missing_host/Results/Survey_9_missing_", i,  "_R2_FINAL.csv", sep = ""), row.names = FALSE)
  }
  
  if(pass.no < 1000)   # continue sampling till 1000
  {
    # Rerun the model with a modified sampler
    results.2 <- mclapply(1:nc, function(i){
      simulate.fun(n.species, (1000 - pass.no)/pass.no *1.2 * 21000, i, sample.function2, keep.k = n.species)  # set the second par: the number of repetitions to a smaller value if you want to try the code
    }, mc.cores = nc)
    
    result.current.2 <- do.call(rbind, results.2) %>%
      data.frame()
    
    result.sens.real.2 <- result.current.2 %>%
      apply(2, as.numeric) %>%
      as.data.frame() %>%
      select(-contains("mu_H"))
    
    write.csv(rbind(result.sens.real, result.sens.real.2), paste("5_Missing_host/Results/Survey_9_missing_", i,  "_R2_FINAL.csv", sep = ""), row.names = FALSE)
  }
}