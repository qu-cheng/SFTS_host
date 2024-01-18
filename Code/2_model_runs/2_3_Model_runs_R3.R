source("0_Library_and_functions.R")

Animal_par <- read.csv("2_model_runs/Animal_par.csv")

Study.ID <- read.csv("../Data/Seroprevalence_data_from_lit.csv") %>%
  select(SurveyID, References) %>% 
  unique() %>%
  mutate(References = gsub(" ", "_", References))

R2.names <- grep("_R2", list.files("2_model_runs/Results"), value = TRUE)
print(R2.names)


R3.surveyID <- c(1:9)[!grepl("_FINAL", R2.names)]    # need to run R3

print(R3.surveyID)

nc = parallel::detectCores()
print(nc)
options('mc.cores' = nc)

for(survey.i in R3.surveyID){
  print(survey.i)
  Animal_prev <- read.csv("../Data/Seroprevalence_data_from_lit.csv") %>%
    filter(SurveyID == survey.i) %>%
    rowwise() %>%
    mutate(p = n_positive/n,
           prev_low = binconf(n_positive, n, method = "wilson")[2],
           prev_high = binconf(n_positive, n, method = "wilson")[3],
           Host.species = factor(Host.species, levels = c("Goat/sheep", "Cattle", "Poultry", "Dog", "Pig", "Rodent", "Hedgehog", "Weasel", "Hare"))) %>%
    arrange(Host.species)
  
  n.species = nrow(Animal_prev)      # sheep, cattle, dogs, pigs, chickens, rodents
  
  #========================== Round 3 =========================
  # in the round of analyses, we first estimated the upper bound for chi_max * beta1 (denoted as CB_high here) from all passes in R2, then sampled chi_max from [0, 1.2*CB_high/beta1], after we sampled beta1 from its range in Round 2
  
  R2 <- read.csv(paste("2_model_runs/Results/Survey_", survey.i, "_", Study.ID$References[survey.i], "_R2.csv", sep = ""))
  
  beta.range <- select(R2, matches("beta[1-9]")) %>% apply(2, range) %>% t() %>% as.data.frame() %>% rename(beta.min = V1, beta.max = V2)
  Chi_max_high <- 1.2*CB_high_R3(R2)   # use the function beta_for_R2() to obtain the new ranges for betas
  
  # Rerun the model with a modified sampler
  results <- mclapply(1:nc, function(i){
    simulate.fun(n.species, 21000, i, sample.function3, keep.k = n.species)  # set the second par: the number of repetitions to a smaller value if you want to try the code
  }, mc.cores = nc)
  
  result.current<- do.call(rbind, results) %>%
    data.frame()
  
  result.sens.real <- result.current %>%
    apply(2, as.numeric) %>%
    as.data.frame() %>%
    select(-contains("mu_H"))
  
  pass.no <- get_pass_no(result.sens.real)
  
  if(pass.no >= 1000)
  {
    write.csv(rbind(result.sens.real, result.sens.real.2), paste("2_model_runs/Results/Survey_", survey.i, "_", Study.ID$References[survey.i], "_R3_FINAL.csv", sep = ""), row.names = FALSE)
  }

  if(pass.no < 1000)   # continue sampling till 1000
  {
    # Rerun the model with a modified sampler
    results.2 <- mclapply(1:nc, function(i){
      simulate.fun(n.species, (1000 - pass.no)/pass.no *1.2 * 21000, i, sample.function3, keep.k = n.species)  # set the second par: the number of repetitions to a smaller value if you want to try the code
    }, mc.cores = nc)
    
    result.current.2 <- do.call(rbind, results.2) %>%
      data.frame()
    
    result.sens.real.2 <- result.current.2 %>%
      apply(2, as.numeric) %>%
      as.data.frame() %>%
      select(-contains("mu_H"))
    
    write.csv(rbind(result.sens.real, result.sens.real.2), paste("2_model_runs/Results/Survey_", survey.i, "_", Study.ID$References[survey.i], "_R3_FINAL.csv", sep = ""), row.names = FALSE)
  }
}



