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
  
  #============= Round 1 ==================
  ######## !!!!! WARNING: TAKES A LONG TIME ###########
  results <- mclapply(1:nc, function(i){
    simulate.fun(n.species, 21000, i, sample.function1, keep.k = 2)  # set the second par: the number of repetitions to a smaller value if you want to try the code; we used 21000 here to get 21000*48 = ~ 1 million simulations
  }, mc.cores = nc)
  
  result.current<- do.call(rbind, results) %>%
    data.frame()
  
  result.sens.real <- result.current %>%
    apply(2, as.numeric) %>%
    as.data.frame() %>%
    select(-contains("mu_H"))
  
  print(get_pass_no(result.sens.real)) 
  
  write.csv(result.sens.real, paste("5_Missing_host/Results/Survey_9_missing_", i,  "_R1.csv", sep = ""), row.names = FALSE)
}
