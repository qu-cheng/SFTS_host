source("0_Library_and_functions.R")

Animal_par <- read.csv("2_model_runs/Animal_par.csv")

Animal_prev <- read.csv("../Data/Seroprevalence_data_from_lit.csv") %>%
  rowwise() %>%
  mutate(p = n_positive/n,
         prev_low = binconf(n_positive, n, method = "wilson")[2],
         prev_high = binconf(n_positive, n, method = "wilson")[3],
         Host.species = factor(Host.species, levels = c("Goat/sheep", "Cattle", "Poultry", "Dog", "Pig", "Rodent", "Hedgehog", "Weasel", "Hare"))) %>%
  arrange(SurveyID, Host.species)

final.names <- grep("_FINAL", list.files("2_model_runs/Results"), value = TRUE)
print(final.names)





#===================== take the first 1000 passes for each survey ===========================
all.passes <- NULL
for(survey.i in 1:9)
{
  print(survey.i)
  current.result <- read.csv(paste("2_model_runs/Results/", final.names[survey.i], sep = ""))   # read final result
  
  current.Animal_prev <- Animal_prev %>%
    filter(SurveyID == survey.i) 
  n.species <- nrow(current.Animal_prev)
  
  pass.ID <- get_pass(current.result, current.Animal_prev)
  Location.pass <- current.result[pass.ID,]    # get passes
  
  Location.pass.long <- Location.pass %>%
    mutate(ID = 1:n()) %>%
    pivot_longer(cols = !c(mu_T, phi, R0, ID),
                 names_to = c(".value", "Host"),
                 names_pattern = "(.*)(.)") %>%
    left_join(data.frame(Host = as.character(1:n.species),
                         HostName = current.Animal_prev$Host.species)) %>%
    dplyr::select(-Host, -alpha, -rho_H) %>%
    left_join(Animal_par %>%
                rename(HostName = Animal) %>%
                mutate(mu_H = 1/(life_span_mean*365)) %>%
                select(HostName, mu_H)) %>%
    mutate(SurveyID = survey.i) %>%
    as.data.table()
  
  all.passes <- rbindlist(list(all.passes, Location.pass.long)) 
}


all.passes %>%
  filter(ID <= 1000) %>%
  write.csv("2_model_runs/Results/Passes_1000.csv", row.names = FALSE)   # save the first 1000 passes







#===================== resample with the likelihood weighting ===========================
set.seed(12345678)
all.passes <- read.csv("2_model_runs/Results/Passes_1000.csv") %>%
  left_join(Animal_prev %>%
              select(SurveyID, Host.species, n_positive, n) %>%
              rename(HostName = Host.species)) %>%
  mutate(ll = dbinom(n_positive, prob = R, size = n))

passes.ll <- all.passes %>%
  group_by(SurveyID, ID) %>%
  summarise(ll = prod(ll))

# sample with the likelihood
passes.final.ID <- NULL
for(i in 1:9)
{
  current.passes <- passes.ll %>%
    filter(SurveyID == i)
  
  current.samples <- data.frame(SurveyID = i, ID = sample(current.passes$ID, 10000, prob = current.passes$ll, replace = TRUE))
  
  passes.final.ID <- rbind(passes.final.ID, current.samples)
}

passes.final.ID %>%
  write.csv("2_model_runs/Results/Passes_1000_resampledID.csv", row.names = FALSE)   # save the resampled ID


Host.count <- unique(all.passes %>%      # get the number of surveyed species for each survey
                       select(SurveyID, HostName)) %>%
  group_by(SurveyID) %>%
  tally()


HostName.ID <- unique(all.passes %>% 
                        select(SurveyID, HostName)) %>%
  group_by(SurveyID) %>%
  mutate(HostID = 1:n()) %>%
  ungroup()


passes.final.ID %>%
  left_join(Host.count) %>%
  uncount(n, .id = "HostID") %>%
  left_join(HostName.ID) %>%
  select(-HostID) %>%
  ungroup() %>%
  left_join(all.passes) %>%
  #  mutate(HostName = factor(HostName, levels = c("Goat/sheep", "Cattle", "Poultry", "Dog", "Pig", "Rodent", "Hedgehog", "Weasel", "Hare"))) %>%
  group_by(SurveyID, HostName) %>%
  dplyr::summarize(R0.median = median(R0), R0.low = quantile(R0, 0.025), R0.high = quantile(R0, 0.975),
                   R0_i = median(R_i), R0_i.low = quantile(R_i, 0.025), R0_i.high = quantile(R_i, 0.975),
                   R0.report = paste(round(R0.median, 2), " (", round(R0.low, 2), ", ", round(R0.high, 2), ")", sep = ""),
                   R_i.report = paste(round(R0_i, 2), " (", round(R0_i.low, 2), ", ", round(R0_i.high, 2), ")", sep = "")) %>%
  ungroup() %>%
  select(SurveyID, R0.report, HostName, R_i.report) %>%
  write.csv("2_model_runs/Results/Resampled_species_R0i.csv", row.names = FALSE)
