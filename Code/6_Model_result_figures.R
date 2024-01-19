source("0_Library_and_functions.R")

Animal_prev <- read.csv("../Data/Seroprevalence_data_from_lit.csv") %>%
  rowwise() %>%
  mutate(p = n_positive/n,
         prev_low = binconf(n_positive, n, method = "wilson")[2],
         prev_high = binconf(n_positive, n, method = "wilson")[3],
         Host.species = factor(Host.species, levels = c("Goat/sheep", "Cattle", "Poultry", "Dog", "Pig", "Rodent", "Hedgehog", "Weasel", "Hare"))) %>%
  arrange(SurveyID, Host.species)


#============================================================================
#                         Figure 3: Overall R0
#============================================================================
all.passes <- read.csv("2_model_runs/Results/Passes_1000.csv")                     # read in 1000 passes for each survey
passes.final.ID <- read.csv("2_model_runs/Results/Passes_1000_resampledID.csv")    # read the resampled weights

overall.R0 <- passes.final.ID %>%
  left_join(unique(all.passes[, c("R0", "SurveyID", "ID")])) %>%
  select(SurveyID, ID, R0) %>%
  mutate(SurveyID = as.factor(SurveyID),
         SurveyID = factor(SurveyID, levels = 9:1)) %>%
  group_by(SurveyID) %>%
  summarise(R0.median = median(R0),
            R0.low = quantile(R0, 0.025),
            R0.high = quantile(R0, 0.975)) %>%
  mutate(R0.format = paste(round(R0.median, 2), " (", round(R0.low, 2), ",", round(R0.high, 2),")", sep = ""),
         R0 = 3.5)   # for adding the text in the figure to the correct place


Fig3A <- passes.final.ID %>%
  left_join(all.passes %>%
              select(SurveyID, ID, R0) %>%
              unique()) %>%
  mutate(SurveyID = as.factor(SurveyID),
         SurveyID = factor(SurveyID, levels = 9:1)) %>%
  ggplot(aes(x = R0, y = SurveyID, fill = as.factor(SurveyID))) +
  geom_density_ridges(stat = "binline", bins = 100) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_y_discrete(breaks = 1:9, expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  guides(fill = FALSE) +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_log10() +
  xlab(expression("Overall R"[0])) +
  scale_fill_hue(direction = -1) +
  ylab("Survey ID") +
  geom_text(data = overall.R0, aes(x = R0, y = SurveyID, label = R0.format), nudge_y = 0.15, size = 3) 






#========== Overall R0 with animal prev ============
Fig3B <- Animal_prev %>%
  select(SurveyID, Host.species, prev, prev_low, prev_high) %>%
  left_join(passes.final.ID %>%
              left_join(all.passes %>%
                          select(SurveyID, ID, R0) %>%
                          unique()) %>%
              group_by(SurveyID) %>%
              summarise(R0.low = quantile(R0, 0.025), R0.high = quantile(R0, 0.975), R0 = median(R0))) %>%
  ggplot(aes(x = R0, y = prev, col = as.factor(SurveyID), label = SurveyID)) +
  geom_point() +
  geom_errorbar(aes(xmin = R0.low, xmax = R0.high)) +
  geom_errorbar(aes(ymin = prev_low*100, ymax = prev_high*100)) +
  #  geom_text_repel(col = "black", size = 4) +
  theme_cowplot() +
  facet_wrap(~Host.species, scales = "free") +
  ylab("Seroprevalence (%)") +
  labs(color = "Survey ID", shape = "Survey ID") +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  xlab(expression("Overall R"[0])) +
  scale_x_continuous(limits = c(0.9, 5.5)) +
  scale_y_continuous(limits = c(0, 100)) +
  guides(fill = FALSE, col = FALSE)

save_plot("../Figures/Fig3_overall_R0.pdf", plot_grid(Fig3A, Fig3B, rel_widths = c(0.45, 1.1), labels = c("A", "B"), label_size = 14), base_width = 9, base_height = 6)




#============================================================================
#                           Figure 4: Species-level R0i
#============================================================================
# estimate R0 for the larger gray circles
node.set1 <- passes.final.ID %>%
  left_join(all.passes %>%
              select(R0, ID, SurveyID) %>%
              unique()) %>%
  group_by(SurveyID) %>%
  summarise(R0 = median(R0)) %>%
  rename(name = SurveyID) %>%
  mutate(name = paste("Survey", name),
         R_host = "Overall")

node.set1.center.coord <- node.set1 %>%
  circleProgressiveLayout(sizecol = 2, sizetype='radius') %>%
  mutate(ID = 1:9, 
         SurveyID = paste("Survey", ID, sep = " "))


# estimate R0i for the colored circles
Host.count <- unique(all.passes %>%      # get the number of surveyed species for each survey
                       select(SurveyID, HostName)) %>%
  group_by(SurveyID) %>%
  tally()

HostName.ID <- unique(all.passes %>%    # get a lookup table for Host ID and HostName
                        select(SurveyID, HostName)) %>%
  group_by(SurveyID) %>%
  mutate(HostID = 1:n()) %>%
  ungroup()

node.set2.coord <- passes.final.ID %>%
  left_join(Host.count) %>%
  uncount(n, .id = "HostID") %>%
  left_join(HostName.ID) %>%
  select(-HostID) %>%
  ungroup() %>%
  left_join(all.passes) %>%
  group_by(SurveyID, HostName) %>%
  dplyr::summarize(R0 = median(R_i)) %>%
  ungroup() %>%
  mutate(name = paste(SurveyID, HostName, sep = "_")) %>%
  split(.$SurveyID) %>% 
  map(~circleProgressiveLayout(.x$R0, sizetype='radius') %>% mutate(R_host = .x$HostName)) %>% 
  imap_dfr(~circleLayoutVertices(.x, npoints=50, idcol = 4) %>% mutate(category = paste("Survey", .y, sep = " ")))

Location.all <- unique(node.set2.coord$category)

level2.all <- NULL
for(i in 1:length(Location.all))
{
  current.sets <- node.set2.coord %>%
    filter(category == Location.all[i])
  
  polygons <- df_to_SpatialPolygons(current.sets, keys = "category", c("x", "y"),CRS())
  poly.pad <- boundingcircle(as.owin(polygons))
  
  level2.center.point <- c(mean(poly.pad$xrange), mean(poly.pad$yrange))
  level2.radius <- diff(poly.pad$xrange)/2
  
  poly.pad.df <- data.frame(x = poly.pad$bdry[[1]]$x, y = poly.pad$bdry[[1]]$y)
  
  scale.factor <- node.set1.center.coord$radius[i]/level2.radius
  
  current.sets <- current.sets %>%
    mutate(x.new = (x - level2.center.point[1])*scale.factor + node.set1.center.coord$x[i],
           y.new = (y - level2.center.point[2])*scale.factor + node.set1.center.coord$y[i])
  
  # current.sets <- data.frame(x.new = level2.center.point[1] + node.set1.center.coord$x[i],
  #                            y.new = level2.center.point[2] + node.set1.center.coord$y[i],
  #                            id = i)
  
  level2.all <- rbind(level2.all, current.sets)
}


node.set1 %>%
  circleProgressiveLayout(sizecol = 2, sizetype='radius') %>%
  mutate(SurveyID = 1:9) %>%
  circleLayoutVertices(npoints = 50, idcol = 4) %>%
  ggplot() +
  ggplot2::geom_polygon(aes(x = x, y = y, group = id), fill = "gray") +
  theme_void() +
  coord_equal() +
  geom_polygon(data = level2.all, aes(x = x.new, y = y.new, group = interaction(category, id), fill = id), col = "black") +
  scale_fill_manual(breaks = c("Goat/sheep", "Cattle", "Poultry", "Dog", "Pig", "Rodent", "Hedgehog", "Weasel", "Hare"), values = ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value[1:9]) +
  geom_node_label(data = node.set1.center.coord, aes(x = x, y = y, label = SurveyID), alpha = 0.7, label.padding = unit(0.2, "lines"), size = 6) +
  labs(fill = "") +
  theme(text = element_text(size = 16))
ggsave("../Figures/Fig4_species_R0i.pdf", width = 8, height = 6)











#============================================================================
#                   Determinants of  R0i - figure
#============================================================================
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



# plot the scaling parameter against R0i
all.change.par <- unique(sens.all.location$change.par)
delta <- c(seq(0.1, 0.9, 0.1), 2:10)
for(i in 1:9)
{
  sens.all.location %>%
    filter(SurveyID == i) %>%
    group_by(SurveyID, change.par, change.value, Host) %>%
    dplyr::summarize(pass.prop = mean(pass), R0_host = median(R0_host[pass])) %>%
    filter(pass.prop > 0) %>%
    mutate(change.par = factor(change.par, levels = c(paste("mu_H_", c("Goat/sheep", "Cattle", "Poultry", "Dog", "Pig", "Rodent", "Hedgehog", "Weasel", "Hare"), sep = ""),
                                                      "mu_T", 
                                                      paste("beta_", c("Goat/sheep", "Cattle", "Poultry", "Dog", "Pig", "Rodent", "Hedgehog", "Weasel", "Hare"), sep = ""),
                                                      
                                                      paste("gamma_H_", c("Goat/sheep", "Cattle", "Poultry", "Dog", "Pig", "Rodent", "Hedgehog", "Weasel", "Hare"), sep = ""),
                                                      paste("v_H_", c("Goat/sheep", "Cattle", "Poultry", "Dog", "Pig", "Rodent", "Hedgehog", "Weasel", "Hare"), sep = ""),
                                                      paste("sigma_H_", c("Goat/sheep", "Cattle", "Poultry", "Dog", "Pig", "Rodent", "Hedgehog", "Weasel", "Hare"), sep = ""),
                                                      "phi",
                                                      paste("chi_H_", c("Goat/sheep", "Cattle", "Poultry", "Dog", "Pig", "Rodent", "Hedgehog", "Weasel", "Hare"), sep = "")),
                               
                               ordered = TRUE,
                               labels = c(expression(mu[Goat/sheep]), expression(mu[Cattle]), expression(mu[Poultry]), expression(mu[Dog]), expression(mu[Pig]), expression(mu[Rodent]), expression(mu[Hedgehog]), expression(mu[Weasel]), expression(mu[Hare]),
                                          expression(mu["T"]), 
                                          expression(beta[Goat/sheep]), expression(beta[Cattle]), expression(beta[Poultry]), expression(beta[Dog]), expression(beta[Pig]), expression(beta[Rodent]), expression(beta[Hedgehog]), expression(beta[Weasel]), expression(beta[Hare]),
                                          expression(gamma[Goat/sheep]), expression(gamma[Cattle]), expression(gamma[Poultry]), expression(gamma[Dog]), expression(gamma[Pig]), expression(gamma[Rodent]), expression(gamma[Hedgehog]), expression(gamma[Weasel]), expression(gamma[Hare]),
                                          expression(nu[Goat/sheep]), expression(nu[Cattle]), expression(nu[Poultry]), expression(nu[Dog]), expression(nu[Pig]), expression(nu[Rodent]), expression(nu[Hedgehog]), expression(nu[Weasel]), expression(nu[Hare]),
                                          expression(sigma[Goat/sheep]), expression(sigma[Cattle]), expression(sigma[Poultry]), expression(sigma[Dog]), expression(sigma[Pig]), expression(sigma[Rodent]), expression(sigma[Hedgehog]), expression(sigma[Weasel]), expression(sigma[Hare]),
                                          expression(phi),
                                          expression(chi[Goat/sheep]), expression(chi[Cattle]), expression(chi[Poultry]), expression(chi[Dog]), expression(chi[Pig]), expression(chi[Rodent]), expression(chi[Hedgehog]), expression(chi[Weasel]), expression(chi[Hare])))) %>%
    ggplot(aes(x = change.value, y = R0_host, size = pass.prop, col = Host))+
    #  geom_vline(xintercept = delta, col = "gray80") +
    geom_point() +
    facet_wrap(~change.par, scales = "free", labeller = label_parsed) +
    theme_cowplot() +
    ggtitle(paste("Survey", i)) +
    xlab("Scaling factor for the parameter value") +
    ylab(expression("R"["0i"])) +
    labs(size = "Passing\nproportion") +
    theme(strip.background = element_blank()) +
    scale_color_manual(breaks = c("Goat.sheep", "Cattle", "Poultry", "Dog", "Pig", "Rodent", "Hedgehog", "Weasel", "Hare"), labels = c("Goat/sheep", "Cattle", "Poultry", "Dog", "Pig", "Rodent", "Hedgehog", "Weasel", "Hare"), values = ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value[1:9])+
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    scale_x_continuous(breaks = scales::pretty_breaks(4), limits = c(0, 10))
  
  
  ggsave(paste("../Figures/FigS3_sensitivity_survey", i, ".pdf", sep = ""), width = 14, height = 10)
}









#============================================================================
#                   Determinants of  R0i - RF figure
#============================================================================
RF.result <- read.csv("3_Parameter_permutation/3_6_RF_importance.csv")

Fig.sens <- RF.result %>%
  group_by(variable) %>%
  summarise(importance.median = median(importance), low = quantile(importance, 0.025), high = quantile(importance, 0.975)) %>%
  ggplot(
    aes(
      x = reorder(variable, importance.median),
      y = importance.median,
      fill = importance.median
    )
  ) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.15) +
  coord_flip() +
  ylab("Permutation importance score") +
  xlab("") +
  guides(fill = "none") +
  theme_cowplot() +
  scale_fill_continuous_tableau() +
  scale_x_discrete(labels = c("beta" = expression(beta[i]),
                              "chi_H" = expression(chi[i]),
                              "phi" = expression(phi),
                              "mu_H" = expression(mu[i]),
                              "sigma_H" = expression(sigma[i]),
                              "mu_T" = expression(mu[T]),
                              "gamma_H" = expression(gamma[i]),
                              "v_H" = expression(nu[i])))

save_plot("../Figures/Fig5_RF_importance.pdf", Fig.sens, base_height = 4, base_width = 8)






#============================================================================
#                          Optimal interventions
#============================================================================
all.result <- read.csv("4_Optimal_intervention/4_2_R0_decrease.csv")

# unify the range
all.rank <- all.result %>%
  mutate(change.value = ifelse(change.par == "muT", 1/change.value, change.value),
         change.value = ifelse(grepl("sigma", change.par), 1/change.value, change.value)) %>%
  filter(change.value == 0.1) %>%
  group_by(SurveyID) %>%
  mutate(rank = order(order(new.R0))) %>%
  select(-change.value, -new.R0)



change.par.include <- all.rank %>%
  filter(rank <= 5) %>%
  pull(change.par) %>%
  unique()

all.result %>%
  mutate(change.value = ifelse(change.par == "muT", 1/change.value, change.value),
         change.value = ifelse(grepl("sigma", change.par), 1/change.value, change.value))%>%
  left_join(all.rank) %>% 
  filter(rank <= 5) %>%
  mutate(change.par = ifelse(change.par == "muT", "1/mu[T]", change.par),
         change.par = ifelse(change.par == "beta_Goat/sheep", "beta[Goat/sheep]", change.par),
         change.par = ifelse(change.par == "beta_Poultry", "beta[Poultry]", change.par),
         change.par = ifelse(change.par == "beta_Rodent", "beta[Rodent]", change.par),
         change.par = ifelse(change.par == "beta_Cattle", "beta[Cattle]", change.par),
         change.par = ifelse(change.par == "chi_H_Goat/sheep", "chi[Goat/sheep]", change.par),
         change.par = ifelse(change.par == "chi_H_Poultry", "chi[Poultry]", change.par),
         change.par = ifelse(change.par == "chi_H_Rodent", "chi[Rodent]", change.par),
         change.par = ifelse(change.par == "chi_H_Cattle", "chi[Cattle]", change.par),
         change.par = ifelse(change.par == "sigma_H_Goat/sheep", "1/sigma[Goat/sheep]", change.par),
         change.par = ifelse(change.par == "sigma_H_Poultry", "1/sigma[Poultry]", change.par),
         change.par = ifelse(change.par == "sigma_H_Rodent", "1/sigma[Rodent]", change.par),
         change.par = ifelse(change.par == "sigma_H_Cattle", "1/sigma[Cattle]", change.par),
         greeklabel = ifelse(change.value == 0.1, change.par, NA)) %>%
  ggplot(aes(x = change.value, y = new.R0, col = change.par, label = greeklabel)) +
  geom_line() +
  geom_text_repel(parse = TRUE, force_pull = 10, nudge_x = 0.1, nudge_y = 0.1, size = 3) +
  facet_wrap(~SurveyID) +
  scale_x_reverse(breaks = seq(0, 1, 0.2), limits = c(1, -0.15)) +
  theme_cowplot() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  guides(col = "none") +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, face = "bold")) +
  ylab(expression("Overall R"[0])) +
  xlab("Scaling parameter") +
  scale_color_tableau(palette = "Tableau 20")
ggsave("../Figures/FigS4_optimal_intervention.pdf", width = 8, height = 7.5)










#============================================================================
#                              Missing host
#============================================================================
Animal_par <- read.csv("2_model_runs/Animal_par.csv")
all.passes <- read.csv("2_model_runs/Results/Passes_1000.csv")                     # read in 1000 passes for each survey
passes.final.ID <- read.csv("2_model_runs/Results/Passes_1000_resampledID.csv")    # read the resampled weights

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


# get the R0s when using the full dataset
original.species.R0 <- passes.final.ID %>%
  filter(SurveyID == 9) %>%
  left_join(all.passes %>%
              filter(SurveyID == 9) %>%
              select(ID, R0, R_i, HostName) %>%
              spread(HostName, R_i)) %>%
  select(-SurveyID) %>%
  gather(HostName, "R0", R0:Weasel) %>%
  group_by(HostName) %>%
  summarise(R0.median = median(R0), R0.low = quantile(R0, 0.025), R0.high = quantile(R0, 0.975)) %>%
  mutate(Without.Host = "Full dataset",
         HostName = ifelse(HostName == "R0", "Overall", HostName))



# get the R0s when missing each host species
file.names <- grep("_FINAL", list.files("5_Missing_host/Results"), value = TRUE)

all.result <- NULL
for(i in 1:9)
{
  Animal_prev <- Animal_prev_all %>%
    ungroup() %>%
    slice(-i) %>%
    mutate(Host.species = as.character(Host.species))
  
  current.species <- data.frame(Host = paste("R", 1:8, sep = ""), Host.species = Animal_prev$Host.species)
  
  n.species = nrow(Animal_prev) 
  
  current.result <- read.csv(paste("5_Missing_host/Results/", file.names[i], sep = ""))
  current.pass <- current.result[get_pass(current.result, Animal_prev),] %>%
    mutate(ID = 1:n()) %>%
    filter(ID <= 1000)   # take the first 1000 passes
  
  # resample with the likkelihood
  passes.ll <- current.pass  %>%
    select(ID, R1:R8) %>%
    gather("Host", "R", R1:R8) %>%
    left_join(current.species) %>%
    left_join(Animal_prev[, c("Host.species", "n_positive", "n")]) %>%
    mutate(ll = dbinom(n_positive, prob = R, size = n)) %>%
    group_by(ID) %>%
    summarise(ll = prod(ll))
  
  current.samples <- data.frame(ID = sample(passes.ll$ID, 10000, prob = passes.ll$ll, replace = TRUE))
  Host.name.lookup <- data.frame(HostName.ID = c("R0", paste("R_i", 1:8, sep = "")), HostName = c("R0", Animal_prev$Host.species))
  
  # get summary statistics for the R0s
  current.result <- current.samples %>%
    left_join(current.pass %>%
                select(ID, R0:R_i8)) %>%
    gather(HostName.ID, "R0", R0:R_i8) %>%
    group_by(HostName.ID) %>%
    summarise(R0.median = median(R0), R0.low = quantile(R0, 0.025), R0.high = quantile(R0, 0.975)) %>%
    ungroup() %>%
    mutate(Without.Host = paste("w/o", tolower(Animal_prev_all$Host.species[i]))) %>%
    left_join(Host.name.lookup) %>%
    select(-HostName.ID)
  
  all.result <- rbind(all.result, current.result)
}



all.result <- all.result %>%
  mutate(HostName = ifelse(HostName == "R0", "Overall", HostName)) %>%
  rbind(original.species.R0) %>%
  mutate(Without.Host = factor(Without.Host, levels = c("Full dataset", "w/o goat/sheep", "w/o cattle", "w/o poultry", "w/o dog", "w/o pig", "w/o rodent", "w/o hedgehog", "w/o weasel", "w/o hare")),
         HostName = factor(HostName, levels = c("Overall", "Goat/sheep", "Cattle", "Poultry", "Dog", "Pig", "Rodent", "Hedgehog", "Weasel", "Hare"))) 


Fig6A <- all.result %>%
  filter(HostName == "Overall") %>%
  ggplot(aes(x = HostName, y = R0.median, ymin = R0.low, ymax = R0.high, col = Without.Host)) +
  geom_point(position = position_dodge(0.6)) +
  geom_errorbar(width = 0.2, position = position_dodge(0.6)) +
  theme_cowplot() +
  scale_color_manual(values = c("black", ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value[1:9])) +
#  xlab(expression("Overall "*R[0])) +
  ylab(expression(R[0])) +
  xlab("") +
  guides(col = FALSE) +
  ggtitle(expression("A.Overall "*R[0])) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, face = "bold"),
        legend.position = c(0.27, 0.95))+
  ylim(0, 2.1)


Fig6B <- all.result %>%
  filter(HostName != "Overall") %>%
  ggplot(aes(x = HostName, y = R0.median, ymin = R0.low, ymax = R0.high, col = Without.Host)) +
  geom_point(position = position_dodge(0.6)) +
  geom_errorbar(width = 0.2, position = position_dodge(0.6)) +
  theme_cowplot() +
  scale_color_manual(values = c("black", ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value[1:9])) +
  #  xlab(expression("Overall "*R[0])) +
  ylab(expression(R["0i"])) +
  xlab(expression(R["0i"]*" for")) +
  ggtitle(expression("B.Species-level "*R["0i"])) +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, face = "bold"),
        legend.position = c(0.05, 0.95))+
  guides(col = guide_legend(ncol = 5)) +
  labs(col = "") +
  ylim(0, 2.1)

save_plot("../Figures/Fig6_missing_host.pdf", plot_grid(Fig6A, Fig6B, rel_widths = c(0.2, 0.9)), base_width = 10, base_height = 5)
