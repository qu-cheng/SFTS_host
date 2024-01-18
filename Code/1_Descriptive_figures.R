source("0_Library_and_functions.R")

#================== Figure 2 ===========================
#===== map ======
register_stadiamaps("68605a8f-cdbf-4488-a1c7-a865c1965715", write = FALSE)
map <- get_stadiamap(c(bottom = 30, left = 112.5, top = 38, right = 123.5), source="stadia", col = "bw", zoom = 7)

Location.dat <- read.csv("../Data/StudyLocation.csv") %>%
  mutate(Ref = factor(Ref, levels = unique(Ref)))
Location.dat <- Location.dat[sample(1:nrow(Location.dat), nrow(Location.dat)),]   # to introduce a random pattern between those identical locations and to avoid 


jitter.width = 0.7
Fig2A <- ggmap(map) + 
  geom_quasirandom(data = Location.dat, aes(x = Long, y = Lat, col = Ref), size = 5, width = jitter.width) +
  geom_text(data = Location.dat, aes(x = Long, y = Lat, label = LocationID), col = "white",
            position = position_quasirandom(width = jitter.width), size = 2.5) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(text = element_text(size = 10),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        legend.position = "top",
        legend.text = element_text(size = 8), 
        legend.background = element_blank(),
        legend.key=element_blank(), 
        legend.spacing.x = unit(-0.05, "cm"),
        legend.box.spacing = unit(0, "cm")) +
  labs(col = "") +
  guides(col = guide_legend(override.aes = list(size = 2), byrow = TRUE, keyheight = unit(0.2, "cm")))


#===== Heatmap  =======
prev.data <- read.csv("../Data/Seroprevalence_data_from_lit.csv") %>%
  rowwise() %>%
  mutate(p = n_positive/n,
         prev_low = binconf(n_positive, n, method = "wilson")[2],
         prev_high = binconf(n_positive, n, method = "wilson")[3])

table(prev.data$Host.species)


# heatmap version
Fig2B <- prev.data %>%
  mutate(Host.species = factor(Host.species, levels = (c("Goat/sheep", "Cattle", "Poultry", "Dog", "Pig", "Rodent", "Hedgehog", "Weasel", "Hare"))),
         Survey.ID = as.factor(SurveyID),
         Survey.ID = factor(Survey.ID, levels = 9:1)) %>%
  ggplot(aes(x = Host.species, y = Survey.ID, fill = prev)) +
  geom_tile(col = "white", linewidth = 0.8) +
  theme_cowplot() +
  #  coord_fixed() +
  ylab("Survey ID") +
  xlab("") +
  scale_fill_continuous_tableau() +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0), position = "top") +
  labs(fill = "Seroprevalence rate (%)") +
  theme(text = element_text(size = 10),
        legend.position = "bottom",
        panel.border = element_rect(color = "black"),
        legend.key.width = unit(0.6, "cm"),
        legend.key.height = unit(0.2, "cm"),
        axis.line = element_line(color = NA), 
        panel.background = element_rect(fill = "gray70"),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 0, size = 9),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) 


save_plot("../Figures/Fig2.pdf",plot_grid(Fig2A, Fig2B, rel_widths = c(0.5, 0.4), labels = c("A", "B"), label_size = 14), base_width = 8, base_height = 4.3)





#================== Figure S2 ===========================
read.csv("../Data/Species_N_Npos.csv") %>%
  mutate(prev = Npos/N*100,
         prev.low = binconf(Npos, N)[,2]*100,
         prev.high = binconf(Npos, N)[,3]*100) %>%
  arrange(-prev) %>%
  mutate(Species = factor(Species, levels = Species)) %>%
  ggplot(aes(x = Species, y = prev)) +
  geom_errorbar(aes(ymin = prev.low, ymax = prev.high), width = 0.2) +
  geom_point(aes(col = as.factor(Number.of.studies), size = N)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_viridis_d() +
  labs(color = "No. of studies", size = "Sample size", y = "Prevalence rate (%)")
ggsave("../Figures/FigS2.pdf", width = 7, height = 5)








