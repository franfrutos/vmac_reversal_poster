source("scripts/functions.R") # This will load functions and data into the workspace

if (!require(pacman)) {
  install.packages(pacman)
  library(pacman)
}

p_load(dplyr,# Data wrangling
       lme4,# LMMs
       lmerTest, # P-values in summary
       marginaleffects,# Model predictions and conditional effects
       hypr,# Set the contrast matrix
       sjPlot,# Tables
       afex,# ANOVA
       Rmisc,# Function to averaged the data in within-subject designs
       #effectsize, # Cohen's d
       ggplot2, # Plotting
       papaja
)

sum_RTs <- filter_data(raw_all, f.absent = F, sd_filter = NULL, fixed = F) %>%
  filter(!Phase %in% c("Awareness", "Practice")) %>%
  dplyr::summarise(RT = mean(rt), .by = c(Experiment, ID))

ggplot(sum_RTs, aes(RT, color = Experiment, fill = Experiment)) +
  geom_histogram(aes(y = ..density..), alpha = .3, position = position_dodge()) +
  geom_density(alpha = .3)+
  geom_vline(xintercept = 1000, linetype = "dashed")+
  labs(x = "Response times (ms)")+
  theme_Publication()

sum_RTs$ID[which(sum_RTs$RT > 1000)] # Participants with overall RTs higher than 1000 ms
