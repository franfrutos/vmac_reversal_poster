# Installing (if needed) and loading the packages ----
if (!require(pacman)) {
  install.packages(pacman)
  library(pacman)
}

p_load(here, dplyr, tidyr, afex)

setwd(here::here())

# Experiment 1----
raw <- read.csv("Input/data_e1.csv") %>%
  mutate(Phase = ifelse(Phase == "Reward", "Acquisition", Phase))

# Check 
# There are participants with less observations than the max number of observations?
check <- raw %>%
  group_by(ID, Phase) %>%
  dplyr::count() %>%
  ungroup() %>%
  complete(ID, Phase) %>% filter(Phase %in% c("Acquisition", "Reversal")) %>%
  filter(n < 288 | is.na(n))

length(unique(check$ID)) # N = 0 participants with missing data

length(unique(raw$ID)) # N = 92 participants

# Participants with less than .7 of accuracy?

acc_ex <- raw %>%
  filter(Phase %in% c("Acquisition", "Reversal")) %>%
  group_by(ID) %>%
  dplyr::summarise(mean_ACC = mean(correct)) %>% filter(mean_ACC < .7) %>% pull(ID)

length(unique(acc_ex)) # 1 participants with less than .7 of accuracy that will be excluded from the analysis

# How many participants?
raw %>%
  filter(!ID %in% acc_ex, Phase %in% c("Acquisition", "Reversal")) %>%
  group_by(ID, Phase) %>%
  dplyr::count() %>%
  group_by(Phase) %>%
  dplyr::count() # 91 participants with complete observations and ACC > .7

# Report trial exclusions:
report_exclusions(raw, acc_ex)

# Sample characteristic:
cat(paste0("Mean age: ", round(mean(quests$Edad), 2), ", SD age: ",
       round(sd(quests$Edad), 2), "\nNumber of self-identified as female: ",
       length(quests$Sexo[quests$Sexo == 2])))

# Experiments 2a and 2b ----
raw2 <- read.csv("Input/data_e2.csv") %>%
  mutate(Phase = ifelse(Phase == "Reward", "Acquisition", Phase),
         Singleton = case_when(
           Singleton == "High" & Phase == "Reversal"~"Low",
           Singleton == "Low" & Phase == "Reversal"~"High",
           T~Singleton
         ))

check <- raw2 %>%
  group_by(ID, Phase) %>%
  dplyr::count() %>%
  ungroup() %>%
  complete(ID, Phase) %>% filter(Phase %in% c("Acquisition", "Reversal")) %>%
  filter(n < 288 | is.na(n))

length(unique(check$ID)) # N = 0 participants with missing data

# How many participants?
acc_ex2 <- raw2 %>%
  filter(Phase %in% c("Acquisition", "Reversal")) %>%
  group_by(ID) %>%
  dplyr::summarise(mean_ACC = mean(correct)) %>%
  filter(mean_ACC < .7) %>%
  pull(ID) # All participants above cut-off

# Excluded participants with mean RTs above 1000 ms (see non_preregistered_exclusions.R)
acc_ex2 <- c(553218, 871382) 

raw2 %>%
  filter(!ID %in% acc_ex2, Phase %in% c("Acquisition", "Reversal")) %>%
  group_by(Experiment, ID, Phase) %>%
  dplyr::count() %>%
  group_by(Experiment, Phase) %>%
  dplyr::count() # 61 and 62 participants with complete observations and ACC > .7

# Report trial exclusions for experiment 2A and 2B:
print("Experiment 2A:")
report_exclusions(raw2[which(raw2$Experiment == "A"),], acc_ex2)

print("Experiment 2B:")
report_exclusions(raw2[which(raw2$Experiment == "B"),], acc_ex2)

# TODO:
# Sample characteristic for experiment 2A:
cat(paste0("Mean age: ", round(mean(quests$Edad), 2), ", SD age: ",
           round(sd(quests$Edad), 2), "\nNumber of self-identified as female: ",
           length(quests$Sexo[quests$Sexo == 2])))

# Sample characteristic for experiment 2B:
cat(paste0("Mean age: ", round(mean(quests$Edad), 2), ", SD age: ",
           round(sd(quests$Edad), 2), "\nNumber of self-identified as female: ",
           length(quests$Sexo[quests$Sexo == 2])))

# Merged dataset ----
raw_all <- bind_rows(
  raw %>%
    mutate(Experiment = "Experiment 1") %>% 
    select(ID, Block_num, trial_num, Phase, Singleton, rt, correct, Experiment),
  raw2 %>% mutate(Experiment = ifelse(Experiment == "A", "Experiment 2A", "Experiment 2B")) %>%
    select(ID, Block_num, trial_num, Phase, Singleton, rt, correct, Experiment)
)


  
