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
       ggplot2, # Plotting
       papaja
)

eyeD <- read.csv("Input/df_le2.csv") %>% 
  filter(fixationTimeout != 1, trialPropGoodSamples > .25) %>%
  mutate(
    timeOnD = case_when(
      distractLoc == 1 ~ timeOnLoc_1,
      distractLoc == 2 ~ timeOnLoc_2,
      distractLoc == 3 ~ timeOnLoc_3,
      distractLoc == 4 ~ timeOnLoc_4,
      distractLoc == 5 ~ timeOnLoc_5,
      TRUE ~ timeOnLoc_6  # Default value if none of the conditions match
    )
  )

group_df <- read.csv("Input/summary_Le.csv") %>% mutate(ID = subNum)

eyeD <- left_join(eyeD, group_df[, c("ID", "revGroup")], by = "ID")

eyeD$distance <- eyeD$targetLoc - eyeD$distractLoc
eyeD$distance <- ifelse(abs(eyeD$distance) > 3, 6 - abs(eyeD$distance), abs(eyeD$distance))
raw_distance <- eyeD %>% filter(phase == 2) %>%
  dplyr::summarise(P = mean(omissionTrial), .by = c("ID", "distance", "Singleton"))

eyeD %>% filter(phase == 2) %>%
  dplyr::summarise(P = mean(omissionTrial), .by = c("ID", "distance", "Singleton")) %>%
  Rmisc::summarySEwithin(idvar = "ID", measurevar = "P", withinvars = c("distance", "Singleton")) %>%
  ggplot(aes(x = distance, color = Singleton, y = P)) +
  geom_point() +
  geom_line(aes(group = Singleton)) +
  geom_errorbar(aes(ymin = P - ci, ymax = P + ci), width = .06) +
  geom_point(data = raw_distance, alpha = .15, position = position_jitter(.1)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Target-distractor distance", y = "Probabilitie of fixating the distractor", color = "Distractor Type") +
  theme_Publication()

# Contrast for VMAC
eyeD$distractType <- ifelse(eyeD$distractType == 1, "High", "Low")
eyeD$Singleton <- factor(eyeD$distractType, levels = c("High", "Low"))
contrasts(eyeD$Singleton) <- contr.sum(2)/2
colnames(contrasts(eyeD$Singleton)) <- "VMAC"

# Contrast for Groups
eyeD<- eyeD %>% mutate(Group = case_when(
  revGroup == 1 ~ "Rev",
  revGroup == 2 ~ "Control",
  T~"RevFB"
))

eyeD$Group <- factor(eyeD$Group, levels = c("Rev", "Control", "RevFB"))

Gcontrast <- hypr(
  Control_Rev = Control ~ Rev,
  Rev_RevFB = Rev ~ RevFB,
  levels = c("Rev", "Control", "RevFB")
)

contrasts(eyeD$Group) <- Gcontrast


# Training:
fit_Le <- glmer(
  omissionTrial~Singleton*scale(block)*Group + (Singleton*scale(block)|ID),
                data = eyeD[eyeD$phase == 2,], control = glmerControl(optimizer = "bobyqa"),
                family=binomial()
  )

summary(fit_Le)


# Test:
fit_LeTest <- glmer(
  omissionTrial~Singleton*scale(block)*Group + (Singleton*scale(block)|ID),
                data = eyeD[eyeD$phase == 3,], glmerControl(optimizer = "bobyqa"),
                family=binomial()
  )

summary(fit_LeTest)

# Model predictions for plot: 
TrainPreds <- comparisons(
  fit_Le,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(block = unique,
                     Group = unique,
                     ID = unique),
  re.form = NULL, by = c("Group", "block")) %>% mutate(Phase = "Training")


TestPreds <- comparisons(
  fit_LeTest,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(block = unique,
                     Group = unique,
                     ID = unique),
  re.form = NULL, by = c("Group", "block")) %>% mutate(Phase = "Test")

TTpredPlots <- rbind(TrainPreds, TestPreds)

TTpredPlots$Phase <- factor(TTpredPlots$Phase, levels = c("Training", "Test"))
TTpredPlots$Group <- ifelse(TTpredPlots$Group == "Rev", "Revaluation only",
                            ifelse(TTpredPlots$Group == "Control", "No revaluation",
                                   "Revaluation + Feedback"))
TTpredPlots$Group <- factor(TTpredPlots$Group, levels = c("No revaluation", "Revaluation only",
                                                          "Revaluation + Feedback"))

Rmisc::summarySEwithin(
  data = effs,
  measurevar = "VMAC",
  withinvars = c("Phase"),
  betweenvars = "Group",
  idvar = "ID"
)
# Raw data for plot: 

eyeD$Block <- create_epochs(eyeD$block)
effs <- eyeD %>% filter(phase != 1) %>%
  mutate(Phase = ifelse(phase == 2, "Training", "Test")) %>%
  group_by(Group, ID, Phase, Block, Singleton) %>%
  dplyr::summarise(Prop = mean(omissionTrial)) %>%
  spread(Singleton, Prop) %>%
  mutate(VMAC = High - Low)

raw_sum <- Rmisc::summarySEwithin(
      data = effs,
      measurevar = "VMAC",
      withinvars = c("Phase", "Block"),
      betweenvars = "Group",
      idvar = "ID"
    )

# Adjusting test:
raw_sum$Phase <- factor(raw_sum$Phase, levels = c("Training", "Test"))
raw_sum$Group <- ifelse(raw_sum$Group == "Rev", "Revaluation only",
                            ifelse(raw_sum$Group == "Control", "No revaluation",
                                   "Revaluation + Feedback"))
raw_sum$Group <- factor(raw_sum$Group, levels = c("No revaluation", "Revaluation only",
                                                          "Revaluation + Feedback"))
raw_sum$estimate <- raw_sum$VMAC
raw_sum$block <- raw_sum$Block

# Adjusting epoch:
raw_sum$Block <- as.numeric(raw_sum$Block) * 2 - .5


# Plot with data: 
ggplot(data = TTpredPlots,
       aes(
         y = estimate*100,
         x = block,
         color = Group,
         fill = Group
         )) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low*100, ymax = conf.high*100),
              alpha = .3,
              color = NA) +
  facet_grid( ~ Phase, scales = "free_x", space = "free_x") +
  geom_point(
    data = raw_sum,
    aes(x = as.numeric(Block)),
    position = position_dodge(.6)) +
  geom_errorbar(
    data = raw_sum,
    aes(
      x = as.numeric(Block),
      ymin = estimate*100 - se*100,
      ymax = estimate*100 + se*100
    ),
    position = position_dodge(.6),
    width = .01
  ) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  #scale_y_continuous(breaks = seq(-20, 60, 20)) +
  scale_color_brewer(palette = "Set1", name = "Group") +
  scale_fill_brewer(palette = "Set1", name = "Group") +
  labs(y = "Distractor difference (%)", x = "Block") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_Publication(text_size = 10) + # text size is adjusted for DPI
  theme(legend.spacing.x = unit(.2, 'cm'))

ggsave(
  "plots/conditional_effs_Le.png",
  height = 12,
  width = 15,
  dpi = 1200,
  units = "cm"
) 

# Comparing phases:
eyeDc <- eyeD %>% filter((block > 14 & phase == 2) | (block < 3 & phase == 3))
eyeDc$Phase <- factor(eyeDc$phase, levels = 2:3)
contrasts(eyeDc$Phase) <- rev(contr.sum(2)/2)
colnames(contrasts(eyeDc$Phase)) <- "Phase"


fit_LeComp <- glmer(
  omissionTrial~Singleton*Phase + (Singleton|ID),
  data = eyeDc[eyeDc$Group == "Rev",], glmerControl(optimizer = "bobyqa"),
  family=binomial()
)

summary(fit_LeComp)

comparisons(
  fit_LeComp,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(Phase = unique,
                     ID = NA),
  re.form = NA)

# Re-analysis on fixation time
