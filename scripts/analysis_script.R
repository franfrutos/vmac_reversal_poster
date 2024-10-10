# Loading packages----
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

# Contrast for Singleton distractor condition (VMAC and AC effects):
HcRep <- hypr(
  VMAC = High ~ Low,
  AC = Low ~ Absent,
  levels = c("High", "Low", "Absent")
)

HcRepRev <- hypr(
  VMAC = High ~ Low, # Negative VMAC effect when there is a reversal
  AC = Low ~ Absent, # New low is compared against absent
  levels = c("High", "Low", "Absent")
)

# RT analysis: ----
# Acquisition model: 
d_RT <- filter_data(raw, f.absent = F) %>% filter(Phase == "Acquisition")
d_RT$log_RT <- log(d_RT$rt)
d_RT$Singleton <- factor(d_RT$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d_RT$Singleton) <- HcRep

f.l <-
  lmer(
    log_RT ~ Singleton * scale(Block_num) + (Singleton * scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d_RT
  )

summary(f.l) 

f.l2 <-
  lmer(
    log_RT ~ Singleton * scale(Block_num) + (Singleton + scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d_RT
  )

summary(f.l2) # Maximal model

f.l2_log <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) + (Singleton + scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d_RT
  )

summary(f.l2_log) # Log-maximal model

anova(f.l2, f.l2_log) ## |AIC_Log - AIC_Linear| = 200.3

# Reversal model:
d_RT_rev <- filter_data(raw, f.absent = F) %>% filter(Phase == "Reversal")
d_RT_rev$log_RT <- log(d_RT_rev$rt)
d_RT_rev$Singleton <- factor(d_RT_rev$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d_RT_rev$Singleton) <- HcRepRev

fRev.l <-
  lmer(
    log_RT ~ Singleton * scale(Block_num) + (Singleton * scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d_RT_rev
  )

summary(fRev.l) 

fRev.l2 <-
  lmer(
    log_RT ~ Singleton * scale(Block_num) + (Singleton + scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d_RT_rev
  )

summary(fRev.l2) 

fRev.l2_log <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) + (Singleton + scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d_RT_rev
  )

summary(fRev.l2_log) 

anova(fRev.l2, fRev.l2_log)

# ACC analysis----
# Acquisition model: 
d_acc <- filter_data(raw_q, f.absent = F, acc = F) %>% filter(Phase == "Rewarded")
d_acc$Singleton <- factor(d_acc$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d_RT$Singleton) <- HcRep

f.acc <-
  glmer(
    correct ~ Singleton * scale(Block_num) + (Singleton * scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    data = d_acc
  )

summary(f.acc) 

f.acc2 <-
  glmer(
    correct ~ Singleton * scale(Block_num) + (Singleton + scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    data = d_acc
  )

summary(f.acc2)

f.acc3 <-
  glmer(
    correct ~ Singleton * scale(Block_num) + (Singleton | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    data = d_acc
  )

summary(f.acc3) # Maximal model


# Reversal model:

d_acc_rev <- filter_data(raw_q, f.absent = F, acc = F) %>% filter(Phase == "Reversal")
d_acc_rev$Singleton <- factor(d_acc_rev$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d_RT$Singleton) <- HcRepRev

fRev.acc <-
  glmer(
    correct ~ Singleton * scale(Block_num) + (Singleton * scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    data = d_acc_rev
  )

summary(fRev.acc) 

fRev.acc2 <-
  glmer(
    correct ~ Singleton * scale(Block_num) + (Singleton + scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    data = d_acc_rev
  )

summary(fRev.acc2)

fRev.acc3 <-
  glmer(
    correct ~ Singleton * scale(Block_num) + (Singleton | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    data = d_acc_rev
  )

summary(fRev.acc3) # Maximal model


# Plots ----

# Model predictions:

# Returns a list with raw averaged data and model preds
preds <- get_predictions(d_RT, d_RT_rev, f.l2_log, fRev.l2_log)

 # Plot predictions

ggplot(data = preds[["mod"]],
       aes(
         y = estimate,
         x = Block_num,
         color = Singleton,
         fill = Singleton
       )) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .3,
              color = NA) +
  facet_wrap( ~ Phase, scales = "free_x") +
  geom_point(data = preds[["raw"]],
             aes(x = as.numeric(Block_num)),
             position = position_dodge(.5)) +
  geom_errorbar(
    data = preds[["raw"]],
    aes(
      x = as.numeric(Block_num),
      ymin = estimate - se,
      ymax = estimate + se
    ),
    position = position_dodge(.5),
    width = .01
  ) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  scale_y_continuous(breaks = seq(600, 900, 50)) +
  scale_color_brewer(palette = "Set1", name = "Singleton type") +
  scale_fill_brewer(palette = "Set1", name = "Singleton type") +
  labs(y = "Response time (ms)", x = "Block") +
  theme_Publication(text_size = 10) + 
  theme(legend.spacing.x = unit(.2, 'cm'))

ggsave(
  "plots/predictions_RT.png",
  height = 12,
  width = 15,
  dpi = 1200,
  units = "cm"
)


# Model comparisons:
# Returns a list with raw averaged VMAC and ACC effects and the conditional effect of
# Singleton as a function of block

comps <- get_comparisons(d_RT, d_RT_rev, f.l2_log, fRev.l2_log)

ggplot(data = comps[["mod"]],
       aes(
         y = estimate,
         x = Block_num,
         color = Effect,
         fill = Effect
       )) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .3,
              color = NA) +
  facet_wrap( ~ Phase, scales = "free_x") +
  geom_point(
    data = comps[["raw"]],
    aes(x = as.numeric(Block_num)),
             position = position_dodge(.5)) +
  geom_errorbar(
    data = comps[["raw"]],
    aes(
      x = as.numeric(Block_num),
      ymin = estimate - se,
      ymax = estimate + se
    ),
    position = position_dodge(.5),
    width = .01
  ) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  #scale_y_continuous(breaks = seq(-20, 60, 20)) +
  scale_color_brewer(palette = "Set1", name = "Contrast") +
  scale_fill_brewer(palette = "Set1", name = "Contrast") +
  labs(y = "Contrast (ms)", x = "Block") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_Publication(text_size = 10) + # text size is adjusted for DPI
  theme(legend.spacing.x = unit(.2, 'cm'))

ggsave(
  "plots/conditional_effs_RT.png",
  height = 12,
  width = 15,
  dpi = 1200,
  units = "cm"
) 
