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
       effectsize, # Cohen's d
       ggplot2 # Plotting
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

# RT ----
# Get data for analysis: 
d2 <- filter_data(raw2, f.absent = F, experiment = "e2")
d2$Epoch <- create_epochs(d2$Block_num, 2)
d2$log_RT <- log(d2$rt)
d2$Singleton <- factor(d2$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d2$Singleton) <- HcRep

# Experiment 2A
f_r <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) * (Singleton+scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Acquisition" & d2$Experiment == "A", ]
  )

summary(f_r)

contrasts(d2$Singleton) <- HcRepRev

f_rev <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) * (Singleton+scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Reversal" & d2$Experiment == "A", ]
  )

summary(f_rev)

contrasts(d2$Singleton) <- HcRep

# Experiment 2B
f_rb <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) + (Singleton+scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Acquisition" & d2$Experiment == "B" & d2$Block_num > 6, ]
  )

f_rb <-
  lmer(
    log_RT ~ Singleton + (Singleton | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Acquisition" & d2$Experiment == "B" & d2$Block_num > 6, ]
  )

summary(f_rb)

contrasts(d2$Singleton) <- HcRepRev

f_revb <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) + (Singleton+scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Reversal" & d2$Experiment == "B", ]
  )

summary(f_revb)

# Accuracy (do properly)----
d_acc <- filter_data(raw2, f.absent = F, acc = F, experiment = "e2") %>% filter(Phase != "Awareness")
d_acc$Singleton <- factor(d_acc$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d_acc$Singleton) <- HcRep

# Trick to improve performance (evaluating likelihood as binomial(p, n), which gave same results)
d_opt <- d_acc %>%
  dplyr::summarise(correct = mean(correct),
                   n = n(),
                   .by = c(Experiment, ID, Phase, Block_num, Singleton))
  
f_r_acc <-
  glmer(
    correct ~ Singleton * scale(log(Block_num)) + (Singleton * scale(log(Block_num)) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n,
    data = d_opt[which(d_opt$Phase == "Acquisition" & d_opt$Experiment == "A"), ]
  ) # Maximal model

summary(f_r_acc) # Singular fit

# Dropping interaction from RE
f_r_acc <-
  glmer(
    correct ~ Singleton * scale(log(Block_num)) + (Singleton + scale(log(Block_num)) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n,
    data = d_opt[which(d_opt$Phase == "Acquisition" & d_opt$Experiment == "A"), ]
  ) 

summary(f_r_acc) # Singular fit

# Dropping Block from RE
f_r_acc <-
  glmer(
    correct ~ Singleton * scale(log(Block_num)) + (Singleton | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n,
    data = d_opt[which(d_opt$Phase == "Acquisition" & d_opt$Experiment == "A"), ]
  ) 

summary(f_r_acc) # Maximal feasible model

contrasts(d_acc$Singleton) <- HcRepRev

f_rev_acc <-
  glmer(
    correct ~ Singleton * scale(log(Block_num)) + (1 | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    data = d_acc[which(d_acc$Phase == "Reversal" & d_acc$Experiment == "A"), ]
  )

summary(f_rev_acc)

contrasts(d_acc$Singleton) <- HcRep

f_r_accb <-
  glmer(
    correct ~ Singleton * scale(log(Block_num)) + (scale(log(Block_num)) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    data = d_acc[which(d_acc$Phase == "Acquisition" & d_acc$Experiment == "B"), ]
  )

summary(f_r_accb)

contrasts(d_acc$Singleton) <- HcRepRev

f_rev_accb <-
  glmer(
    correct ~ Singleton * scale(log(Block_num)) + (1 | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    data = d_acc[which(d_acc$Phase == "Reversal" & d_acc$Experiment == "B"), ]
  )

summary(f_rev_accb)

plot_predictions(
  f_rev_acc,
  newdata = insight::get_data(f_r_acc),
  condition = c("Singleton"),
  re.form = NA
)

plot_comparisons(
  f_r_accb,
  variables = list(Singleton = "revsequential"),
  condition = c("Block_num"),
  re.form = NA
)


# Is there persistence at the beginning of the reversal phase compared to the end of acquisition? ----

d2[d2$Block_num %in% 11:14, ] %>%
  mutate(Experiment = ifelse(Experiment == "A", "Experiment 2A", "Experiment 2B"),
         Phase = ifelse(Phase == "Rewarded", "Acquisition", Phase),
         Singleton = case_when(
           T~Singleton
         )) %>%
  dplyr::summarise(RT = mean(rt), .by = c(University, Experiment, ID, Phase, Singleton)) %>%
  mutate(RT =RT)-> d_aov

d_aov$Experiment <- factor(d_aov$Experiment, levels = c("Experiment 2A", "Experiment 2B"))

afex::aov_ez(data = d_aov %>% filter(Singleton != "Absent", Experiment == "Experiment 2A"), id = "ID", dv = "RT",
             within = c("Phase", "Singleton"), anova_table = list(es = "pes")) # Significant interaction

afex::aov_ez(data = d_aov %>% filter(Singleton != "Absent", Experiment == "Experiment 2B"), id = "ID", dv = "RT",
             within = c("Phase", "Singleton"), anova_table = list(es = "pes")) # Non-significant interaction


# Comparison of experiments plot
d_aov %>%
  mutate(ID = as.factor(paste(ID, Experiment))) %>%
  filter(Singleton != "Absent") %>%
  bind_rows(.,
    filter_data(raw, f.absent = F) %>% mutate(Experiment = "Experiment 1") %>%
      filter(Block_num %in% 11:14) %>%
      mutate(
             Phase = ifelse(Phase == "Rewarded", "Acquisition", Phase),
             Singleton = case_when(
               Phase == "Reversal" & Singleton == "Low"~"High",
               Phase == "Reversal" & Singleton == "High"~"Low",
               T~Singleton
             )) %>%
      dplyr::summarise(RT = mean(rt), .by = c(Experiment, ID, Phase, Singleton)) %>%
  filter(Singleton != "Absent") %>%
  mutate(
    ID = as.factor(paste(ID)),
    Singleton = as.factor(Singleton),
    Phase = as.factor(Phase),
    #Experiment = factor(Experiment, levels = c("Experiment 1", "Experiment 2A","Experiment 2B")),
    RT = RT
  )) -> merge_data

merge_data$Experiment <- factor(merge_data$Experiment, levels = c("Experiment 1", "Experiment 2A","Experiment 2B"))
merge_data$Singleton <- factor(merge_data$Singleton, levels = c("High", "Low"))


afex::aov_ez(data = merge_data, id = "ID", dv = "RT",
               within = c("Phase", "Singleton"), between = "Experiment", anova_table = list(es = "pes"))%>%
  afex_plot(., x="Phase", trace = "Singleton", panel = "Experiment", error = "CMO", mapping = c("color", "shape"),
            dodge = .1,
            data_alpha = 0.25,
            data_arg = list(
              position = 
                ggplot2::position_jitterdodge(
                  jitter.width = .2, 
                  jitter.height = 0, 
                  dodge.width = 0.1  ## needs to be same as dodge
                ))) +
  theme_Publication() +
  labs(y = "Reponse times (ms)") +
  scale_color_brewer(palette = "Set1")

ggsave(
  "plots/reversal_exps.png",
  height = 15,
  width = 20,
  dpi = 1200,
  units = "cm"
) 

# Plots ----
# Experiment 2A
# Predictions:
preds <- get_predictions(d2[d2$Phase == "Acquisition" & d2$Experiment == "A", ],
                         d2[d2$Phase == "Reversal" & d2$Experiment == "A", ],
                         f_r, f_rev)
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
  "plots/preds2A_RT.png",
  height = 15,
  width = 20,
  dpi = 1200,
  units = "cm"
) 

# Conditional effects: 
comps <- get_comparisons(d2[d2$Phase == "Acquisition" & d2$Experiment == "A", ],
                         d2[d2$Phase == "Reversal" & d2$Experiment == "A", ],
                         f_r, f_rev)

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
    data = comps[["raw"]] %>% mutate(Phase = ifelse(Phase == "Acquisition", "Acquisition",
                                                    "Reversal")),
    aes(x = as.numeric(Block_num)),
    position = position_dodge(.5)) +
  geom_errorbar(
    data = comps[["raw"]] %>% mutate(Phase = ifelse(Phase == "Acquisition", "Acquisition",
                                                    "Reversal")),
    aes(
      x = as.numeric(Block_num),
      ymin = estimate - se,
      ymax = estimate + se
    ),
    position = position_dodge(.5),
    width = .01
  ) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  scale_y_continuous(breaks = seq(-20, 60, 20)) +
  scale_color_brewer(palette = "Set1", name = "Contrast") +
  scale_fill_brewer(palette = "Set1", name = "Contrast") +
  labs(y = "Contrast (ms)", x = "Block") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_Publication(text_size = 10) + # text size is adjusted for DPI
  theme(legend.spacing.x = unit(.2, 'cm'))

ggsave(
  "plots/cond2A_RT.png",
  height = 12,
  width = 15,
  dpi = 1200,
  units = "cm"
) 

# Experiment 2B
# Predictions
predsB <- get_predictions(d2[d2$Phase == "Acquisition" & d2$Experiment == "B", ],
                         d2[d2$Phase == "Reversal" & d2$Experiment == "B", ],
                         f_rb, f_revb)

ggplot(data = predsB[["mod"]],
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
  geom_point(data = predsB[["raw"]],
             aes(x = as.numeric(Block_num)),
             position = position_dodge(.5)) +
  geom_errorbar(
    data = predsB[["raw"]],
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
  "plots/preds2B_RT.png",
  height = 15,
  width = 20,
  dpi = 1200,
  units = "cm"
) 

comps <- get_comparisons(d2[d2$Phase == "Acquisition" & d2$Experiment == "B", ],
                         d2[d2$Phase == "Reversal" & d2$Experiment == "B", ],
                         f_rb, f_revb)

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
    data = comps[["raw"]] %>% mutate(Phase = ifelse(Phase == "Acquisition", "Acquisition",
                                                    "Reversal")),
    aes(x = as.numeric(Block_num)),
    position = position_dodge(.5)) +
  geom_errorbar(
    data = comps[["raw"]] %>% mutate(Phase = ifelse(Phase == "Acquisition", "Acquisition",
                                                    "Reversal")),
    aes(
      x = as.numeric(Block_num),
      ymin = estimate - se,
      ymax = estimate + se
    ),
    position = position_dodge(.5),
    width = .01
  ) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  scale_y_continuous(breaks = seq(-20, 60, 20)) +
  scale_color_brewer(palette = "Set1", name = "Contrast") +
  scale_fill_brewer(palette = "Set1", name = "Contrast") +
  labs(y = "Contrast (ms)", x = "Block") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_Publication(text_size = 10) + # text size is adjusted for DPI
  theme(legend.spacing.x = unit(.2, 'cm'))

ggsave(
  "plots/cond2B_RT.png",
  height = 12,
  width = 15,
  dpi = 1200,
  units = "cm"
) 

# Awareness ----
d_aw <- raw2[raw2$Phase == "Awareness" & raw2$Singleton != "Absent" & !raw2$ID %in% acc_ex2,]
d_aw$response_aw <- ifelse(d_aw$response == "c", 1, 0)
d_aw$signal <- ifelse(d_aw$Singleton == "High", .5, -.5)
d_aw$Experiment <- factor(d_aw$Experiment, levels = unique(d_aw$Experiment))
contrasts(d_aw$Experiment) <- contr.sum(2)/2

fit_aw <- glmer(response_aw~signal*Experiment+(signal|ID),data=d_aw,family=binomial(link="probit"),
                control = glmerControl(optimizer = "bobyqa"))


get_roc <- function(d, c) {
  fas <- seq(0, 1, 0.001)
  hts <- pnorm(d + qnorm(fas))
  return(list(dp = tibble(hts, fas), c = tibble(fas=pnorm(c),
                                                hts=pnorm(d+c))))
}

get_roc_probit <- function(model) {
  c_s <- predictions(model, datagrid(ID = unique, Experiment = unique, signal = unique),
                     type = "link", re.form = NULL, by = c("Experiment", "signal"))[c(1, 3),]
  
  d_s <- predictions(model, datagrid(ID = unique, Experiment = unique, signal = unique),
                     type = "link", re.form = NULL, by = c("Experiment", "signal"),
                     hypothesis = "sequential")[c(1, 3),]
  
  roc0 <- get_roc(d_s$estimate[1], c_s$estimate[1])
  roc1 <- get_roc(d_s$estimate[2], c_s$estimate[2])
  
  # dprime
  roc0$dp$upp <- get_roc(d_s$conf.high[1], 0)$dp$hts
  roc0$dp$low <- get_roc(d_s$conf.low[1], 0)$dp$hts
  
  roc1$dp$upp <- get_roc(d_s$conf.high[2], 0)$dp$hts
  roc1$dp$low <- get_roc(d_s$conf.low[2], 0)$dp$hts
  
  rocs_d <- rbind(roc0$dp, roc1$dp)
  rocs_d$Experiment <- as.factor(rep(c("Implicit reversal", "Explicit reversal"), each = nrow(rocs_d)/2))
  
  # c 
  
  print(pnorm(-c_s$conf.low))
  print(pnorm(-c_s$conf.high))
  #roc0$c$upp <- get_roc(d_s$estimate[1], c_s$conf.high[1])$c$
  #roc0$c$low <- get_roc(d_s$estimate[1], c_s$conf.low[1])
  
  #roc1$c$upp <- get_roc(d_s$estimate[2], c_s$conf.high[2])
  #roc1$c$low <- get_roc(d_s$estimate[2], c_s$conf.low[2])
  
  rocs_c <- rbind(roc0$c, roc1$c)
  rocs_c$Experiment <- as.factor(c("Implicit reversal", "Explicit reversal"))
  
  # Raw data:
  raw_aw <- d_aw %>% 
    filter(Singleton != "Absent") %>%
    dplyr::summarise(Prop = mean(response_aw), .by = c(ID, Experiment, signal)) %>% drop_na() %>%
    spread(signal, Prop) %>%
    mutate(hts = `0.5`,
           fas = (`-0.5`))
  
  raw_aw$Experiment <- ifelse(raw_aw$Experiment== "A", "Explicit reversal", "Implicit reversal")
  
  plot <- ggplot(rocs_d, aes(1-fas, 1-hts, color = Experiment)) +
    geom_point(data = raw_aw, alpha = .3) +
    geom_line()+
    geom_ribbon(aes(ymin = 1-upp, ymax = 1-low, fill = Experiment), color = NA,
                alpha = .3)+
    geom_point(data = rocs_c, aes(1-fas, 1-hts),
               size = 2) +
    geom_line(data = rocs_d %>% filter(((fas < pnorm(c_s$conf.high[1]) & fas > pnorm(c_s$conf.low[1])) & Experiment == "Implicit reversal") |
                                         ((fas < pnorm(c_s$conf.high[2]) & fas > pnorm(c_s$conf.low[2])) & Experiment == "Explicit reversal")),
              linewidth = 3, alpha = .3, lineend = "round") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    scale_color_manual(values = c("darkred", "darkblue")) +
    scale_fill_manual(values = c("darkred", "darkblue")) +
    coord_cartesian(xlim = c(0.045,.958)) +
    theme_Publication() +
    # labs(x = "False alarm rate", y = "Hit rate", title = "ROC curve for hierarchical EVSDT probit",
    #      subtitle = "Shaded segments represents 95% Confidence Intervals") +
    labs(x = "False alarm rate", y = "Hit rate") +
    theme(aspect.ratio = 1)
  
  return(plot)
}


get_roc_probit(fit_aw) 

ggsave("plots/plots_roc.png", height = 15, width = 15, dpi = 900,
       units = "cm")

