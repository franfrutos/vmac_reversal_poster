source("scripts/load_data.R") # Automatically load data

if (!require(pacman)) {
  install.packages(pacman)
  library(pacman)
}

p_load(dplyr, Hmisc, grid, ggthemes)

# Function to filter data:
filter_data <- function(data,
                        f.absent = T,
                        acc = T,
                        sd_filter = NULL,
                        fixed = T,
                        two_t = F,
                        phase = "Rewarded",
                        experiment = "e1"
) {
  
  filt_trial <- NA
  if (two_t) {
    filt_trial <- c(seq(1, 577, 24), seq(2, 577, 24))
  }
  
  out <- data %>%
    filter(!(trial_num %in% filt_trial))
  
  if (experiment == "e1") {
    out <- out[which(!out$ID %in% acc_ex,),]
  } else {
    out <- out[which(!out$ID %in% acc_ex2),]
  }

  if (acc)
    out <- out[which(out$correct == 1),]
  
  
  if (f.absent)
    out <- out[which(out$Singleton != "Absent"),]
  
  if (fixed)
    out <- out[which(out$rt > 150 & out$rt < 1800),]
  
  if (is.numeric(sd_filter) & ifelse(is.null(sd_filter), F, sd_filter %in% 1:3)) {
    out <- out %>%
      group_by(ID) %>%
      dplyr::mutate(high_rt = mean(rt, na.rm = T) + sd(rt, na.rm = T)*sd_filter,
                    low_rt = mean(rt, na.rm = T) - sd(rt, na.rm = T)*sd_filter) %>%
      ungroup()  %>%
      filter(rt > low_rt, rt < high_rt) %>%
      select(-c(high_rt, low_rt))
  }
  return(out)
}

# Function to collapse blocks in different epochs.
create_epochs <- function(blocks, epoch = 2) {
  vapply(blocks, function(x, e = epoch) {
    ceiling(x / e)
  }, FUN.VALUE = numeric(1))
}

# Theme used in the plots
theme_Publication <-
  function(base_size = 12,
           base_family = "sans",
           text_size = 11) {
    (
      theme_foundation(base_size = base_size, base_family = base_family)
      + theme(
        plot.title = element_text(
          face = "bold",
          size = rel(1.2),
          hjust = 0.5
        ),
        text = element_text(size = text_size),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.border = element_rect(colour = NA),
        axis.title = element_text(face = "bold", size = rel(1)),
        axis.title.y = element_text(angle = 90, vjust = 2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(),
        panel.grid.major = element_line(colour = "#f0f0f0"),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(colour = NA),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size = unit(0.2, "cm"),
        #legend.margin = margin(0, "cm"),
        #legend.title = element_text(face="italic"),
        plot.margin = unit(c(10, 5, 5, 5), "mm"),
        strip.background = element_rect(colour = NA, fill = NA),
        strip.text = element_text(face = "bold")
      )
    )
    
  }

# Function to transform epochs of raw data to the same scale as the plot
transform_epoch <- function(valor = 2) {
  return(mean(((valor - 1) * valor + 1):(valor * valor)))
}

transform_epoch <- function(epochs, length_v) {
  n <- length(epochs)
  paso <- length_v / n
  ajuste <- sapply(1:n, function(i) {
    mean(seq((i - 1) * paso + 1, i * paso))
  })
  return(ajuste)
}

# Functions to plot

# This function gets predictions for both models presented in the main text to plot (Figure 2)
# d indicate raw data fed in a model for each phase. m indicate models for each phase
# Epoch is set to TRUE, which use epochs of 2 for raw data
# returns a list with 1) averaged raw data by Phase, Block (or Epoch) and Singleton
# and 2) model predictions for both models

get_predictions <- function(d1, d2, m1, m2, epoch = T) {
  library(marginaleffects)
  library(dplyr)
  
  if (epoch) {
    d1$Block_num <- create_epochs(d1$Block)
    d2$Block_num <- create_epochs(d2$Block)
  }
  raw_data <- rbind(
    d1 %>% group_by(ID, Phase, Block_num, Singleton) %>% dplyr::summarise(estimate = mean(rt)) %>%
      Rmisc::summarySEwithin(
        data = .,
        measurevar = "estimate",
        withinvars = c("Phase", "Block_num", "Singleton"),
        idvar = "ID"
      ) %>% mutate(Phase = "Acquisition"),
    d2 %>% group_by(ID, Phase, Block_num, Singleton) %>% dplyr::summarise(estimate = mean(rt)) %>%
      Rmisc::summarySEwithin(
        data = .,
        measurevar = "estimate",
        withinvars = c("Phase", "Block_num", "Singleton"),
        idvar = "ID"
      ) %>% mutate(Phase = "Reversal")
  )
  raw_data$Singleton <-
    factor(raw_data$Singleton, levels = c("High", "Low", "Absent"))
  raw_data$Block_num <- as.numeric(raw_data$Block_num) * 2 - .5
  
  model_preds <- rbind(
    predictions(
      m1,
      newdata = datagrid(
        Singleton = unique,
        Block_num = seq(1, 12, .01),
        ID = NA
      ),
      re.form = NA,
      transform = \(x) exp(x + (sigma(m1) ^ 2) / 2),
    ) %>% mutate(Phase = "Acquisition"),
    predictions(
      m2,
      newdata = datagrid(
        Singleton = unique,
        Block_num = seq(13, 24, .01),
        ID = NA
      ),
      re.form = NA,
      transform = \(x) exp(x + (sigma(m2) ^ 2) / 2)
    ) %>% mutate(Phase = "Reversal")
  )
  
  model_preds$Singleton <-
    factor(model_preds$Singleton, levels = c("High", "Low", "Absent"))
  
  return(list(raw = raw_data,
              mod = model_preds))
}

# Intermediate function to get VMAC and AC effects in get_comparisons()
get_raw_effect <- function(d, epoch = T) {
  if (epoch)
    d$Block_num <- create_epochs(d$Block_num)
  effs <- d %>%
    group_by(ID, Phase, Block_num, Singleton) %>%
    dplyr::summarise(RT = mean(rt)) %>%
    spread(Singleton, RT) %>%
    mutate(VMAC = High-Low,
           AC = Low-Absent) %>%
    ungroup() %>%
    drop_na()
  
  return(
    rbind(
      Rmisc::summarySEwithin(
        data = effs,
        measurevar = "VMAC",
        withinvars = c("Phase", "Block_num"),
        idvar = "ID"
      ) %>%
        dplyr::rename("estimate" = "VMAC") %>% mutate(Block_num = as.numeric(Block_num),
                                                      Effect = "VMAC"),
      Rmisc::summarySEwithin(
        data = effs,
        measurevar = "AC",
        withinvars = c("Phase", "Block_num"),
        idvar = "ID"
      ) %>% dplyr::rename("estimate" = "AC") %>%
        mutate(Block_num = as.numeric(Block_num),
               Effect = "Attentional Capture")
    )
  )
  
}

# This function gets conditional effects for each contrasts in both models presented in the main text to plot (Figure 3)
# d indicate raw data fed in a model for each phase. m indicate models for each phase
# Epoch is set to TRUE, which use epochs of 2 for raw data
# returns a list with 1) averaged raw data by Phase, Block (or Epoch) and Singleton
# and 2) model predictions for both models

get_comparisons <- function(d1, d2, m1, m2, epoch = T) {
  raw_data <- rbind(get_raw_effect(d1, epoch),
                    get_raw_effect(d2, epoch))
  raw_data$Effect <-
    factor(raw_data$Effect, levels = c("VMAC", "Attentional Capture"))
  raw_data$Block_num[raw_data$Phase == "Reversal"] <-
    raw_data$Block_num[raw_data$Phase == "Reversal"] + 6
  raw_data$Block_num <- as.numeric(raw_data$Block) * 2 - .5
  
  
  model_comps <- bind_rows(
    comparisons(
      m1,
      variables = list(Singleton = "revsequential"),
      newdata = datagrid(Block_num = seq(1, 12, .01),
                         ID = NA),
      re.form = NA,
      transform_pre = \(hi, lo) exp(hi + (sigma(m1) ^ 2) / 2) - exp(lo +
                                                                      (sigma(m1) ^ 2) / 2),
    ) %>% mutate(
      Effect = rep(
      c("VMAC", "Attentional Capture"), each = nrow(.) / 2
    ), Phase = "Acquisition"),
    comparisons(
      m2,
      variables = list(Singleton = c("revsequential")),
      newdata = datagrid(Block_num = seq(13, 24, .01),
                         ID = NA),
      re.form = NA,
      transform_pre = \(hi, lo) exp(hi + (sigma(m2) ^ 2) / 2) - exp(lo +
                                                                      (sigma(m2) ^ 2) / 2),
    ) %>% 
      mutate(
         Effect = rep(
           c("VMAC", "Attentional Capture"), each = nrow(.) / 2
         ),
        Phase = "Reversal")
  )
  
   model_comps$Effect <-
     factor(model_comps$Effect, levels = c("VMAC", "Attentional Capture"))
  
  return(list(raw = raw_data,
              mod = model_comps))
}

sim_lme4 <- function(dat, model){
  library(lme4)
  library(lmerTest)
  dat$log_RT <- doSim(model)
  f <- model@call[["formula"]]
  tmp <- lmer(formula = f, data = dat,
              control = lmerControl(optimizer = "bobyqa"))
  return(list(main=summary(tmp)$coefficients[,5][2] < .05,# Main VMAC
              inter = summary(tmp)$coefficients[,5][5] < .05, # VMAC Block interaction
              conv = check_convergence(tmp)))
}

# Report trial exclusions: 
report_exclusions <- function(d, acc) {
  incorrect <- (1-nrow(d[which(d$correct == 1 & !d$ID %in% acc),])/
     nrow(d[which(!d$ID %in% acc),]))*100
  outliers <- (1-nrow(d[which(d$correct == 1 & !d$ID %in% acc & (d$rt > 150 & d$rt < 1800)),])/
      nrow(d[which(d$correct == 1 & !d$ID %in% acc),]))*100 
  
  cat(paste0("Incorrect responses: ", round(incorrect, 2), "%",
             "\nOutliers RTs: ", round(outliers, 2), "%"))
}

