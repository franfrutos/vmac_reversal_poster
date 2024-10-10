if (!require(pacman)) {
  install.packages(pacman)
  library(pacman)
}


p_load(dplyr, lme4, lmerTest, marginaleffects, hypr, simr, pbapply, future, parallel, performance)

source("scripts/functions.R")

# Setting model contrasts with hypr package
HcRepRev <- hypr(
  VMAC = High ~ Low,
  AC = High ~ Absent,
  levels = c("High","Low", "Absent")
)

# Model in reversal: 
# Reversal:
d_RT_rev <- filter_data(raw_q, f.absent = F) %>% filter(Phase == "Reversal")
d_RT_rev$log_RT <- log(d_RT_rev$rt)
d_RT_rev$Singleton <- factor(d_RT_rev$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d_RT_rev$Singleton) <- HcRepRev

fRev.l2_log <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) + (Singleton + scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d_RT_rev
  )

summary(fRev.l2_log) # Power function model

ncls <- detectCores()
steps <- seq(20, 100, 20)
nsim <- 1e3
seed <- 123
li <- list(Exp1 = list(), SESOI = list())
for (i in rev(names(li))) {
  print(paste("Effect:",i)) 
  for (j in 1:length(steps)) {
    print(paste("N:", steps[j], "Participants"))
    set.seed(seed)
    model <- extend(fRev.l2_log, along = "ID", n = steps[j]) #extend the model to include n participants
    if (i == "SESOI") {
      fixef(model)[5] <- fixef(model)[5] *.7
      fixef(model)[2] <- fixef(model)[2] *.7
    }
    dat <- getData(model)
    cl <- makeCluster(ncls)
    clusterExport(cl, c("sim_lme4", "doSim", "lmer", "dat", "model", "fRev.l2_log", "check_convergence"))
    power <- pbreplicate(n = nsim, expr = sim_lme4(dat, model), simplify = T, cl = cl)
    stopCluster(cl)
    li[[i]][[paste(j)]]$VMAC <- binom::binom.confint(sum(unlist(power[1,])), nsim, method = "exact")
    li[[i]][[paste(j)]]$Int <- binom::binom.confint(sum(unlist(power[2,])), nsim, method = "exact")
    li[[i]][[paste(j)]]$Cov <- binom::binom.confint(sum(unlist(power[3,])), nsim, method = "exact")
    print(li[[i]][[j]])
  }
}


power <- do.call(rbind,lapply(li[["SESOI"]], FUN = function(x) return(do.call(rbind, x)))) 

power$Effect <- rep(c("VMAC", "VMAC*Block interaction", "Check"), length.out = nrow(power))
power$N <- rep(steps, each = 3)

ggplot(power[power$Effect != "Check",], aes(x = N, y = mean, color = Effect)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = steps, labels = steps) +
  coord_cartesian(ylim = c(0, 1))+
  scale_y_continuous(labels = paste0(seq(0, 1, .1)*100, "%"), breaks = seq(0, 1, .1)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Effect), alpha = 1/3, color = NA)+
  geom_hline(yintercept = .80, linetype = "dashed") +
  labs(y = "Power", x = "Sample size")+
  scale_fill_manual(values=c("darkblue", "darkorange")) +
  scale_color_manual(values=c("darkblue", "darkorange")) +
  theme_Publication()


ggsave(
  "plots/power_plot.png",
  height = 12,
  width = 15,
  dpi = 1200,
  units = "cm"
)
