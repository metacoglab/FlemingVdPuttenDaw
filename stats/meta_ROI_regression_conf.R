# Regression models for ROI data in Fleming, van der Putten & Daw
# Focuses on single-trial betas aligned at confidence onset
# For more details of setting lme4 options see https://cran.r-project.org/web/packages/lme4/vignettes/lmerperf.html
#
# Steve Fleming 2016

rm(list=ls())
require(R.matlab) 
require(lme4)
require(car)
require(doBy)
require(optimx)
require(stargazer)
options(contrasts = c("contr.treatment", "contr.poly")) # This is R defaults but set it anyway to be safe
source("~/Dropbox/Research/Metacognition/stateactionexpt/github/stats/errorBar.R")

dataDir = "~/Dropbox/Research/Metacognition/stateactionexpt/github/fmri/roidata_conf"
subjects =c(seq(12,19), seq(23,28), seq(30,37))
rois = c("pMFC", "union_46", "union_FPL", "union_FPm", "vmPFC", "rVS_FSL_structAtlas")

# Extract out data from mat files
split_diff_pval = NA
for (r in 1:length(rois)) {
  bigData = NULL
  
  for (s in 1:length(subjects)){
    
    ## EXTRACT DATA FILE
    setwd(dataDir)
    DATA = readMat(paste('ROI_', rois[r], '_sub', subjects[s], '.mat',sep=""))
    
    dat = DATA$ROI
    
    precoh = NULL
    postcoh = NULL
    logRT = NULL
    conf = NULL
    logConfRT = NULL
    accuracy = NULL
    BOLD = NULL
    
    BOLD = scale(dat[,,1]$Y[,1])
    precoh = dat[,,1]$precoh[,1] - 2
    postcoh = dat[,,1]$postcoh[,1] - 2
    logRT = scale(dat[,,1]$logRT[,1])
    conf = dat[,,1]$conf[,1]
    QSR = pmax(1-(1 - conf)^2, 1-(0 - conf)^2)
    logConfRT = scale(dat[,,1]$logConfRT[,1])
    accuracy = dat[,,1]$acc[,1]
    
    QSR = scale(QSR)
    rawConf = conf
    conf = scale(conf)
    subj = rep(s, length(accuracy))
    subData = data.frame("subj"=subj, "BOLD"=BOLD, "precoh"=precoh, "postcoh"=postcoh, "logRT"=logRT, "logConfRT"=logConfRT, "conf"=conf, "accuracy"=accuracy, "QSR"=QSR, "rawConf"=rawConf)
    
    bigData = rbind(bigData, subData)
    
  }
  
  bigData$subj <- factor(bigData$subj)
  bigData$accuracy <- factor(bigData$accuracy)

  # Fit regression models

  # Split confidence above and below 0.5
  bigData$rawConf_low = rep(0, length(bigData$rawConf))
  bigData$rawConf_low[bigData$rawConf <= 0.5] = bigData$rawConf[bigData$rawConf <= 0.5]
  bigData$rawConf_high = rep(0, length(bigData$rawConf))
  bigData$rawConf_high[bigData$rawConf > 0.5] = bigData$rawConf[bigData$rawConf > 0.5]
  confModel_split = lmer(BOLD ~ rawConf_low + rawConf_high + logRT + (1 + rawConf_low + rawConf_high + logRT | subj), data=bigData, control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
  fix <- fixef(confModel_split)
  fix.se <- sqrt(diag(vcov(confModel_split)))
  barX <- barplot(fix, beside = TRUE, ylab = "Regression betas", names.arg = names(fix), cex.names = 0.5)
  errorBar(barX, fix, fix.se)
  betas <- c(fix, fix.se)
  # write.csv(betas, file = paste('regression_conf_split_', rois[r], '.csv', sep=""))
  
  print(paste("ROI analysis for ", rois[r], sep=""))
  print(summary(confModel_split))
  print(Anova(confModel_split, type=3))
  
  # Contrasts between high and low
  lambda1 <- c(0, 1, -1, 0)
  out = print(esticon(confModel_split,lambda1))
  split_diff_pval[r] = out$`Pr(>|X^2|)`
  
  # Store as individual model for plotting in table
  eval(parse(text = paste("confModel_split", r, "= confModel_split", sep="")))
  
}

stargazer(confModel_split1, confModel_split2, confModel_split3, confModel_split4, confModel_split5, confModel_split6, star.cutoffs = c(0.05, 0.01, 0.001), keep.stat = c("n"), dep.var.labels = c("BOLD"), column.labels = c("pMFC", "area 46", "FPl", "FPm", "vmPFC", "v. striatum"), covariate.labels = c("Confidence <= 0.5", "Confidence > 0.5", "log(RT)"))

