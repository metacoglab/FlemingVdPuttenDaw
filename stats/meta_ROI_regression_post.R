# Regression models for ROI data in Fleming, van der Putten & Daw
# Focuses on single-trial betas aligned at post-decision onset
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
source("~/Dropbox/Published/FlemingNN2018/Code&analysis/FlemingVdPuttenDaw/stats/errorBar.R")

dataDir = "~/Dropbox/Published/FlemingNN2018/Code&analysis/FlemingVdPuttenDaw/fmri/roidata_post"
modelDir = "~/Dropbox/Published/FlemingNN2018/Code&analysis/FlemingVdPuttenDaw/regressors"
subjects =c(seq(12,19), seq(23,28), seq(30,37))
rois = c("pMFC", "union_46", "union_FPL", "union_FPm", "vmPFC", "rVS_FSL_structAtlas")

# Extract out data from mat files
# This is a bit inefficient as only the ROI signal changes from loop to loop
prepost_pval = NA
for (r in 1:length(rois)) {
  bigData = NULL
  
  for (s in 1:length(subjects)){
    
    ## EXTRACT DATA FILES
    setwd(dataDir)
    DATA = readMat(paste('ROI_', rois[r], '_sub', subjects[s], '.mat',sep=""))
    dat = DATA$ROI
    
    setwd(modelDir)
    modelDATA = readMat(paste('ideal_dt', subjects[s], '_loglikPost.mat',sep=""))
    loglikpre_cor = modelDATA$model.loglikPre.cor
    loglikpost_cor = modelDATA$model.loglikPost.cor
    loglikpre_err = modelDATA$model.loglikPre.err
    loglikpost_err = modelDATA$model.loglikPost.err
    
    precoh = NULL
    postcoh = NULL
    logRT = NULL
    conf = NULL
    logConfRT = NULL
    accuracy = NULL
    BOLD = NULL
    
    BOLD = scale(dat[,,1]$Y[,1])
    precoh = dat[,,1]$precoh[,1]
    postcoh = dat[,,1]$postcoh[,1]
    logRT = scale(dat[,,1]$logRT[,1])
    conf = dat[,,1]$conf[,1]
    QSR = pmax(1-(1 - conf)^2, 1-(0 - conf)^2)
    logConfRT = scale(dat[,,1]$logConfRT[,1])
    accuracy = dat[,,1]$acc[,1]
    
    # Assign log-odds to trials
    logoddspre = NULL
    logoddspost = NULL
    for (t in 1:length(precoh)) {
      if (accuracy[t] == 1){
        logoddspre[t] = loglikpre_cor[precoh[t], postcoh[t]]
        logoddspost[t] = loglikpost_cor[precoh[t], postcoh[t]]
      } else {
        logoddspre[t] = loglikpre_err[precoh[t], postcoh[t]]
        logoddspost[t] = loglikpost_err[precoh[t], postcoh[t]]
      }
      
    }
    
    QSR = scale(QSR)
    rawConf = conf
    conf = scale(conf)
    subj = rep(s, length(accuracy))
    subData = data.frame("subj"=subj, "BOLD"=BOLD, "precoh"=precoh-2, "postcoh"=postcoh-2, "logRT"=logRT, "logConfRT"=logConfRT, "conf"=conf, "accuracy"=accuracy, "QSR"=QSR, "rawConf"=rawConf, "logoddspre"=logoddspre, "logoddspost"=logoddspost)
    
    bigData = rbind(bigData, subData)
    
  }
  
  bigData$subj <- factor(bigData$subj)
  bigData$accuracy <- factor(bigData$accuracy)

  # Fit regression models

  ## COHERENCE MODEL
  postcohModel = lmer(BOLD ~ precoh*accuracy + postcoh*accuracy + logRT + (1 + precoh*accuracy + postcoh*accuracy + logRT | subj), data=bigData, control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
  fix <- fixef(postcohModel)
  fix.se <- sqrt(diag(vcov(postcohModel)))
  barX <- barplot(fix, beside = TRUE, ylab = "Regression betas", names.arg = names(fix), cex.names = 0.5)
  errorBar(barX, fix, fix.se)
  betas <- c(fix, fix.se)
  # write.csv(betas, file = paste('regression_coh_', rois[r], '.csv', sep=""))

  # P-values and estimates
  print(paste("ROI analysis for ", rois[r], sep=""))
  print(summary(postcohModel))
  print(Anova(postcohModel, type=3))

  # Store as individual model for plotting in table
  eval(parse(text = paste("postcohModel", r, "= postcohModel", sep="")))
  
  ## LOGODDS REGRESSION FROM FIT OF BAYESIAN MODEL
  logoddsModel = lmer(BOLD ~ logoddspre + logoddspost + logRT + (1 + logoddspre + logoddspost + logRT | subj), data=bigData, control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
  fix <- fixef(logoddsModel)
  fix.se <- sqrt(diag(vcov(logoddsModel)))
  barX <- barplot(fix, beside = TRUE, ylab = "Regression betas", names.arg = names(fix), cex.names = 0.5)
  errorBar(barX, fix, fix.se)
  betas <- c(fix, fix.se)
  
  # P-values and estimates
  print(summary(logoddsModel))
  print(Anova(logoddsModel, type=3))

  # Store as individual model for plotting in table
  eval(parse(text = paste("logoddsModel", r, "= logoddsModel", sep="")))
  
}

# Plot all regressions in one big table
stargazer(postcohModel1, postcohModel2, postcohModel3, postcohModel4, postcohModel5, postcohModel6, star.cutoffs = c(0.05, 0.01, 0.001), keep.stat = c("n"), dep.var.labels = c("BOLD"), column.labels = c("pMFC", "area 46", "FPl", "FPm", "vmPFC", "v. striatum"), covariate.labels = c("accuracy", "pre-decision coherence", "post-decision coherence", "pre*accuracy", "post*accuracy", "log(RT)"), order = c(2, 1, 3, 5, 6, 4), report = ('vc*s'),  digits = 2, digits.extra = 1)
