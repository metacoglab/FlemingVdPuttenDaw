# Behavioural regression models for Fleming, van der Putten & Daw
# For more details of setting lme4 options see https://cran.r-project.org/web/packages/lme4/vignettes/lmerperf.html
#
# Steve Fleming 2016 stephen.fleming@ucl.ac.uk 

rm(list=ls())
require(R.matlab) 
require(lme4)
require(car)
require(optimx)
require(stargazer)
options(contrasts = c("contr.treatment", "contr.poly")) # This is R defaults but set it anyway to be safe
source("~/Dropbox/Published/FlemingNN2018/Code&analysis/analysis/errorBar.R")

for (d in 1:2) {
  
  if (d == 1) {
    dataset = "2"
  } else {
    dataset = "3"
  }
  
  if (dataset == "2") {
    dataDir = "~/Dropbox/Published/FlemingNN2018/Code&analysis/task/locFullData/"
    filePrefix = "fMRI_pilotData_sub_"
    suffix = ""
    subjects = c(seq(12,28), seq(30,37))
  } else {
    dataDir = "~/Dropbox/Published/FlemingNN2018/Code&analysis/task/locFullData/"
    filePrefix = "fMRI_pilotData_sub_"
    suffix = "_fMRI"
    subjects =c(seq(12,19), seq(23,28), seq(30,37))
  }
  
  # Extract out data from mat files
  bigData = NULL
  for (s in 1:length(subjects)){
    
    ## EXTRACT DATA FILE
    setwd(dataDir)
    DATA = readMat(paste(filePrefix,subjects[s],suffix,'_2.mat',sep=""))
    
    dat = DATA$locDATA
    
    precoh = NULL
    postcoh = NULL
    RT = NULL
    conf = NULL
    confRT = NULL
    response = NULL
    dir = NULL
    precoh = dat[,,1]$dots.coherence[1,]
    postcoh = dat[,,1]$post.coherence[1,]
    RT = dat[,,1]$reaction.time.button[1,]
    conf = dat[,,1]$mouse.response[1,]
    confRT = dat[,,1]$reaction.time.mouse[1,]
    response = dat[,,1]$button.response[1,] - 1
    response[response == 0] = -1
    dir = dat[,,1]$dots.direction[1,]/360
    dir[dir == 0.5] = -1
    accuracy = response == dir
    
    logRT = scale(log(RT))
    logConfRT = scale(log(confRT))
    subj = rep(s, length(accuracy))
    subData = data.frame("subj"=subj, "dir"=dir, "precoh"=precoh, "postcoh"=postcoh, "conf"=conf, "logConfRT"=logConfRT, "response"=response, "dir"=dir, "logRT"=logRT, "accuracy"=accuracy)
    
    bigData = rbind(bigData, subData)
    
  }
  
  # Fit regression models
  setwd("~/Dropbox/Published/FlemingNN2018/Code&analysis/analysis/")
  bigData$subj <- factor(bigData$subj)
  
  ## Accuracy
  accModel = glmer(accuracy ~ precoh + (1 + precoh |subj), data=bigData, family="binomial")
  fix <- fixef(accModel)
  fix.se <- sqrt(diag(vcov(accModel)))
  barX <- barplot(fix, beside = TRUE, ylab = "Regression betas for accuracy", names.arg = names(fix), cex.names = 0.5)
  errorBar(barX, fix, fix.se)
  betas <- c(fix, fix.se)
  write.csv(betas, file = 'regression_accuracy.csv')
  
  print(paste("Performance analysis for dataset ", dataset, sep=""))
  print(summary(accModel))
  hist(resid(accModel))
  
  # Store as individual model for plotting in table
  eval(parse(text = paste("accModel", d, "= accModel", sep="")))
  
  # stargazer(accModel, star.cutoffs = c(0.05, 0.01, 0.001), keep.stat = c("n"), dep.var.labels = c("Behav accuracy", "fMRI accuracy"), covariate.labels = c("pre-decision coherence"), report = ('vc*sp'))
  
  ## Confidence
  bigData_correct <- bigData[bigData$accuracy == 1, ]
  bigData_error <- bigData[bigData$accuracy == 0, ]
  
  confModel_correct = lmer(conf ~ precoh*postcoh + logRT + (1 +  precoh*postcoh + logRT |subj), data=bigData_correct
                           , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
  fix <- fixef(confModel_correct)
  fix.se <- sqrt(diag(vcov(confModel_correct)))
  barX <- barplot(fix, beside = TRUE, ylab = "Regression betas for confidence", names.arg = names(fix), cex.names = 0.5, ylim = c(-2,2))
  errorBar(barX, fix, fix.se)
  betas <- c(fix, fix.se)
  write.csv(betas, file = 'regression_correct.csv')
  
  # P-values and estimates
  print(paste("Confidence-correct analysis for dataset ", dataset, sep=""))
  print(summary(confModel_correct))
  print(Anova(confModel_correct, type=3))
  hist(resid(confModel_correct))
  
  # Store as individual model for plotting in table
  eval(parse(text = paste("confModel_correct", d, "= confModel_correct", sep="")))
  
  confModel_error = lmer(conf ~ precoh*postcoh + logRT + (1 +  precoh*postcoh + logRT|subj), data=bigData_error
                         , control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "bobyqa", starttests = FALSE, kkt = FALSE)))
  fix <- fixef(confModel_error)
  fix.se <- sqrt(diag(vcov(confModel_error)))
  barX <- barplot(fix, beside = TRUE, ylab = "Regression betas for confidence", names.arg = names(fix), cex.names = 0.5, ylim = c(-2,2))
  errorBar(barX, fix, fix.se)
  betas <- c(fix, fix.se)
  write.csv(betas, file = 'regression_error.csv')
  
  # P-values and estimates
  print(paste("Confidence-error analysis for dataset ", dataset, sep=""))
  print(summary(confModel_error))
  print(Anova(confModel_error, type=3))
  hist(resid(confModel_error))
  
  # Store as individual model for plotting in table
  eval(parse(text = paste("confModel_error", d, "= confModel_error", sep="")))
  
}

stargazer(accModel1, accModel2, star.cutoffs = c(0.05, 0.01, 0.001), keep.stat = c("n"), dep.var.labels = c("Accuracy"), column.labels = c("Behaviour", "fMRI"), covariate.labels = c("pre-decision coherence"), report = ('vc*s'),  digits = 2, digits.extra = 1)
stargazer(confModel_correct1, confModel_error1, confModel_correct2, confModel_error2, star.cutoffs = c(0.05, 0.01, 0.001), keep.stat = c("n"), dep.var.labels = c("Confidence"), column.labels = c("Behav. correct", "Behav. error","fMRI correct","fMRI error"), covariate.labels = c("pre-decision coherence", "post-decision coherence", "log(RT)", "pre*post coherence"), report = ('vc*s'),  digits = 2, digits.extra = 1)

