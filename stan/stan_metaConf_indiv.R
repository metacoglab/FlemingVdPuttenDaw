rm(list=ls())
require(R.matlab)
require(rstan)
require(shinystan)
require(parallel)
require(loo)
source('~/Dropbox/Utils/stan/stanUtilities.R')
baseDir = "~/Dropbox/Research/Metacognition/stateactionexpt/analysis/stan/modelfits/indiv/"

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

models = c('ideal', 'weighted', 'choiceweighted', 'choicebias', 'mapping', 'ideal_dt', 'weighted_dt', 'choiceweighted_dt', 'choicebias_dt', 'mapping_dt')
dataset = readline(prompt = "Dataset? 1=behav pilot, 2=behav, 3=fmri ")

if (dataset == "1") {
  dataDir = "~/Dropbox/Research/Metacognition/stateactionexpt/task/locData/"
  filePrefix = "locData_sub_"
  suffix = ""
  subjects = seq(3,17)
  Ntrials <- 900
} else if (dataset == "2") {
  dataDir = "~/Dropbox/Research/Metacognition/stateactionexpt/task/locFullData/"
  filePrefix = "fMRI_pilotData_sub_"
  suffix = ""
  subjects = c(seq(12,28), seq(30,37))
  Ntrials <- 900
} else {
  dataDir = "~/Dropbox/Research/Metacognition/stateactionexpt/task/locFullData/"
  filePrefix = "fMRI_pilotData_sub_"
  suffix = "_fMRI"
  subjects =c(seq(12,19), seq(23,28), seq(30,37))
  Ntrials <- 360
}

Nchains <- 3
Niter <- 4000
Nburn <- 1000

a_matrix <- NA
d_matrix <- NA
conf_matrix <- NA
theta1_matrix <- NA
theta2_matrix <- NA
rt_matrix <- NA
allCoh <- NA

# Get data for each subject and fit
for (subi in 1:length(subjects)){
  
  ## EXTRACT DATA FILE
  setwd(dataDir)
  DATA = readMat(paste(filePrefix,subjects[subi],suffix,'_2.mat',sep=""))
  
  dat = DATA$locDATA
  
  precoh = NULL
  precoh_index = NULL
  postcoh = NULL
  postcoh_index = NULL
  conf = NULL
  response = NULL
  rt = NULL
  dir = NULL
  precoh = dat[,,1]$dots.coherence[1,]
  postcoh = dat[,,1]$post.coherence[1,]
  conf = dat[,,1]$mouse.response[1,]
  response = dat[,,1]$button.response[1,] - 1
  rt = dat[,,1]$reaction.time.button
  dir = dat[,,1]$dots.direction[1,]/360
  dir[dir == 0.5] = -1
  
  a_matrix <- response
  d_matrix <- dir
  conf_matrix <- conf
  theta1_matrix <- precoh
  theta2_matrix <- postcoh
  rt_matrix <- rt
  allCoh <- unique(precoh)
  
  ## Handle ragged matrices due to missed trials (only an issue for fMRI data)
  goodTrials <- NA
  
  missed <- conf_matrix == "NaN"
  goodTrials <- Ntrials - sum(missed)
  a_matrix <- c(a_matrix[!missed], rep(0,sum(missed)))
  d_matrix <- c(d_matrix[!missed], rep(0,sum(missed)))
  conf_matrix <- c(conf_matrix[!missed], rep(0,sum(missed)))
  theta1_matrix <- c(theta1_matrix[!missed], rep(0,sum(missed)))
  theta2_matrix <- c(theta2_matrix[!missed], rep(0,sum(missed)))
  rt_matrix <- c(rt_matrix[!missed], rep(0,sum(missed)))

  ## FIT MODELS
  for (m in 1:length(models)) {
    
    setwd("~/Dropbox/Research/Metacognition/stateactionexpt/analysis/stan/")
    
    outputPath <- paste(baseDir, models[m], "/", sep="")
    data <- list(N=Ntrials, ST=goodTrials, a=a_matrix, d=d_matrix, conf=conf_matrix, theta1=theta1_matrix, theta2=theta2_matrix, rt=rt_matrix, coh=allCoh)
    
    # Fit model
    seed <- 12345
    # Compile and initialize parameter values
    if (models[m] == 'ideal') {
      parameters <- c("k1", "m")
    } else if (models[m] == 'choicebias') {
      parameters <- c("k1", "m", "w")
    } else if (models[m] == 'weighted' | models[m] == 'choiceweighted') {
      parameters <- c("k1", "m", "w1", "w2")
    } else if (models[m] == 'mapping') {
      parameters <- c("k1", "m", "gamma")
    } else if (models[m] == 'ideal_dt') {
      parameters <- c("k1", "m", "brt")      
    } else if (models[m] == 'choicebias_dt') {
      parameters <- c("k1", "m", "w", "brt")
    } else if (models[m] == 'weighted_dt' | models[m] == 'choiceweighted_dt') {
      parameters <- c("k1", "m", "w1", "w2", "brt")
    } else if (models[m] == 'mapping_dt') {
      parameters <- c("k1", "m", "gamma", "brt")
    }
    
    # Compile and run
    meta_fit <- stan(paste("metaConf_fit_indiv_",models[m],'.stan',sep=""), seed = seed,data = data, pars = parameters, warmup = Nburn, iter = Niter, chains = Nchains)
    
    # Extract posterior means and quantiles
    samp <- diagnoseMCMC(meta_fit)
    rhats <- summary(meta_fit)$summary[,"Rhat"]
    params = summary(meta_fit)$summary
    
    # Store MAP and upper/lower CI for each subject separately
    sub_param <- NULL
    sub_param$k1 <- density(samp$k1)$x[which(density(samp$k1)$y ==
                                                      max(density(samp$k1)$y))]
    sub_param$m <- density(samp$m)$x[which(density(samp$m)$y ==
                                                    max(density(samp$m)$y))]
    sub_param$k1_sd = params["k1",3]
    sub_param$m_sd = params["m",3]
    if (models[m] == 'choicebias' | models[m] == 'choicebias_dt') {
      sub_param$w <- density(samp$w)$x[which(density(samp$w)$y ==
                                                      max(density(samp$w)$y))]
      sub_param$w_sd = params["w",3]        
    }
    if (models[m] == 'weighted' | models[m] == 'choiceweighted' | models[m] == 'weighted_dt' | models[m] == 'choiceweighted_dt') {
      sub_param$w1 <- density(samp$w1)$x[which(density(samp$w1)$y ==
                                               max(density(samp$w1)$y))]
      sub_param$w1_sd = params["w1",3]       
      sub_param$w2 <- density(samp$w2)$x[which(density(samp$w2)$y ==
                                               max(density(samp$w2)$y))]
      sub_param$w2_sd = params["w2",3]        
    }
    if (models[m] == 'mapping' | models[m] == 'mapping_dt') {
      sub_param$gamma <- density(samp$gamma)$x[which(density(samp$gamma)$y ==
                                                              max(density(samp$gamma)$y))]
      sub_param$gamma_sd = params["gamma",3]    
    }
    if (models[m] == 'ideal_dt' | models[m] == 'weighted_dt' | models[m] == 'choicebias_dt' | models[m] == 'choiceweighted_dt' | models[m] == 'mapping_dt') {
      sub_param$brt <- density(samp$brt)$x[which(density(samp$brt)$y ==
                                                       max(density(samp$brt)$y))]
      sub_param$brt_sd = params["brt",3]    
    }

    fitParams <- list(sub_param=sub_param, rhats=rhats)
    save(fitParams, file = paste(outputPath, "fitParams_indiv", subjects[subi], '_dataset', dataset, ".RData",sep=""))
    save(meta_fit, file = paste(outputPath, "fitObject_indiv", subjects[subi], '_dataset', dataset, ".RData",sep=""))
    
    # Write out .mat file for easy access in matlab
    writeMat(paste(outputPath, models[m], '_sub', subjects[subi], '_dataset', dataset,'_fit.mat',sep=""), sub_param=sub_param, rhats=rhats)
    
  }
}

setwd("~/Dropbox/Research/Metacognition/stateactionexpt/analysis/stan/")
