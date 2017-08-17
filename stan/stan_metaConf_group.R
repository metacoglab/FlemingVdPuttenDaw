rm(list=ls())
require(R.matlab)
require(rstan)
require(parallel)
source('~/Dropbox/Utils/stan/stanUtilities.R')
baseDir = "~/Dropbox/Research/Metacognition/stateactionexpt/github/stan/modelfits/"

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

models = c('ideal') # selection from one or more of: ideal, weighted, accweighted, mapping
dataset = readline(prompt = "Dataset? 1=behav pilot, 2=behav, 3=fmri ")

if (dataset == "2") {
  dataDir = "~/Dropbox/Research/Metacognition/stateactionexpt/github/data/"
  filePrefix = "fMRI_pilotData_sub_"
  suffix = ""
  subjects = c(seq(12,28), seq(30,37))
  Ntrials <- 900
} else {
  dataDir = "~/Dropbox/Research/Metacognition/stateactionexpt/github/data/"
  filePrefix = "fMRI_pilotData_sub_"
  suffix = "_fMRI"
  subjects =c(seq(12,19), seq(23,28), seq(30,37))
  Ntrials <- 360
}

Nchains <- 3
Niter <- 3000
Nburn <- 1000

# GET DATA
a_matrix <- matrix(0, length(subjects), Ntrials)
d_matrix <- matrix(0, length(subjects), Ntrials)
conf_matrix <- matrix(0, length(subjects), Ntrials)
theta1_matrix <- matrix(0, length(subjects), Ntrials)
theta2_matrix <- matrix(0, length(subjects), Ntrials)
allCoh = matrix(0, length(subjects), 3)

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
  dir = NULL
  precoh = dat[,,1]$dots.coherence[1,]
  postcoh = dat[,,1]$post.coherence[1,]
  conf = dat[,,1]$mouse.response[1,]
  response = dat[,,1]$button.response[1,] - 1
  dir = dat[,,1]$dots.direction[1,]/360
  dir[dir == 0.5] = -1
  
  a_matrix[subi,] <- response
  d_matrix[subi,] <- dir
  conf_matrix[subi,] <- conf
  theta1_matrix[subi,] <- precoh
  theta2_matrix[subi,] <- postcoh
  allCoh[subi,] <- unique(precoh)
}

## Handle ragged matrices due to missed trials (only an issue for fMRI data)
goodTrials <- NA
for (subi in 1:length(subjects)) {
  missed <- conf_matrix[subi,] == "NaN"
  goodTrials[subi] <- Ntrials - sum(missed)
  a_matrix[subi,] <- c(a_matrix[subi,!missed], rep(0,sum(missed)))
  d_matrix[subi,] <- c(d_matrix[subi,!missed], rep(0,sum(missed)))
  conf_matrix[subi,] <- c(conf_matrix[subi,!missed], rep(0,sum(missed)))
  theta1_matrix[subi,] <- c(theta1_matrix[subi,!missed], rep(0,sum(missed)))
  theta2_matrix[subi,] <- c(theta2_matrix[subi,!missed], rep(0,sum(missed)))
}

## FIT MODELS
for (m in 1:length(models)) {
  
  setwd("~/Dropbox/Research/Metacognition/stateactionexpt/analysis/stan/")
  
  outputPath <- paste(baseDir, models[m], "/", sep="")
  data <- list(Ns=length(subjects), N=Ntrials, ST=goodTrials, a=a_matrix, d=d_matrix, conf=conf_matrix, theta1=theta1_matrix, theta2=theta2_matrix, coh=allCoh)
  
  # Fit model
  seed <- 12345
  # Compile model object
  # Fit model
  seed <- 12345
  # Compile and initialize parameter values
  if (models[m] == 'ideal') {
    parameters <- c("mu_k1", "sd_k1", "k1", "mu_m", "sd_m", "m","log_lik")
    group_param <- c("mu_k1", "sd_k1","mu_m","sd_m")
  } else if (models[m] == 'weighted' || models[m] == 'accweighted') {
    parameters <- c("mu_k1", "sd_k1", "k1", "mu_m", "sd_m", "m", "alpha_w", "beta_w", "w", "log_lik")
    group_param <- c("mu_k1", "sd_k1", "mu_m", "sd_m", "alpha_w", "beta_w")
  } else if (models[m] == 'mapping') {
    parameters <- c("mu_k1", "sd_k1", "k1", "mu_m", "sd_m", "m", "sd_gamma", "mu_gamma", "gamma", "log_lik")
    group_param <- c("mu_k1", "sd_k1", "mu_m", "sd_m", "sd_gamma", "mu_gamma")
  } else if (models[m] == 'weightedmapping' || models[m] == 'accweightedmapping') {
    parameters <- c("mu_k1", "sd_k1", "k1", "mu_m", "sd_m", "m", "alpha_w", "beta_w", "w", "sd_gamma", "mu_gamma", "gamma", "log_lik")
    group_param <- c("mu_k1", "sd_k1", "mu_m", "sd_m", "alpha_w", "beta_w", "sd_gamma", "mu_gamma")
  }
  
  # Compile and run
  meta_fit <- stan(paste("metaConf_fit_group_",models[m],'.stan',sep=""), seed = seed,data = data, pars = parameters, warmup = Nburn, iter = Niter, chains = Nchains)
  
  # Extract posterior means and quantiles
  samp <- diagnoseMCMC(meta_fit)
  rhats <- summary(meta_fit)$summary[,"Rhat"]

  # Extract posterior means and quantiles for group-level parameters
  MAP <- NULL
  upper_CI <- NULL
  lower_CI <- NULL
  q95 <- c(0.025,0.975)
  for (i in 1:length(group_param)){
    # Group-level parameters
    sample_vector = as.numeric(samp[group_param[i]][[1]])
    MAP[parameters[i]] <- density(sample_vector)$x[which(density(sample_vector)$y ==
                                                           max(density(sample_vector)$y))]
    upper_CI[parameters[i]] <- quantile(sample_vector,probs=q95[2])
    lower_CI[parameters[i]] <- quantile(sample_vector,probs=q95[1])
    
  }
  
  sub_param <- NULL
  for (subi in 1:length(subjects)) {
    sub_param$k1[subi] <- density(samp$k1[,subi])$x[which(density(samp$k1[,subi])$y ==
                                                            max(density(samp$k1[,subi])$y))]
    sub_param$m[subi] <- density(samp$m[,subi])$x[which(density(samp$m[,subi])$y ==
                                                          max(density(samp$m[,subi])$y))]
    if (models[m] == 'weighted' | models[m] == 'weightedmapping' | models[m] == 'accweighted' | models[m] == 'accweightedmapping') {
      sub_param$w[subi] <- density(samp$w[,subi])$x[which(density(samp$w[,subi])$y ==
                                                            max(density(samp$w[,subi])$y))]
    }
    if (models[m] == 'mapping' | models[m] == 'weightedmapping' | models[m] == 'accweightedmapping') {
      sub_param$gamma[subi] <- density(samp$gamma[,subi])$x[which(density(samp$gamma[,subi])$y ==
                                                                    max(density(samp$gamma[,subi])$y))]
    }
  }
  
  # Put in list and write to disk
  fitParams <- list(sub_param=sub_param, group_param_names=group_param, MAP=MAP, upper_CI=upper_CI, lower_CI=lower_CI,
                    rhats=rhats)
  save(fitParams, file = paste(outputPath, "fitParams_group", dataset, ".RData",sep=""))
  save(meta_fit, file = paste(outputPath, "fitObject_group", dataset, ".RData",sep=""))
  
  # Write out .mat file for easy access in matlab
  writeMat(paste(outputPath, models[m], '_dataset', dataset,'_fit.mat',sep=""), group_param_names=group_param, MAP=MAP, upper_CI=upper_CI, lower_CI=lower_CI, sub_param=sub_param, rhats=rhats)
  
  # Write out subject level parameters as CSV for matlab
  write.csv(subjects, paste(outputPath, "subParams_subjects_dataset", dataset, ".csv",sep=""))
  write.csv(sub_param$k1, paste(outputPath, "subParams_k1_dataset", dataset, ".csv",sep=""))
  write.csv(sub_param$m, paste(outputPath, "subParams_m_dataset", dataset, ".csv",sep=""))
  if (models[m] == 'weighted' | models[m] == 'weightedmapping' | models[m] == 'accweighted' | models[m] == 'accweightedmapping') {
    write.csv(sub_param$w, paste(outputPath, "subParams_w_dataset", dataset, ".csv",sep=""))
  }
  if (models[m] == 'mapping' | models[m] == 'weightedmapping' | models[m] == 'accweightedmapping') {
    write.csv(sub_param$gamma, paste(outputPath, "subParams_gamma_dataset", dataset, ".csv",sep=""))
  }
}

