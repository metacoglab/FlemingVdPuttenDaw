rm(list=ls())
require(rstan)
baseDir = "~/Dropbox/Research/Metacognition/stateactionexpt/github/stan/modelfits/"

models = c('ideal', 'weighted', 'accweighted', 'mapping') # selection from: ideal, weighted, accweighted, mapping

# Loop over and load models, write out variance parameters
sumrhats = NA
rhat <- matrix(list(), 1, 8)
n_eff <- matrix(list(), 1, 8)
j = 1
for (d in 1:2) {
  
  for (m in 1:length(models)) {
    
    if (d == 1) {
      dataset = 2
    } else {
      dataset = 3
    }
    
    outputPath <- paste(baseDir, models[m], "/", sep="")
    
    # Load fitted model
    load(paste(outputPath, "fitObject_group", dataset, ".RData", sep=""))
    
    # Load vector of fitted subject IDs
    subjects = read.csv(paste(outputPath, "subParams_subjects_dataset", dataset, ".csv",sep=""))
    subjectCode = subjects$x
    subjectID = subjects$X
    
    # Get summary of fit object amd rhats as matrix, each rhat[[i]] entry is one dataset/model
    params = summary(meta_fit)$summary
    n_eff[[j]] = params[,"n_eff"]
    rhat[[j]] = params[,"Rhat"]
    sumrhats[j] = sum(rhat[[j]] > 1.1)
    j=j+1
    
    # Extract out subject specific parameter estimates and Rhats () for each model type
    k1_sd = NA
    m_sd = NA
    w_sd = NA
    gamma_sd = NA
    for (i in 1:length(subjectID)) {
      tempInd = paste("k1[", subjectID[i], "]", sep="")
      k1_sd[i] = params[tempInd,3]
      tempInd = paste("m[", subjectID[i], "]", sep="")
      m_sd[i] = params[tempInd,3]
      if (models[m] == 'weighted' | models[m] == 'weightedmapping' | models[m] == 'accweighted' | models[m] == 'accweightedmapping') {
        tempInd = paste("w[", subjectID[i], "]", sep="")
        w_sd[i] = params[tempInd,3]        
      }
      if (models[m] == 'mapping' | models[m] == 'weightedmapping' | models[m] == 'accweightedmapping') {
        tempInd = paste("gamma[", subjectID[i], "]", sep="")
        gamma_sd[i] = params[tempInd,3]    
      }
    }
    
    # Write out subject level parameters as CSV for matlab
    write.csv(k1_sd, paste(outputPath, "subParams_k1_sd_dataset", dataset, ".csv",sep=""))
    write.csv(m_sd, paste(outputPath, "subParams_m_sd_dataset", dataset, ".csv",sep=""))
    if (models[m] == 'weighted' | models[m] == 'weightedmapping' | models[m] == 'accweighted' | models[m] == 'accweightedmapping') {
      write.csv(w_sd, paste(outputPath, "subParams_w_sd_dataset", dataset, ".csv",sep=""))
    }
    if (models[m] == 'mapping' | models[m] == 'weightedmapping' | models[m] == 'accweightedmapping') {
      write.csv(gamma_sd, paste(outputPath, "subParams_gamma_sd_dataset", dataset, ".csv",sep=""))
    }

  }
}