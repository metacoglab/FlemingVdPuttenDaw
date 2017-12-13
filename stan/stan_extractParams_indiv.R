rm(list=ls())
require(rstan)
baseDir = "~/Dropbox/Research/Metacognition/stateactionexpt/analysis/stan/modelfits/indiv/"

models = c('ideal', 'weighted', 'choiceweighted', 'choicebias', 'mapping', 'ideal_dt', 'weighted_dt', 'choiceweighted_dt', 'choicebias_dt', 'mapping_dt')

# Loop over and load models, write out variance parameters
rhat <- matrix(list(), 1, 10)
n_eff <- matrix(list(), 1, 10)
j = 1
for (d in 1:2) {
  if (d == 1) {
    dataset = 2
    subjects = c(seq(12,28), seq(30,37))
  } else {
    dataset = 3
    subjects =c(seq(12,19), seq(23,28), seq(30,37))
  }
  
  eval(parse(text=paste("sumrhats", d, " <- matrix(list(), length(models), length(subjects))", sep="")))
  eval(parse(text=paste("ndraws", d, " <- matrix(list(), length(models), length(subjects))", sep="")))
  
  for (m in 1:length(models)) {
    
    for (s in 1:length(subjects)) {

      outputPath <- paste(baseDir, models[m], "/", sep="")
      
      # Load fitted model
      load(paste(outputPath, "fitObject_indiv", subjects[s], "_dataset", dataset, ".RData", sep=""))
      
      # Get summary of fit object amd rhats as matrix, each rhat[[i]] entry is one dataset/model
      params = summary(meta_fit)$summary
      n = dim(params)
      n_eff[[j]] = params[1:n[1]-1,"n_eff"]
      rhat[[j]] = params[1:n[1]-1,"Rhat"]
      eval(parse(text=paste("sumrhats", d, "[m,s] = sum(rhat[[j]] > 1.1)", sep="")))
      
      list_of_draws = extract(meta_fit)
      eval(parse(text=paste("ndraws", d, "[m,s] = dim(list_of_draws$k1)", sep="")))
      
      j=j+1

    }
    
  }
}

# Print subjects/models with out of bounds rhats:
which(sumrhats1!=0,arr.ind = T)
which(sumrhats2!=0,arr.ind = T)

which(ndraws1<9000,arr.ind = T)
which(ndraws2<9000,arr.ind = T)