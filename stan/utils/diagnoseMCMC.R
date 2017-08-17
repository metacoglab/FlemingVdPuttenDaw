diagnoseMCMC <- function(fit, tplot=F){
  #Provide diagnostics for a stanfit object. This function checks for
  #divergence errors, values of Rhat > 1.01 (potential scale reduction
  #factor: Gelman and Rubin 92), and, optionally, traceplots for all
  #parameters. Returns extracted samples.
  
  samplerParams <- get_sampler_params(fit, inc_warmup=FALSE)
  
  sumDiv <- c()
  for (i in 1:length(samplerParams)) {
    sumDiv <- c(sumDiv, sum(samplerParams[[1]][,"divergent__"]))
  }
  
  if (any(as.logical(sumDiv))) {
    print("divergence errors !!! try reparameterizing your model and/or decreasing step size.")
  }
  
  rhats <- summary(fit)$summary[,"Rhat"]
  
  if (any(rhats>1.1, na.rm = TRUE)) {
    print("Rhats > 1.1 !!! try longer chains.")
  }
  
  numSamps <- length(extract(fit, pars = "lp__")[[1]])
  effPer <- summary(fit)$summary[,"n_eff"] / numSamps
  
  if (any(effPer<0.1)) {
    print("params with effective samples < 10 percent of total !!! try longer chains or thinning.")
  }
  
  if (tplot) {
    traceplot(fit, ask = T, inc_warmup=F)
  }
  
  samp <- extract(fit)
  return(samp)
  
}