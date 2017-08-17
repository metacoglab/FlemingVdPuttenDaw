data {
  int Ns;
  int N;
  int ST[Ns];
  int a[Ns,N];
  vector[N] d[Ns];
  vector[N] theta1[Ns];
  vector[N] theta2[Ns];
  vector[3] coh[Ns];
  vector[N] conf[Ns];
}
parameters {
  real<lower=0> mu_k1;
  real<lower=0,upper=10> sd_k1;
  real mu_m;
  real<lower=0,upper=10> sd_m;
  real<lower=0> alpha_w;
  real<lower=0> beta_w;
  real<lower=0, upper=1> w[Ns];
  real<lower=0> k1[Ns];
  real m[Ns];
  vector[N] x1[Ns];
  vector[N] x2[Ns];
}
transformed parameters {
vector<lower=0,upper=1>[N] model_conf[Ns];
for (s in 1:Ns) {
  real muT;
  real varT;
  real temp_varT[3];
  real kTheta[3];
  for (c in 1:3) {
    kTheta[c] <- k1[s] * coh[s,c];
  }
  muT <- sum(kTheta)/3;
  for (c in 1:3) {
    temp_varT[c] <- pow((kTheta[c] - muT), 2);
  }
  varT <- sum(temp_varT)/3 + 1;
   
  for (i in 1:N) {
    if (i <= ST[s]) {   
        real loglikdir_pre;
        real loglikdir_post;
        real loglikC_pre;
        real loglikC_post;
        real loglikC;
        loglikdir_pre <- (2 * muT * x1[s,i])/varT;
        loglikdir_post <- (2 * muT * x2[s,i])/varT;
        if (a[s,i] == 1) {
          loglikC_pre <- loglikdir_pre;
          loglikC_post <- loglikdir_post;
        } else {
          loglikC_pre <- -loglikdir_pre;
          loglikC_post <- -loglikdir_post;
        }  
        loglikC <- (2*w[s])*loglikC_pre + (2-(2*w[s]))*loglikC_post;
        model_conf[s,i] <- inv_logit(loglikC);
          
      // if missing data enter arbitrary confidence value
      } else {
        model_conf[s,i] <- 0.5;
      }
  }
}
}
model {
  // hyperparameters
  mu_k1 ~ normal(0, 10);
  mu_m ~ normal(0, 10);

  for (s in 1:Ns) {
    // priors on k
    k1[s] ~ normal(mu_k1, sd_k1);
    m[s] ~ normal(mu_m, sd_m);
    w[s] ~ beta(alpha_w, beta_w);
    
    for (i in 1:N) {
      if (i <= ST[s]) {   
        x1[s,i] ~ normal(d[s,i].*theta1[s,i]*k1[s], 1);
        x2[s,i] ~ normal(d[s,i].*theta2[s,i]*k1[s], 1);
        a[s,i] ~ bernoulli_logit(100*(x1[s,i]-m[s]));
        conf[s,i] ~ normal(model_conf[s,i], 0.025);
      }
    }
  }
}