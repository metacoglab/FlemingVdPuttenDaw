data {
  int N;
  int ST;
  int a[N];
  vector[N] d;
  vector[N] theta1;
  vector[N] theta2;
  vector[3] coh;
  vector[N] conf;
  vector[N] rt;
}
parameters {
  real<lower=0> k1;
  real m;
  real brt;
  real<lower=0> gamma;
  vector[N] x1;
  vector[N] x2;
}
transformed parameters {
vector<lower=0,upper=1>[N] model_conf;

  real muT;
  real varT;
  real temp_varT[3];
  real kTheta[3];
  
  for (c in 1:3) {
    kTheta[c] = k1 * coh[c];
  }
  muT = sum(kTheta)/3;
  for (c in 1:3) {
    temp_varT[c] = pow((kTheta[c] - muT), 2);
  }
  varT = sum(temp_varT)/3 + 1;
   
  for (i in 1:N) {
    if (i <= ST) {   
        real loglikdir_pre;
        real loglikdir_post;
        real loglikC_pre;
        real loglikC_post;
        real loglikC;
        real temp_conf;
        real LO_pi_conf;
        loglikdir_pre = (2 * muT * x1[i])/varT;
        loglikdir_post = (2 * muT * x2[i])/varT;
        if (a[i] == 1) {
          loglikC_pre = loglikdir_pre;
          loglikC_post = loglikdir_post;
        } else {
          loglikC_pre = -loglikdir_pre;
          loglikC_post = -loglikdir_post;
        }  
        loglikC = loglikC_pre + loglikC_post + brt*log(rt[i]);
        temp_conf = inv_logit(loglikC);
        
        LO_pi_conf = gamma*(log(temp_conf/(1-temp_conf)));
        model_conf[i] = inv_logit(LO_pi_conf);
      // if missing data enter arbitrary confidence value
      } else {
        model_conf[i] = 0.5;
      }
  }
}
model {

    // priors on k
    k1 ~ normal(0, 10);
    m ~ normal(0, 10);
    gamma ~ normal(1, 1);
    brt ~ normal(0, 10);
    
    for (i in 1:N) {
      if (i <= ST) {   
        x1[i] ~ normal(d[i].*theta1[i]*k1, 1);
        x2[i] ~ normal(d[i].*theta2[i]*k1, 1);
        a[i] ~ bernoulli_logit(100*(x1[i]-m));
        conf[i] ~ normal(model_conf[i], 0.025);
    }
  }
}
