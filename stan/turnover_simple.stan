data {
  int N;  // sample size of taxa
  int T;  // sample size of temporal units
  int D;  // number of indiv-level predictors
  int U;  // number of group-level predictors
  int P;  // number of floral phases

  int sight[N, T];  // observed presence
  row_vector[D] x[N];  // matrix of indiv-level covariates
  row_vector[U] u[T];  // matrix of group-level covariates

  int phase[T];  // plant phase
}
parameters {
  vector[T] intercept;  // intercept varies by what point in time
  real intercept_mu;  // mean intercept
  real<lower=0> sigma;  // variance of varying-intercept
  
  matrix[T, D] beta_std;
  real beta_mu[D];
  real<lower=0> beta_sigma[D];
  
  matrix[T, U] alpha_std;  // effect of group-level covariates
  real alpha_mu[U];
  real<lower=0> alpha_sigma[U];

  vector[P] eff_phase;  // effect of plant-phase (mean 0)
  real<lower=0> scale_phase;  // variance in plant-phase effect
}
transformed parameters {
  matrix[N, T] pred; 
  matrix[T, D] beta;  // effect of indiv-level covariates
  matrix[T, U] alpha; 
  
  // noncentered parameterization with each beta_mu independent
  for(d in 1:D) {
    for(t in 1:T) {
      beta[t, d] <- beta_mu[d] + beta_sigma[d] * beta_std[t, d];
    }
  }
  for(d in 1:U) {
    for(t in 1:T) {
      alpha[t, d] <- alpha_mu[d] + alpha_sigma[d] * alpha_std[t, d];
    }
  }
  
  // assemble predictor w/ intercept + effects of indiv-level covariates
  for(t in 1:T) {
    for(n in 1:N) {
      // non-centered parameterization following Betacourt and Girolami
      pred[n, t] <- inv_logit(intercept[t] + (beta[t, ] * x[n]'));
    }
  }
}
model {
  intercept_mu ~ normal(-2, 2);
  sigma ~ cauchy(0, 1);
  
  for(t in 1:T) {  // intercept is unique for each time unit
    intercept[t] ~ normal(intercept_mu + eff_phase[phase[t]] + 
        alpha[t, ] * u[t]', sigma);
  }
  
  // change to multivariate normal prior?
  to_vector(beta_std) ~ normal(0, 1);
  beta_sigma ~ cauchy(0, 1);

  // diet
  //  carn (intercept)
  //  herb
  beta_mu[1] ~ normal(0.5, 1);
  //  inse
  beta_mu[2] ~ normal(0.5, 1);
  //  omni
  beta_mu[3] ~ normal(0.5, 1);
  
  // loco
  //  arboreal (intercept)
  //  digi
  beta_mu[4] ~ normal(0, 1);
  //  foss
  beta_mu[5] ~ normal(-0.5, 1);
  //  plan
  beta_mu[6] ~ normal(0, 1);
  //  scan
  beta_mu[7] ~ normal(0.5, 1);
  //  ungu
  beta_mu[8] ~ normal(0, 1);

  //mass
  beta_mu[9] ~ normal(0, 1);


  // change to multivariate normal prior?
  to_vector(alpha_std) ~ normal(0, 1);
  alpha_mu ~ normal(0, 1);
  alpha_sigma ~ cauchy(0, 1);

  eff_phase ~ normal(0, scale_phase);
  scale_phase ~ cauchy(0, 1);

  for(n in 1:N) {
    sight[n, ] ~ bernoulli(pred[n, ]);
  }
}
generated quantities {
  int sight_tilde[N, T];  // observed presence
  
  for(n in 1:N) {
    for(t in 1:T) {
      sight_tilde[n, t] <- bernoulli_rng(pred[n, t]);
    }
  }
}
