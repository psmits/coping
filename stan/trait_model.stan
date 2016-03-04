data {
  int<lower=2> K;  // possible categories
  int<lower=0> N;  // sample size
  int<lower=1> D;  // number of indiv-level predictors
  int<lower=1> U;  // number of group-level predictors
  int<lower=1> C;
  int<lower=1> P;
  int<lower=1,upper=K> y[N];  // state
  vector[D] x[N];  // matrix of indiv-level predictors
  int cohort[N];
  vector[U] u[C];  // matrix of group-level predictors
  int phase[P];
}
transformed data {
  row_vector[D] zeroes;  // for mixing params and constants
  vector[1] z;

  zeroes <- rep_row_vector(0, D);
  z[1] <- 0.0;
}
parameters {
  vector[K - 1] intercept_raw[C];
  vector[K - 1] intercept_mu;  // group-level intercept
  vector<lower=0>[K - 1] sigma;
  
  matrix[K - 1, D] beta_raw;  // makes identifiable
  
  matrix[K - 1, U] alpha;  // effects of group-level covariates
  vector[C] eff_phase;
  real<lower=0> scale_phase;
}
transformed parameters {
  matrix[K, D] beta;  // makes identifiable
  vector[K] intercept[C];
  vector[K] hold[N];

  beta <- append_row(beta_raw, zeroes);
  for(c in 1:C) {
    intercept[c] <- append_row(intercept_raw[c], z);
  }
  
  // assemble the length K vector of effects
  for(n in 1:N) {
    for(k in 1:K) {
      hold[n][k] <- intercept[cohort[n]][k] + beta[k] * x[n];
    }
  }
}
model {
  intercept_mu ~ normal(0, 5);  
  sigma ~ cauchy(0, 1); 
  eff_phase ~ normal(0, scale_phase);
  scale_phase ~ cauchy(0, 1);
  for(k in 1:(K - 1)) {
    for(c in 1:C) {
      intercept_raw[c][k] ~ normal(intercept_mu[k] + eff_phase[phase[c]] + 
          alpha[k] * u[c], sigma[k]);
    }
  }

  to_vector(beta_raw) ~ normal(0, 1);
  to_vector(alpha) ~ normal(0, 1);

  // final sampling statement
  for(n in 1:N) {
    y[n] ~ categorical_logit(hold[n]);
  }
}
generated quantities {
  int y_tilde[N];

  // generate posterior predictive datasets
  // softmax function because categorical rng has no logit form
  for(n in 1:N) {
    y_tilde[n] <- categorical_rng(softmax(hold[n]));
  }
}
