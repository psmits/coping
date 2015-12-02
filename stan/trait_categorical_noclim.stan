data {
  int<lower=2> K;  // possible categories
  int<lower=0> N;  // sample size
  int<lower=1> D;  // number of indiv-level predictors
  int<lower=1> C;
  int<lower=1,upper=K> y[N];  // state
  vector[D] x[N];  // matrix of indiv-level predictors
  int cohort[N];
}
transformed data {
  row_vector[D] zeroes;  // for mixing params and constants

  zeroes <- rep_row_vector(0, D);
}
parameters {
  matrix[K - 1, D] beta_raw[C];  // makes identifiable
  matrix[K - 1, D] beta_mu;
  matrix<lower=0>[K - 1, D] sigma;
}
transformed parameters {
  matrix[K, D] beta[C];  // makes identifiable

  for(c in 1:C) { // mixing constants and params
    beta[c] <- append_row(beta_raw[c], zeroes);
  }
}
model {
  vector[K] hold[N];

  // priors for reg coefs 
  // for each response, vary by geo-unit
  //  also group level predictor
  for(k in 1:(K - 1)) {
    beta_mu[k] ~ normal(0, 1);  
    // can make this MVN, but would need for all response types
    sigma[k] ~ cauchy(0, 1);
    for(c in 1:C) {
      beta_raw[c][k] ~ normal(beta_mu[k], sigma[k]);
    }
  }

  // assemble the length K vector of individual-level predictors for each n
  for(n in 1:N) {
    for(k in 1:K) {
      hold[n][k] <- beta[cohort[n]][k] * x[n];
    }
  }

  // final sampleing statement
  for(n in 1:N) {
    y[n] ~ categorical_logit(hold[n]);
  }
}
generated quantities {
  int y_tilde[N];
  vector[K] hold[N];

  for(n in 1:N) {
    for(k in 1:K) {
      hold[n][k] <- beta[cohort[n]][k] * x[n];
    }
  }

  // generate posterior predictive datasets
  // softmax function because categorical rng has no logit form
  for(n in 1:N) {
    y_tilde[n] <- categorical_rng(softmax(hold[n]));
  }
}
