data {
  int<lower=2> K;  // possible categories
  int<lower=0> N;  // sample size
  int<lower=1> D;  // number of indiv-level predictors
  int<lower=1> C;
  int<lower=1,upper=K> y[N];  // state
  vector[D] x[N];  // matrix of indiv-level predictors
  int cohort[N];
  real isoval[C];
  real isorang[C];
}
transformed data {
  row_vector[D] zeroes;  // for mixing params and constants
  vector[1] z;

  zeroes <- rep_row_vector(0, D);
  z[1] <- 0.0;
}
parameters {
  vector[K - 1] intercept_raw[C];
  vector[K - 1] intercept_mu;
  matrix<lower=0>[K - 1, D] sigma;
  
  matrix[K - 1, D] beta_raw;  // makes identifiable
  
  vector[K - 1] alpha;  // coef for group-level predictor
  vector[K - 1] gamma;  // coef for group-level predictor
}
transformed parameters {
  matrix[K, D] beta;  // makes identifiable
  vector[K] intercept[C];

  beta <- append_row(beta_raw, zeroes);
  for(c in 1:C) {
    intercept[c] <- append_row(intercept_raw[c], z);
  }
}
model {
  vector[K] hold[N];

  // priors for reg coefs 
  // for each response k, vary by time c
  for(k in 1:(K - 1)) {
    alpha[k] ~ normal(0, 1); 
    gamma[k] ~ normal(0, 1); 
    intercept_mu[k] ~ normal(0, 5);  
    sigma[k] ~ cauchy(0, 1);
    for(c in 1:C) {
      intercept_raw[c][k] ~ normal(intercept_mu[k] + alpha[k] * isoval[c] + 
          gamma[k] * isorang[c], sigma[k]);
      // only include group level predictor for intercept parameter
      // because it effects the baseline occurrence of response
      // assumes covariate effects aren't affected by "climate"
      // effect at time c is drawn from shared mean beta_mu
    }
  }

  to_vector(beta_raw) ~ normal(0, 1);

  // assemble the length K vector of effects
  for(n in 1:N) {
    for(k in 1:K) {
      hold[n][k] <- intercept[cohort[n]][k] + beta[k] * x[n];
    }
  }

  // final sampling statement
  for(n in 1:N) {
    y[n] ~ categorical_logit(hold[n]);
  }
}
generated quantities {
  int y_tilde[N];
  vector[K] hold[N];

  for(n in 1:N) {
    for(k in 1:K) {
      hold[n][k] <- intercept[cohort[n]][k] + beta[k] * x[n];
    }
  }

  // generate posterior predictive datasets
  // softmax function because categorical rng has no logit form
  for(n in 1:N) {
    y_tilde[n] <- categorical_rng(softmax(hold[n]));
  }
}
