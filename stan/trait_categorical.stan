data {
  int<lower=2> K;  // possible categories
  int<lower=0> N;  // sample size
  int<lower=1> D;  // number of predictors
  int<lower=1> C;
  int<lower=1,upper=K> y[N];  // state
  vector[D] x[N];  // matrix of predictors
  int cohort[N];
}
parameters {
  matrix[K, D] beta[C];  // matrix of regression coefficients
  matrix[K, D] beta_mu;
  matrix<lower=0>[K, D] sigma;
}
model {
  for(k in 1:K) {
    beta_mu[k] ~ normal(0, 1);
    sigma[k] ~ cauchy(0, 1);
    for(c in 1:C) {
      beta[c][k] ~ normal(beta_mu[k], sigma[k]);
    }
  }
  for(n in 1:N) {
    y[n] ~ categorical(softmax(beta[cohort[n]] * x[n]));
  }
}
generated quantities {
  int y_tilde[N];

  for(n in 1:N) {
    y_tilde[n] <- categorical_rng(softmax(beta[cohort[n]] * x[n]));
  }
}
