data {
  int<lower=2> K;  // possible categories
  int<lower=0> N;  // sample size
  int<lower=1> D;  // number of indiv-level predictors
  int<lower=1> U;  // number of group-level predictors
  int<lower=1> C;
  int<lower=1,upper=K> y[N];  // state
  vector[D] x[N];  // matrix of indiv-level predictors
  int cohort[N];
  vector[U] u[C];  // matrix of group-level predictors
}
transformed data {
  row_vector[D] zeroes;  // for mixing params and constants
  vector[1] z;
  real log_unif;

  zeroes <- rep_row_vector(0, D);
  z[1] <- 0.0;

  log_unif <- -log(C);
}
parameters {
  vector[K - 1] intercept_raw[C];
  vector[K - 1] intercept_mu_first;  // group-level intercept
  vector[K - 1] intercept_mu_second;  // group-level intercept
  vector<lower=0>[K - 1] sigma;
  
  matrix[K - 1, D] beta_raw;  // makes identifiable
  
  matrix[K - 1, U] alpha;  // effects of group-level covariates
}
transformed parameters {
  matrix[K, D] beta;  // makes identifiable
  vector[K] intercept[C];
  vector[C] lp;
  lp <- rep_vector(log_unif, C);

  beta <- append_row(beta_raw, zeroes);
  for(c in 1:C) {
    intercept[c] <- append_row(intercept_raw[c], z);
  }
  
  for(k in 1:(K - 1)) {
    for(c1 in 1:C) {
      for(c2 in 1:C) {
        if(c1 < c2) {
          lp[c1] <- lp[c1] + normal_log(intercept_raw[c1][k],
            intercept_mu_first[k] + alpha[k] * u[c1], sigma[k]);
        } else {
          lp[c1] <- lp[c1] + normal_log(intercept_raw[c1][k],
            intercept_mu_second[k] + alpha[k] * u[c1], sigma[k]);
        }
      }
    }
  }
}
model {
  vector[K] hold[N];

  intercept_mu_first ~ normal(0, 5);  
  intercept_mu_second ~ normal(0, 5);  
  sigma ~ cauchy(0, 1); 
  increment_log_prob(log_sum_exp(lp));

  to_vector(beta_raw) ~ normal(0, 1);
  to_vector(alpha) ~ normal(0, 1);

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
  int<lower=1,upper=C> s;

  s <- categorical_rng(softmax(lp));

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

