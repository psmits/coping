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

  zeroes <- rep_row_vector(0, D);
}
parameters {
  matrix[K - 1, D] beta_raw[C];  // makes identifiable
  matrix[K - 1, D] beta_mu;
  matrix<lower=0>[K - 1, D] sigma;
  vector[K - 1] alpha;  // coef for group-level predictor
  vector[K - 1] gamma;  // coef for group-level predictor
}
transformed parameters {
  matrix[K, D] beta[C];  // makes identifiable

  for(c in 1:C) { // mixing constants and params
    beta[c] <- append_row(beta_raw[c], zeroes);
  }
}
model {
  vector[K] hold[N];

  // ideas on how to rewrite model
  // intercept is own term
  // vector[K - 1] intercept[C];
  //   varies with time
  // covariate effects constant with time
  // matrix[K - 1, D] beta_raw;
  // matrix[K, D] beta;
  //  x no longer includes column of 1s
  //  D no longer counts intercept

  // priors for reg coefs 
  // for each response k, vary by time c
  for(k in 1:(K - 1)) {
    alpha[k] ~ normal(0, 1); 
    gamma[k] ~ normal(0, 1); 
    beta_mu[k] ~ normal(0, 1);  
    // can make this MVN, but would need for all response types
    sigma[k] ~ cauchy(0, 1);
    for(d in 1:D) {
      for(c in 1:C) {
        if(d == 1) {
          beta_raw[c][k][d] ~ normal(beta_mu[k] + alpha[k] * isoval[c] + 
              gamma[k] * isorang[c], 
              sigma[k]);
          // only include group level predictor for intercept parameter
          // because it effects the baseline occurrence of response
          // assumes covariate effects aren't affected by "climate"
        } else {
          beta_raw[c][k][d] ~ normal(beta_mu[k], sigma[k]);
          // effect at time c is drawn from shared mean beta_mu
        }
      }
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
//generated quantities {
//  int y_tilde[N];
//  vector[K] hold[N];
//
//  for(n in 1:N) {
//    for(k in 1:K) {
//      hold[n][k] <- beta[cohort[n]][k] * x[n];
//    }
//  }
//
//  // generate posterior predictive datasets
//  // softmax function because categorical rng has no logit form
//  for(n in 1:N) {
//    y_tilde[n] <- categorical_rng(softmax(hold[n]));
//  }
//}
