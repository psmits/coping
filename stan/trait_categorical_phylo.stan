data {
  int<lower=2> K;  // possible categories
  int<lower=0> N;  // sample size
  int<lower=1> D;  // number of predictors
  int<lower=1> C;
  int U;  // number of unique species
  int<lower=1,upper=K> y[N];  // state
  vector[D] x[N];  // matrix of predictors
  int cohort[N];
  int id[N];  // vector of species ids
  matrix[U, U] vcv;
}
transformed data {
  matrix[U, U] vcv_inv;
  matrix[U, U] vcv_inv_la;
  row_vector[D] zeroes;  // for mixing params and constants

  vcv_inv <- inverse(vcv);
  vcv_inv_la <- transpose(cholesky_decompose(vcv_inv));
  
  zeroes <- rep_row_vector(0, D);
}
parameters {
  matrix[K - 1, D] beta_raw[C];  // identifiable? sum to one constraint...
  matrix[K - 1, D] beta_mu;
  matrix<lower=0>[K - 1, D] sigma;
  vector[U] phy;
  real<lower=0> sigma_phy;
}
transformed parameters {
  matrix[K, D] beta[C];
  real<lower=0> sig_phy_sq;

  for(c in 1:C) { // mixing constants and params
    beta[c] <- append_row(beta_raw[c], zeroes);
  }
  sig_phy_sq <- sigma_phy^2;
}
model {
  vector[U] v; 
  real sum_of_squares;
  vector[K] hold[N];

  for(k in 1:(K - 1)) {
    beta_mu[k] ~ normal(0, 1);
    sigma[k] ~ cauchy(0, 1);
    for(c in 1:C) {
      beta_raw[c][k] ~ normal(beta_mu[k], sigma[k]);
    }
  }

  sigma_phy ~ cauchy(0,1);
  v <- vcv_inv_la * phy;
  sum_of_squares <- dot_product(v, v);

  // non-constant part of log(det(sigma_phy * vcv) ^ -0.5
  increment_log_prob(-0.5 * N * log(sig_phy_sq));
  // log of kernal of mulinorm
  increment_log_prob(sum_of_squares / (2 * sig_phy_sq));

  // assemble the length K vector of predictors for each n
  for(n in 1:N) {
    for(k in 1:K) {
      if(k != K) {
        hold[n][k] <- beta[cohort[n]][k] * x[n] + phy[id[n]];
      } else if(k == K) {
        hold[n][k] <- beta[cohort[n]][k] * x[n];
      }
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
      if(k != K) {
        hold[n][k] <- beta[cohort[n]][k] * x[n] + phy[id[n]];
      } else if(k == K) {
        hold[n][k] <- beta[cohort[n]][k] * x[n];
      }
    }
  }

  // generate posterior predictive datasets
  // softmax function because categorical rng has no logit form
  for(n in 1:N) {
    y_tilde[n] <- categorical_rng(softmax(hold[n]));
  }
}
