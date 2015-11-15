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

  vcv_inv <- inverse(vcv);
  vcv_inv_la <- transpose(cholesky_decompose(vcv_inv));
}
parameters {
  matrix[K, D] beta[C];  // matrix of regression coefficients
  matrix[K, D] beta_mu;
  matrix<lower=0>[K, D] sigma;
  vector[U] phy[K];  // need for each of K categories
  real<lower=0> sigma_phy[K];  // needed for each of K categories 
}
transformed parameters {
  real<lower=0> sig_phy_sq[K];  // needed for each of K categories 

  for(k in 1:K) {
    sig_phy_sq[k] <- sigma_phy[k]^2;  // needed for each of K categories 
  }
}
model {
  vector[U] v[K];  // needed for each of K categories 
  real sum_of_squares[K];  // needed for each of K categories 
  vector[K] hold[N];

  for(k in 1:K) {
    beta_mu[k] ~ normal(0, 1);
    sigma[k] ~ cauchy(0, 1);
    for(c in 1:C) {
      beta[c][k] ~ normal(beta_mu[k], sigma[k]);
    }
  }

  // needed for each of K categories 
  sigma_phy ~ cauchy(0,1);
  for(k in 1:K) {
    v[k] <- vcv_inv_la * phy[k];
    sum_of_squares[k] <- dot_product(v[k], v[k]);

    // non-constant part of log(det(sigma_phy * vcv) ^ -0.5
    increment_log_prob(-0.5 * N * log(sig_phy_sq[k]));
    // log of kernal of mulinorm
    increment_log_prob(sum_of_squares[k] / (2 * sig_phy_sq[k]));
  }
  for(n in 1:N) {
    for(k in 1:K) {
      hold[n][k] <- beta[cohort[n]][k] * x[n] + phy[k][id[n]];
    }
  }

  for(n in 1:N) {
    y[n] ~ categorical_logit(hold[n]);
  }
}
generated quantities {
  int y_tilde[N];
  vector[K] hold[N];
  
  for(n in 1:N) {
    for(k in 1:K) {
      hold[n][k] <- beta[cohort[n]][k] * x[n] + phy[k][id[n]];
    }
  }

  for(n in 1:N) {
    y_tilde[n] <- categorical_rng(softmax(hold[n]));
  }
}
