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
  vector[U] phy;
  real<lower=0> sigma_phy;
}
transformed parameters {
  real<lower=0> sig_phy_sq;

  sig_phy_sq <- sigma_phy^2;
}
model {
  vector[U] v;
  real sum_of_squares;

  for(k in 1:K) {
    beta_mu[k] ~ normal(0, 1);
    sigma[k] ~ cauchy(0, 1);
    for(c in 1:C) {
      beta[c][k] ~ normal(beta_mu[k], sigma[k]);
    }
  }

  sigma_phy ~ cauchy(0,1);
  v <- vcv_inv_la * phy;
  sum_of_squares <- dot_product(v, v);

  // non-constant part of log(det(sigma_phy * vcv) ^ -0.5
  increment_log_prob(-0.5 * N * log(sig_phy_sq));
  // log of kernal of mulinorm
  increment_log_prob(sum_of_squares / (2 * sig_phy_sq));

  for(n in 1:N) {
    y[n] ~ categorical_logit(beta[cohort[n]] * x[n] + phy[id[n]]);
  }
}
generated quantities {
  int y_tilde[N];

  for(n in 1:N) {
    y_tilde[n] <- categorical_rng(softmax(beta[cohort[n]] * x[n] + phy[id[n]]));
  }
}
