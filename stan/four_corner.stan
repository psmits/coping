data {
  int N;
  int D;
  int U;
  int C;
  int cohort[N];
  int y[N];
  matrix[N, D] x;
  matrix[N, U] u;
}
parameters {
  real intercept[C];
  real intercept_mu;
  real<lower=0> sigma;
  vector[D] beta;
  vector[U] gamma;
  real<lower=0> lambda;
  real<lower=0> phi_d[D];
  real<lower=0> phi_u[U];
}
model {
  // final sampling statement

  for(c in 1:C) {
    intercept[c] ~ normal(intercept_mu, sigma);
  }
  for(d in 1:D) {
    beta[d] ~ normal(0, lambda * phi_d[d]);
  }
  for(uu in 1:U) {
    gamma[uu] ~ normal(0, lambda * phi_u[uu]);
  }

  intercept_mu ~ normal(0, 5);
  sigma ~ cauchy(0, 1);
  lambda ~ cauchy(0, 1);
  phi_d ~ cauchy(0, 1);
  phi_u ~ cauchy(0, 1);

  for(n in 1:N) {
    y[n] ~ bernoulli_logit(intercept[cohort[n]] + x[n] * beta + u[n] * gamma);
  }
}
