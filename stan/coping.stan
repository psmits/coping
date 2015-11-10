data {
  int<lower=1> N;
  int<lower=1> T;
  int D;
  real trait[N];
  int year[N];
  matrix[N, D - 1] cat;
}
parameters {
  // break this up; have diet be varying intercept with mean 0 WITH overall mean.
  real loc[T];  // location for year
  real loc_mu;
  real<lower=0> loc_scale;
  // real mu;
  // real diet_eff[D];
  // real<lower=0> diet_var;

  real exp_scale[T];
  real scale_mu; 
  real<lower=0> scale_scale;

  real skew[T];  // skew for year
  real skew_mu;
  real<lower=0> skew_scale;

  vector[D - 1] beta_cat;
}
transformed parameters {
  real<lower=0> scale[T];  // scale for year

  for(t in 1:T) {
    scale[t] <- exp(exp_scale[t]);
  }
}
model {
  for(t in 1:T - 1) {
    loc[t] ~ normal(loc_mu, loc_scale);
    exp_scale[t] ~ normal(scale_mu, scale_scale);
    skew[t] ~ normal(skew_mu, skew_scale);
  }
  
  scale_mu ~ normal(0, 5);
  scale_scale ~ cauchy(0, 1);
  
  skew_mu ~ normal(5, 10);
  skew_scale ~ cauchy(0, 1);

  beta_cat ~ normal(0, 1);

 // sampling statement/likelihood
  for(n in 1:N) {
    trait[n] ~ skew_normal(loc[year[n]] + cat[n] * beta_cat, 
        scale[year[n]], skew[year[n]]);
  }
}
