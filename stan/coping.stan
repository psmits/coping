data {
  int<lower=1> N;
  int<lower=1> T;
  real trait[N];
  int year[N];
}
transformed data {
}
parameters {
  real loc[T];  // location for year
  real loc_mu;
  real<lower=0> loc_scale;

  real exp_scale[T];
  real scale_mu; 
  real<lower=0> scale_scale;

  real skew[T];  // skew for year
  real skew_mu;
  real<lower=0> skew_scale;

  real beta[3];
  real<lower=0> sigma;
}
transformed parameters {
  real expect[T];
  real<lower=0> scale[T];  // scale for year

  for(t in 1:T) {
    scale[t] <- exp(exp_scale[t]);
  }

  for(t in 1:T) {
    expect[t] <- loc[t] + scale[t] * (skew[t] / sqrt(1 + skew[t])) * sqrt(2 / pi());
  }
}
model {
  for(t in 1:T) {
    loc[t] ~ normal(loc_mu, loc_scale);
    exp_scale[t] ~ normal(scale_mu, scale_scale);
    skew[t] ~ normal(skew_mu, skew_scale);
  }
  loc_mu ~ normal(0, 5);
  loc_scale ~ cauchy(0, 1);
  
  scale_mu ~ normal(0, 5);
  scale_scale ~ cauchy(0, 1);
  
  skew_mu ~ normal(5, 10);
  skew_scale ~ cauchy(0, 1);

  for(t in 2:T) {
    increment_log_prob(normal_log(expect[t], 
      beta[1] + beta[2] * expect[t - 1] + beta[3] * t, 
      sigma));
  }
  beta[1] ~ normal(0, 5);
  beta[2] ~ normal(0, 1);
  beta[3] ~ normal(0, 1);
  sigma ~ cauchy(0, 1);
  
 // sampling statement/likelihood
  for(n in 1:N) {
    trait[n] ~ skew_normal(loc[year[n]], scale[year[n]], skew[year[n]]);
  }
}
