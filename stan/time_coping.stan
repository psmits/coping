data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1> D;
  real trait[N];
  int year[N];
  int diet[N];
}
parameters {
  real start;
  real intercept;
  real alpha;
  real beta;

  vector[D] diet_eff;
  real<lower=0> diet_scale;

  real exp_scale[T];
  real scale_mu; 
  real<lower=0> scale_scale;

  real skew[T];  // skew for year
  real skew_mu;
  real<lower=0> skew_scale;
}
transformed parameters {
  real loc[T];  // location for year
  real<lower=0> scale[T];  // scale for year

  loc[1] <- start;
  for(t in 2:T) {
    loc[t] <- alpha * loc[t - 1] + beta * t + intercept;
  }

  for(t in 1:T) {
    scale[t] <- exp(exp_scale[t]);
  }
}
model {
  start ~ normal(0, 5);
  alpha ~ normal(0, 1);
  intercept ~ normal(0, 5);
  beta ~ normal(0, 1);

  diet_eff ~ normal(0, diet_scale);
  diet_scale ~ cauchy(0, 1);

  exp_scale ~ normal(scale_mu, scale_scale);
  scale_mu ~ normal(0, 5);
  scale_scale ~ cauchy(0, 1);

  skew ~ normal(skew_mu, skew_scale);
  skew_mu ~ normal(10, 5);
  skew_scale ~ cauchy(0, 1);

  for(n in 1:N) {
    trait[n] ~ skew_normal(loc[year[n]] + diet_eff[diet[n]], 
        scale[year[n]], skew[year[n]]);
  }
}
