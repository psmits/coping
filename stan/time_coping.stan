data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1> D;
  real trait[N];
  int year[N];
  int diet[N];
}
parameters {
  // break this up; have diet be varying intercept with mean 0 WITH overall mean.
  real loc[T];  // location for year
  real loc_mu;
  real<lower=0> loc_scale;

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
  real<lower=0> scale[T];  // scale for year

  for(t in 1:T) {
    scale[t] <- exp(exp_scale[t]);
  }
}
model {
  for(d in 1:D) {
    diet_eff[d] ~ normal(0, diet_scale);
  }
  diet_scale ~ cauchy(0, 1);

  for(t in 1:T) {
    loc[t] ~ normal(loc_mu, loc_scale);
    exp_scale[t] ~ normal(scale_mu, scale_scale);
    skew[t] ~ normal(skew_mu, skew_scale);
  }
  loc_mu ~ normal(0, 10);
  loc_scale ~ cauchy(0, 1);

  scale_mu ~ normal(0, 5);
  scale_scale ~ cauchy(0, 1);
  
  skew_mu ~ normal(5, 10);
  skew_scale ~ cauchy(0, 1);


 // sampling statement/likelihood
  for(n in 1:N) {
    //if(year[n] == 1) {
    //  trait[n] ~ skew_normal(start,
    //      scale[year[n]], skew[year[n]]);
    //} else {
      trait[n] ~ skew_normal(loc[year[n]] + diet_eff[diet[n]],
          scale[year[n]], skew[year[n]]);
    //}
  }
}
