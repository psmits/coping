functions {
  int first_capture(int[] y_i) {
    for (k in 1:size(y_i))
      if (y_i[k])
        return k;
    return 0;
  }
  int last_capture(int[] y_i) {
    for (k_rev in 0:(size(y_i) - 1)) {
      int k;
      k <- size(y_i) - k_rev;
      if (y_i[k])
        return k;
    }
    return 0;
  }
  real state_space_log(int[] y, vector intercept, real h, vector p) {
    int ft;
    int lt;
    int S;
    vector[first_capture(y) * (size(y) - last_capture(y) + 1)] lp;
    int i;

    ft <- first_capture(y);
    lt <- last_capture(y);
    S <- size(y);
    i <- 1;

    for(t_first_alive in 1:ft) {
      for (t_last_alive in lt:S) {
        real sl;
        int z[S];

        for(l in 1:S) {
          z[l] <- 0;
        }
        for(a in t_first_alive:t_last_alive) {
          z[a] <- 1;
        }

        for(j in 1:S) {
          sl <- sl + bernoulli_log(z[j], intercept[j] + h);
        }
        for(k in 1:S) {
          sl <- sl + bernoulli_log(y[k], z[k] * p[k]);
        }

        lp[i] <- sl;
        i <- i + 1;
      }
    }
    return log_sum_exp(lp);
  }
}
data {
  int N;  // sample size of taxa
  int T;  // sample size of temporal units
  int D;  // number of indiv-level predictors
  int U;  // number of group-level predictors
  int P;  // number of floral phases
  
  int sight[N, T];  // observed presence
  row_vector[D] x[N];  // matrix of indiv-level covariates
  row_vector[U] u[T];  // matrix of group-level covariates

  int phase[T];  // plant phase
}
parameters {
  vector[T] intercept;  // intercept varies by what point in time
  real intercept_mu;  // mean intercept
  real<lower=0> sigma;  // variance of varying-intercept
  
  vector[P] eff_phase;  // effect of plant-phase (mean 0)
  real<lower=0> scale_phase;  // variance in plant-phase effect

  vector[D] beta;  // effect of indiv-level covariates
  vector[U] alpha;  // effect of group-level covariates

  vector[T] p_norm; 
  real p_mu;
  real<lower=0> p_sigma;
}
transformed parameters {
  vector[N] h;
  vector<lower=0,upper=1>[T] p;  // presence

  // compose the indidual effects
  //   will want to extend to break point
  //   beta in two parts
  //   this is on the complicated side of things because of marginalization
  for(n in 1:N) {
    h[n] <- x[n] * beta;
  }

  for(t in 1:T) {
    p[t] <- inv_logit(p_norm[t]);
  }
}
model {
  for(t in 1:T) {
    p_norm[t] ~ normal(p_mu, p_sigma);
  }
  p_mu ~ normal(0, 1);
  p_sigma ~ cauchy(0, 1);

  for(t in 1:T) {  // intercept is unique for each time unit
    intercept[t] ~ normal(intercept_mu + eff_phase[phase[t]] + 
        u[t] * alpha, sigma);
  }
  
  intercept_mu ~ normal(0, 5);
  sigma ~ cauchy(0, 1);
  eff_phase ~ normal(0, scale_phase);
  scale_phase ~ cauchy(0, 1);
  
  beta ~ normal(0, 1);
  alpha ~ normal(0, 1);
  
  for(n in 1:N) {
    sight[n] ~ state_space(intercept, h[n], p);
  }
}
generated quantities {
  matrix[N, T] y_tilde;  // generate a "true" record

  for(n in 1:N) {
    for(t in 1:T) {
      y_tilde[n, t] <- bernoulli_rng(intercept[t] + h[n]);
    }
  }
}
