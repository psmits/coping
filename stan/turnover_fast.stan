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
  real state_space_log(int[] y, row_vector pred, vector p) {
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

        sl <- bernoulli_log(z[1], pred[1]);
        for(j in 2:S) {
          sl <- sl + bernoulli_log(z[j], (z[j - 1] * pred[j]) + 
              ((1 - z[j - 1]) * pred[j]));
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
  vector[D] x[N];  // matrix of indiv-level covariates
  vector[U] u[T];  // matrix of group-level covariates

  int phase[T];  // plant phase
}
parameters {
  vector[T] inter_std;  // intercept varies by what point in time
  real intercept_mu;  // mean intercept
  real<lower=0> sigma;  // variance of varying-intercept

  vector[D] beta_mu;

  vector[T] p_norm; 
  real p_mu;
  real<lower=0> p_sigma;
}
transformed parameters {
  vector<lower=0,upper=1>[T] p;  // sampling probability
  matrix[N, T] pred; 

  // noncentered parameterization 
  for(t in 1:T) {
    p[t] <- inv_logit(p_mu + p_sigma * p_norm[t]);
  }

  // noncentered parameterization with each beta_mu independent
  // assemble predictor w/ intercept + effects of indiv-level covariates
  for(t in 1:T) {
    for(n in 1:N) {
      // non-centered parameterization following Betacourt and Girolami
      pred[n, t] <- inv_logit(intercept_mu + sigma*inter_std[t] + 
          (beta_mu' * x[n]));
    }
  }
}
model {
  p_mu ~ normal(0, 1);
  p_sigma ~ cauchy(0, 1);

  // diet
  //  carn (intercept)
  //  herb
  beta_mu[1] ~ normal(0.5, 1);
  //  inse
  beta_mu[2] ~ normal(0.5, 1);
  //  omni
  beta_mu[3] ~ normal(0.5, 1);
  
  // loco
  //  arboreal (intercept)
  //  digi
  beta_mu[4] ~ normal(0, 1);
  //  foss
  beta_mu[5] ~ normal(-0.5, 1);
  //  plan
  beta_mu[6] ~ normal(0, 1);
  //  scan
  beta_mu[7] ~ normal(0.5, 1);
  //  ungu
  beta_mu[8] ~ normal(0, 1);

  //mass
  beta_mu[9] ~ normal(0, 1);
 

  inter_std ~ normal(0, 1);
  intercept_mu ~ normal(-2, 2);
  sigma ~ cauchy(0, 1);

  for(n in 1:N) {
    sight[n] ~ state_space(pred[n, ], p);
  }
}
