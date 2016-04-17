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

  matrix[T, D] beta_std;
  real beta_mu[D];
  real<lower=0> beta_sigma[D];

  vector[T] p_norm; 
  real p_mu;
  real<lower=0> p_sigma;
}
transformed parameters {
  vector<lower=0,upper=1>[T] p;  // sampling probability
  matrix[N, T] pred; 
  matrix[T, D] beta;  // effect of indiv-level covariates

  // noncentered parameterization 
  for(t in 1:T) {
    p[t] <- inv_logit(p_mu + p_sigma * p_norm[t]);
  }

  // noncentered parameterization with each beta_mu independent
  for(d in 1:D) {
    for(t in 1:T) {
      beta[t, d] <- beta_mu[d] + beta_sigma[d] * beta_std[t, d];
    }
  }

  // assemble predictor w/ intercept + effects of indiv-level covariates
  for(t in 1:T) {
    for(n in 1:N) {
      // non-centered parameterization following Betacourt and Girolami
      pred[n, t] <- inv_logit(intercept_mu + sigma*inter_std[t] + 
          (beta[t, ] * x[n]));
    }
  }
}
model {
  p_mu ~ normal(0, 1);
  p_sigma ~ cauchy(0, 1);

  // change to multivariate normal prior?
  to_vector(beta_std) ~ normal(0, 1);
  beta_mu ~ normal(0, 1);
  beta_sigma ~ cauchy(0, 1);

  inter_std ~ normal(0, 1);
  intercept_mu ~ normal(0, 5);
  sigma ~ cauchy(0, 1);

  for(n in 1:N) {
    sight[n] ~ state_space(pred[n, ], p);
  }
}
generated quantities {
  matrix[T, N] z_tilde;
  matrix[T, N] y_tilde;
  real log_lik[N];

  for(n in 1:N) {
    z_tilde[1, n] <- bernoulli_rng(pred[n, 1]);
    y_tilde[1, n] <- bernoulli_rng(z_tilde[1, n] * p[1]);
    for(t in 2:T) {
      z_tilde[t, n] <- bernoulli_rng(z_tilde[t - 1, n] * pred[n, t] +
          ((prod(1 - z_tilde[1:(t - 1), n])) * pred[n, t]));
      y_tilde[t, n] <- bernoulli_rng(z_tilde[t, n] * p[t]);
    }
    log_lik[n] <- state_space_log(sight[n], pred[n, ], p);
  }
}
