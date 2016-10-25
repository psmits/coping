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
      k = size(y_i) - k_rev;
      if (y_i[k])
        return k;
    }
    return 0;
  }
  real state_space_lp(int[] y, real phi, row_vector pred, row_vector p) {
    int ft;
    int lt;
    int S;
    vector[first_capture(y) * (size(y) - last_capture(y) + 1)] lp;
    int i;
    int prod_term;

    ft = first_capture(y);
    lt = last_capture(y);
    S = size(y);
    i = 1;

    for(t_first_alive in 1:ft) {
      for (t_last_alive in lt:S) {
        real sl;
        int z[S];

        for(l in 1:S) {
          z[l] = 0;
        }
        for(a in t_first_alive:t_last_alive) {
          z[a] = 1;
        }

        sl = bernoulli_lpmf(z[1] | pred[1]);
        prod_term = 1 - z[1];
        for(j in 2:S) {
          prod_term = prod_term * (1 - z[j - 1]);
          sl = sl + bernoulli_lpmf(z[j] | (z[j - 1] * pred[j - 1]) + 
              prod_term * pred[j - 1]);
        }
        for(k in 1:S) {
          sl = sl + bernoulli_lpmf(y[k] | z[k] * p[k]);
        }

        lp[i] = sl;
        i = i + 1;
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

  int sight[N, T];  // observed presence
  matrix[N, D] x;  // matrix of indiv-level covariates

  vector[N] mass;

  row_vector[U] u[T];  // matrix of group-level covariates
}
parameters {
  corr_matrix[D] Omega;
  vector<lower=0>[D] tau;
  vector[D] beta[T-1];
  
  matrix[U, D] gamma; 

  vector[T] alpha_0;
  real alpha_1;
  real mu;
  real<lower=0> sigma;

  real<lower=0,upper=1> phi;
}
transformed parameters {
  matrix<lower=0,upper=1>[N, T] p;  // sampling probability
  matrix[N, T - 1] pred; 
  
  for(n in 1:N) {
    for(t in 1:T) {
      p[n, t] = inv_logit(alpha_0[t] + alpha_1 * mass[n]);
    }
  }
  
  // assemble predictor w/ intercept + effects of indiv-level covariates
  for(t in 1:(T-1)) {
    for(n in 1:N) {
      pred[n, t] = inv_logit(x[n, ] * beta[t, ]);
    }
  }
}
model {
  tau ~ normal(0, 1);
  Omega ~ lkj_corr(2);

  {
    matrix[D, D] Sigma_beta;
    Sigma_beta = quad_form_diag(Omega, tau);

    for(t in 1:(T - 1)) {
      beta[t] ~ multi_normal(u[t + 1] * gamma, Sigma_beta);
    }
  }

  to_vector(gamma) ~ normal(0, 1);

  alpha_0 ~ normal(mu, sigma);
  mu ~ normal(0, 1);
  sigma ~ normal(0, 1);

  alpha_1 ~ normal(0, 1);

  for(n in 1:N) {
    target += state_space_lp(sight[n], phi, pred[n, ], p[n, ]);
  }
}
