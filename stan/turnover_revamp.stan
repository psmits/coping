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

    // have to go through all possible event histories for each species
    for(t_first_alive in 1:ft) {
      for (t_last_alive in lt:S) {
        real sl;
        int z[S];

        for(j in 1:S) {
          z[j] = 0;
        }
        for(a in t_first_alive:t_last_alive) {
          z[a] = 1;
        }

        // initial conditions
        sl = bernoulli_lpmf(z[1] | phi);
        prod_term = 1 - z[1];

        // z as function of trait state
        {
          vector[S - 1] gg;
          for(j in 2:S) {
            prod_term = prod_term * (1 - z[j - 1]);
            gg[j - 1] = bernoulli_lpmf(z[j] | (z[j - 1] * pred[j - 1]) + 
                prod_term * pred[j - 1]);
          }
          sl = sl + sum(gg);
        }

        // finally, y as function of z and p
        {
          vector[S] hh;
          for(k in 1:S) {
            hh[k] = bernoulli_lpmf(y[k] | z[k] * p[k]);
          }
          sl = sl + sum(hh);
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
  int state[N];  // vector of ecology state

  vector[N] mass;

  matrix[T - 1, U] u;  // matrix of group-level covariates
}
parameters {
  real b; // effect of mass on occurrence
  matrix[U, D] gamma; // effect of group level covariates
  
  matrix[D, T - 1] a_z; // part of non-centering
  cholesky_factor_corr[D] L_Omega;
  vector<lower=0>[D] tau;

  real alpha_0;
  real alpha_1;
  vector[T] alpha_time;
  real<lower=0> sigma;

  real<lower=0,upper=1> phi;
}
transformed parameters {
  matrix[T - 1, D] a;  // effect associated with ecology
  matrix<lower=0,upper=1>[N, T] p;  // sampling probability
  matrix[N, T - 1] pred;  // occurrence probabilty 
  
  // vectorized, non-centered, chol-decom group-level predictors
  // multivariate normal
  a = u * gamma + (diag_pre_multiply(tau, L_Omega) * a_z)';

  // probability of observing
  for(n in 1:N) {
    for(t in 1:T) {
      p[n, t] = inv_logit(alpha_0 + alpha_time[t] + alpha_1 * mass[n]);
    }
  }
  
  // probability of occurring
  for(t in 1:(T-1)) {
    for(n in 1:N) {
      pred[n, t] = inv_logit(a[t, state[n]] + b * mass[n]);
    }
  }
}
model {
  to_vector(a_z) ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(2);
  tau ~ normal(0, 1);
  to_vector(gamma) ~ normal(0, 1);

  alpha_0 ~ normal(0, 5);
  alpha_1 ~ normal(0, 1);
  alpha_time ~ normal(0, sigma);
  sigma ~ normal(0, 1);
  
  b ~ normal(0, 1);

  for(n in 1:N) {
    target += state_space_lp(sight[n], phi, pred[n, ], p[n, ]);
  }
}
generated quantities {
  corr_matrix[D] Omega;
  
  // convert back to a correlation matrix
  Omega = multiply_lower_tri_self_transpose(L_Omega);
}
