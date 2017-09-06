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
  real state_space_lp(int[] y, row_vector origin, row_vector stay, row_vector p) {
    int ft;
    int lt;
    int S;
    int i;
    int prod_term;
    int lp_size;

    S = size(y);
    i = 1;
    ft = first_capture(y);
    lt = last_capture(y);
    lp_size = ft * (S - lt + 1);

    // have to go through all possible event histories for each species
    {
      vector[lp_size] lp;
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
          sl = bernoulli_lpmf(z[1] | origin[1]);
          prod_term = 1 - z[1];

          // z as function of trait state
          {
            vector[S - 1] gg;
            for(j in 2:S) {
              prod_term = prod_term * (1 - z[j - 1]);
              gg[j - 1] = bernoulli_lpmf(z[j] | (z[j - 1] * stay[j - 1]) + 
                  prod_term * origin[j]);
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
}
data {
  int N;  // sample size of taxa
  int T;  // sample size of temporal units
  int D;  // number of indiv-level predictors
  int U;  // number of group-level predictors

  int sight[N, T];  // observed presence
  int state[N];  // vector of ecology state

  vector[N] mass;

  matrix[T, U] ufull;  // matrix of group-level covariates
  matrix[T - 1, U] u;  // matrix of group-level covariates
}
parameters {
  // origination
  real o_b_1; // effect of mass on occurrence
  real o_b_2; // effect of mass on occurrence

  // effect associated with ecology
  matrix[D, T] o_a_z; // part of non-centering
  cholesky_factor_corr[D] o_L_Omega;
  vector<lower=0>[D] o_tau;

  matrix[U, D] o_gamma; // effect of group level covariates

  // survival
  real s_b_1; // effect of mass on occurrence
  real s_b_2; // effect of mass on occurrence

  // effect associated with ecology
  matrix[D, T - 1] s_a_z; // part of non-centering
  cholesky_factor_corr[D] s_L_Omega;
  vector<lower=0>[D] s_tau;

  matrix[U, D] s_gamma; // effect of group level covariates

  // preservation
  real p_b_1;
  real p_b_2;

  // effect associated with ecology
  matrix[D, T] p_a_z; // part of non-centering
  cholesky_factor_corr[D] p_L_Omega;
  vector<lower=0>[D] p_tau;

  matrix[U, D] p_gamma; // effect of group level covariates
}
transformed parameters {
  matrix[T, D] o_a;  // origin: effect associated with ecology
  matrix[T - 1, D] s_a;  // survival: effect associated with ecology
  matrix[T, D] p_a;  // origin: effect associated with ecology
  matrix[N, T] origin;  // occurrence probabilty 
  matrix[N, T - 1] stay;  // occurrence probabilty 
  matrix[N, T] p;  // sampling probability

  // vectorized, non-centered, chol-decom group-level predictors
  o_a = ufull * o_gamma + (diag_pre_multiply(o_tau, o_L_Omega) * o_a_z)'; // or
  s_a = u * s_gamma + (diag_pre_multiply(s_tau, s_L_Omega) * s_a_z)'; // ext
  p_a = ufull * p_gamma + (diag_pre_multiply(p_tau, p_L_Omega) * p_a_z)'; // pr

  // probability of occurring

  for(t in 1:T) {
    for(n in 1:N) {
      origin[n, t] = inv_logit(o_a[t, state[n]] + 
          o_b_1 * mass[n] + o_b_2 * (mass[n] ^ 2));
      p[n, t] = inv_logit(p_a[t, state[n]] + 
          p_b_1 * mass[n] + p_b_2 * (mass[n] ^ 2));
    }
  }
  for(t in 1:(T-1)) {
    for(n in 1:N) {
      stay[n, t] = inv_logit(s_a[t, state[n]] + 
          s_b_1 * mass[n] + s_b_2 * (mass[n] ^ 2));
    }
  }
}
model {
  // effects of ecologies
  // origin
  to_vector(o_a_z) ~ normal(0, 1);
  o_L_Omega ~ lkj_corr_cholesky(2);
  o_tau ~ normal(0, 1);
  to_vector(o_gamma) ~ normal(0, 1);

  // survival
  to_vector(s_a_z) ~ normal(0, 1);
  s_L_Omega ~ lkj_corr_cholesky(2);
  s_tau ~ normal(0, 1);
  to_vector(s_gamma) ~ normal(0, 1);
  
  // survival
  to_vector(p_a_z) ~ normal(0, 1);
  p_L_Omega ~ lkj_corr_cholesky(2);
  p_tau ~ normal(0, 1);
  to_vector(p_gamma) ~ normal(0, 1);

  // effects of mass on orignation, survival, observation
  o_b_1 ~ normal(0, 1);
  o_b_2 ~ normal(0, 1);
  s_b_1 ~ normal(0, 1);
  s_b_2 ~ normal(0, 1);
  p_b_1 ~ normal(0, 1);
  p_b_2 ~ normal(0, 1);

  for(n in 1:N) {
    target += state_space_lp(sight[n], origin[n, ], stay[n, ], p[n, ]);
  }
}
generated quantities {
  corr_matrix[D] o_Omega;
  corr_matrix[D] s_Omega;
  corr_matrix[D] p_Omega;

  // convert back to a correlation matrix
  o_Omega = multiply_lower_tri_self_transpose(o_L_Omega);
  s_Omega = multiply_lower_tri_self_transpose(s_L_Omega);
  p_Omega = multiply_lower_tri_self_transpose(p_L_Omega);
}
