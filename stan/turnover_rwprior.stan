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
  int O;  // number of mammal orders

  int sight[N, T];  // observed presence
  int state[N];  // vector of ecology state

  vector[N] mass;  // each species mass

  int ords[N];  // each species order

  matrix[T, U] ufull;  // matrix of group-level covariates
  matrix[T - 1, U] u;  // matrix of group-level covariates
}
parameters {
  // origination
  real o_b_1; // effect of mass on occurrence
  // effect associated with ecology
  matrix[D, T] o_a_z; // part of non-centering
  cholesky_factor_corr[D] o_L_Omega;
  vector<lower=0>[D] o_tau;

  matrix[T, D] o_inter;
  matrix[U - 1, D] o_gamma; // effect of group level covariates

  //// effect of order
  vector[O] o_ordeff;
  real<lower=0> o_ordscale;

  ////
  // survival
  real s_b_1; // effect of mass on occurrence
  // effect associated with ecology
  matrix[D, T - 1] s_a_z; // part of non-centering
  cholesky_factor_corr[D] s_L_Omega;
  vector<lower=0>[D] s_tau;

  matrix[T - 1, D] s_inter;
  matrix[U - 1, D] s_gamma; // effect of group level covariates

  //// effect of order
  vector[O] s_ordeff;
  real<lower=0> s_ordscale;

  ////
  // preservation
  vector[T] p_timeeff; // time eff
  real p_b_1;  // mass coef
  vector[D] p_funceff; // ecology eff
  real<lower=0> p_funcscale;
}
transformed parameters {
  matrix[T, D] o_a;  // origin: effect associated with ecology
  matrix[T - 1, D] s_a;  // survival: effect associated with ecology
  matrix[N, T] origin;  // occurrence probabilty 
  matrix[N, T - 1] stay;  // occurrence probabilty 

  matrix[N, T] p;  // sampling probability

  // vectorized, non-centered, chol-decom group-level predictors
  o_a = (diag_pre_multiply(o_tau, o_L_Omega) * o_a_z)'; // or
  s_a = (diag_pre_multiply(s_tau, s_L_Omega) * s_a_z)'; // ext
  for(d in 1:D) {
    for(t in 1:T) {
      o_a[t, d] = o_a[t, d] + o_inter[t, d] + o_gamma[1, d] * ufull[t, 2] + 
        o_gamma[2, d] * ufull[t, 3] + o_gamma[3, d] * ufull[t, 4];
    }
    for(t in 1:(T - 1)) {
      s_a[t, d] = s_a[t, d] + s_inter[t, d] + s_gamma[1, d] * u[t, 2] + 
        s_gamma[2, d] * u[t, 3] + s_gamma[3, d] * u[t, 4];
    }
  }

  // probability of occurring and sampling
  for(t in 1:T) {
    for(n in 1:N) {
      origin[n, t] = inv_logit(o_a[t, state[n]] + o_b_1 * mass[n] + o_ordeff[ords[n]]);
      p[n, t] = inv_logit(p_timeeff[t] + p_b_1 * mass[n] + p_funceff[state[n]]);
    }
  }
  for(t in 1:(T-1)) {
    for(n in 1:N) {
      stay[n, t] = inv_logit(s_a[t, state[n]] + s_b_1 * mass[n] + s_ordeff[ords[n]]);
    }
  }
}
model {
  ////
  // origination
  to_vector(o_a_z) ~ normal(0, 1);
  o_L_Omega ~ lkj_corr_cholesky(2);  // really strong!
  o_tau ~ normal(0, 1);
  to_vector(o_gamma) ~ normal(0, 0.5);  // really strong for no eff!
  o_inter[1] ~ normal(0, 1);
  for(ii in 2:T) {
    o_inter[ii] ~ normal(o_inter[ii - 1], 1);
  }

  // priors for order effect on origination
  o_ordeff ~ normal(0, o_ordscale);
  o_ordscale ~ normal(0, 0.5);  // really strong!

  // priors for effects of mass on orignation
  o_b_1 ~ normal(0, 0.5);


  ////
  // survival 
  // by functional group at time w/ random walk prior
  to_vector(s_a_z) ~ normal(0, 1);
  s_L_Omega ~ lkj_corr_cholesky(2);  // really strong for no corr!
  s_tau ~ normal(0, 1);
  to_vector(s_gamma) ~ normal(0, 0.5);  // really strong for no eff!
  s_inter[1] ~ normal(0, 1);
  for(ii in 2:(T - 1)) {
    s_inter[ii] ~ normal(s_inter[ii - 1], 1);
  }

  // priors for order effect on origination
  s_ordeff ~ normal(0, s_ordscale);
  s_ordscale ~ normal(0, 0.5);  // really strong for no eff!

  // priors for effects of mass on survival
  s_b_1 ~ normal(0, 0.5);  // really strong for no eff!


  ////
  // observation
  // by time w/ random walk prior
  p_timeeff[1] ~ normal(0, 1);
  for(ii in 2:T) {
    //p_timeeff[ii] - p_timeeff[ii-1] ~ normal(0, 1);
    p_timeeff[ii] ~ normal(p_timeeff[ii - 1], 1);
  }

  // 
  p_funceff ~ normal(0, p_funcscale);
  p_funcscale ~ normal(0, 0.5);  // really strong for no eff!
  p_b_1 ~ normal(1, 1);  // slope mass; guessing positive


  for(n in 1:N) {
    target += state_space_lp(sight[n], origin[n, ], stay[n, ], p[n, ]);
  }
}
generated quantities {
  corr_matrix[D] o_Omega;
  corr_matrix[D] s_Omega;

  // convert back to a correlation matrix
  o_Omega = multiply_lower_tri_self_transpose(o_L_Omega);
  s_Omega = multiply_lower_tri_self_transpose(s_L_Omega);
}

