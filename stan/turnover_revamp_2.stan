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
  real state_space_lp(int[] y, int I, real phi, row_vector origin, row_vector stay, row_vector p) {
    // idea
    //  indicator for augmented points
    //  they get lp, ft, lt all calculated differently
    //  proceed as normal
    int ft;
    int lt;
    int S;
    int i;
    int prod_term;
    int lp_size;
    int second_start;

    S = size(y);
    i = 1;

    if (I == 1) {
      ft = first_capture(y);
      lt = last_capture(y);
      lp_size = ft * (S - lt + 1);
    } else {
      ft = S;
      lp_size = (S * (S + 1)) / 2;
    }

    {
      vector[lp_size] lp;

      // have to go through all possible event histories for each species
      for(t_first_alive in 1:ft) {
        if(I == 1) {
          second_start = lt;
        } else {
          second_start = t_first_alive;
        }
        
        for (t_last_alive in second_start:S) {
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
              gg[j - 1] = bernoulli_lpmf(z[j] | (z[j - 1] * stay[j - 1]) + 
                  prod_term * origin[j - 1]);
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
  int M;  // total sample + augment
  int D;  // number of indiv-level predictors
  int U;  // number of group-level predictors

  int I[M];

  int sighta[M, T];  // observed presence
  int statea[M];  // vector of ecology state

  vector[N] mass;

  matrix[T - 1, U] u;  // matrix of group-level covariates
}
transformed data {
  int Nplus;

  Nplus = M - N;
}
parameters {
  // general
  vector[Nplus] mass_est;  // estimated masses
  real mu_mass;  // average of mass
  real<lower=0> sigma_mass; // variance of mass
  
  
  // origination
  real o_b_1; // effect of mass on occurrence
  real o_b_2; // effect of mass on occurrence

  // effect associated with ecology
  matrix[D, T - 1] o_a_z; // part of non-centering
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
  real alpha_0;
  real alpha_1;
  vector[T] alpha_time;
  real<lower=0> sigma;

  // initial state
  real<lower=0,upper=1> phi;
}
transformed parameters {
  vector[M] mass_full;  // combined observed and estimated values of mass

  matrix[T - 1, D] o_a;  // origin: effect associated with ecology
  matrix[T - 1, D] s_a;  // survival: effect associated with ecology
  matrix<lower=0,upper=1>[M, T] p;  // sampling probability
  matrix[M, T - 1] origin;  // occurrence probabilty 
  matrix[M, T - 1] stay;  // occurrence probabilty 

  // combine observed and estmated
  mass_full[1:N] = mass;
  mass_full[(N + 1):M] = mass_est;

  // vectorized, non-centered, chol-decom group-level predictors
  o_a = u * o_gamma + (diag_pre_multiply(o_tau, o_L_Omega) * o_a_z)';
  s_a = u * s_gamma + (diag_pre_multiply(s_tau, s_L_Omega) * s_a_z)';

  // probability of observing
  for(n in 1:M) {
    for(t in 1:T) {
      p[n, t] = inv_logit(alpha_0 + alpha_time[t] + alpha_1 * mass_full[n]);
    }
  }

  // probability of occurring
  for(t in 1:(T-1)) {
    for(n in 1:M) {
      origin[n, t] = inv_logit(o_a[t, statea[n]] + 
          o_b_1 * mass_full[n] + o_b_2 * (mass_full[n] ^ 2));
      stay[n, t] = inv_logit(s_a[t, statea[n]] + 
          s_b_1 * mass_full[n] + s_b_2 * (mass_full[n] ^ 2));
    }
  }
}
model {
  // imputation of unobserved masses
  //mass ~ normal(mu_mass, sigma_mass);
  mass_est ~ normal(0, 0.5);
  //mu_mass ~ normal(0, 0.1);
  //sigma_mass ~ normal(0, 1);

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

  // effects of mass on occurrence
  alpha_0 ~ normal(0, 5);
  alpha_1 ~ normal(0, 1);
  alpha_time ~ normal(0, sigma);
  sigma ~ normal(0, 1);

  // effects of mass on observation
  o_b_1 ~ normal(0, 1);
  o_b_2 ~ normal(0, 1);
  s_b_1 ~ normal(0, 1);
  s_b_2 ~ normal(0, 1);

  for(n in 1:M) {
    target += state_space_lp(sighta[n], I[n], phi, origin[n, ], stay[n, ], p[n, ]);
  }
}
generated quantities {
  corr_matrix[D] o_Omega;
  corr_matrix[D] s_Omega;

  // convert back to a correlation matrix
  o_Omega = multiply_lower_tri_self_transpose(o_L_Omega);
  s_Omega = multiply_lower_tri_self_transpose(s_L_Omega);
}
