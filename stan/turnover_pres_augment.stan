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
  real state_space_lp(int[] y, int I, real phi, row_vector pred, row_vector p) {
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
  vector[Nplus] mass_est;  // estimated masses
  real mu_mass;  // average of mass
  real<lower=0> sigma_mass; // variance of mass
  real b_1; // effect of mass on occurrence
  real b_2; // effect of mass on occurrence

  // effect associated with ecology
  matrix[D, T - 1] a_z; // part of non-centering
  cholesky_factor_corr[D] L_Omega;
  vector<lower=0>[D] tau;
  
  matrix[U, D] gamma; // effect of group level covariates

  real alpha_0;
  real alpha_1;
  vector[T] alpha_time;
  real<lower=0> sigma;

  real<lower=0,upper=1> phi;
}
transformed parameters {
  vector[M] mass_full;  // combined observed and estimated values of mass

  matrix[T - 1, D] a;  // effect associated with ecology
  matrix<lower=0,upper=1>[M, T] p;  // sampling probability
  matrix[M, T - 1] pred;  // occurrence probabilty 

  // combine observed and estmated
  mass_full[1:N] = mass;
  mass_full[(N + 1):M] = mass_est;

  // vectorized, non-centered, chol-decom group-level predictors
  a = u * gamma + (diag_pre_multiply(tau, L_Omega) * a_z)';

  // probability of observing
  for(n in 1:M) {
    for(t in 1:T) {
      p[n, t] = inv_logit(alpha_0 + alpha_time[t] + alpha_1 * mass_full[n]);
    }
  }

  // probability of occurring
  for(t in 1:(T-1)) {
    for(n in 1:M) {
      pred[n, t] = inv_logit(a[t, statea[n]] + 
          b_1 * mass_full[n] + b_2 * (mass_full[n] ^ 2));
    }
  }
}
model {
  // imputation of unobserved masses
  mass ~ normal(mu_mass, sigma_mass);
  mass_est ~ normal(mu_mass, sigma_mass);
  mu_mass ~ normal(0, 0.1);
  sigma_mass ~ normal(0, 1);
  
  // effects of ecologies
  to_vector(a_z) ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(2);
  tau ~ normal(0, 1);
  to_vector(gamma) ~ normal(0, 1);

  // effects of mass on occurrence
  alpha_0 ~ normal(0, 5);
  alpha_1 ~ normal(0, 1);
  alpha_time ~ normal(0, sigma);
  sigma ~ normal(0, 1);

  // effects of mass on observation
  b_1 ~ normal(0, 1);
  b_2 ~ normal(0, 1);

  for(n in 1:M) {
    target += state_space_lp(sighta[n], I[n], phi, pred[n, ], p[n, ]);
  }
}
generated quantities {
  corr_matrix[D] Omega;

  // convert back to a correlation matrix
  Omega = multiply_lower_tri_self_transpose(L_Omega);
}
