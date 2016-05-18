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
    int prod_term;

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
        prod_term <- 1 - z[1];
        for(j in 2:S) {
          prod_term <- prod_term * (1 - z[j - 1]);
          sl <- sl + bernoulli_log(z[j], (z[j - 1] * pred[j]) + 
              prod_term * pred[j]);
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

  int sight[N, T];  // observed presence
  matrix[N, D] x;  // matrix of indiv-level covariates
  row_vector[U] u[T];  // matrix of group-level covariates
}
parameters {
  corr_matrix[D] Omega;
  vector<lower=0>[D] tau;
  matrix[U, D] gamma;
  vector[D] beta[T];
  
  real<lower=0> lambda[U];
  vector<lower=0>[U] phi[D];

  vector[T] p_norm; 
  real p_mu;
  real<lower=0> p_sigma;
}
transformed parameters {
  vector<lower=0,upper=1>[T] p;  // sampling probability
  matrix[N, T] pred; 

  for(t in 1:T) {
    p[t] <- inv_logit(p_mu + p_sigma * p_norm[t]);
  }
  
  // assemble predictor w/ intercept + effects of indiv-level covariates
  for(t in 1:T) {
    for(n in 1:N) {
      // non-centered parameterization following Betacourt and Girolami
      pred[n, t] <- inv_logit(x[n, ] * beta[t, ]);
    }
  }
}
model {
  tau ~ cauchy(0, 1);
  Omega ~ lkj_corr(2);

  {
    matrix[D, D] Sigma_beta;
    Sigma_beta <- quad_form_diag(Omega, tau);

    for(t in 1:T) {
      beta[t] ~ multi_normal(u[t] * gamma, Sigma_beta);
    }
  }

  for(k in 1:U) {
    for(j in 1:D) {
      gamma[k, j]  ~ normal(0, lambda[k] * phi[j][k]);
      phi[j][k] ~ cauchy(0, 1);
    }
  }

  p_mu ~ normal(0, 1);
  p_sigma ~ cauchy(0, 1);


  for(n in 1:N) {
    sight[n] ~ state_space(pred[n, ], p);
  }
}
//generated quantities {
//  int true_tilde[N, T]
//  int sight_tilde[N, T];  // observed presence
//
//  {
//    int prod_state;
//
//    for(n in 1:N) {
//      true_tilde[n, 1] <- bernoulli_rng(pred[n, 1]);
//      sight_tilde[n, 1] <- bernoulli_rng(p[1] * true_tilde[n, 1]);
//      prod_state <- (1 - true_tide[n, 1]);
//      for(t in 2:T) {
//        prod_state <- prod_state * (1 - true_tilde[n, t]);
//        true_tilde[n, t] <- bernoulli_rng(true_tilde[n, t - 1] * pred[n, t] + 
//            prod_state * pred[n, t]);
//        sight_tilde[n, t] <- bernoulli_rng(p[t] * true_tilde[n, t - 1]);
//      }
//    }
//  }
//}
