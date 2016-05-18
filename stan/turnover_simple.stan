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
}
transformed parameters {
  matrix[N, T] pred;

  for(t in 1:T) {
    for(n in 1:N) {
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

  to_vector(gamma) ~ normal(0, 5);
  
  {
    int prod_state;

    for(n in 1:N) {
      prod_state <- (1 - sight[n, 1]);
      sight[n, 1] ~ bernoulli(pred[n, 1]);
      for(t in 2:T) {
        prod_state <- prod_state * (1 - sight[n, t - 1]);
        sight[n, t] ~ bernoulli(sight[n, t - 1] * pred[n, t] + 
            prod_state * pred[n, t]);
      }
    }
  }
}
//generated quantities {
//  int sight_tilde[N, T];  // observed presence
//
//  {
//    int prod_state;
//
//    for(n in 1:N) {
//      prod_state <- (1 - sight[n, 1]);
//      sight_tilde[n, 1] <- sight[n, 1];
//      for(t in 2:T) {
//        prod_state <- prod_state * (1 - sight_tilde[n, t - 1]);
//        sight_tilde[n, t] <- bernoulli_rng(sight_tilde[n, t - 1] * pred[n, t] + 
//            prod_state * pred[n, t]);
//      }
//    }
//  }
//}
