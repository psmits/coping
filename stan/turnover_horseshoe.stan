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
  matrix[U, D] gamma_std;
  vector[D] beta[T];

  real<lower=0> lambda[D];
  vector<lower=0>[D] phi[U];
}
transformed parameters {
  matrix[U, D] gamma;
  matrix[N, T] pred;
  
  for(k in 1:U) {
    for(j in 1:D) {
      gamma[k, j] <- gamma_std[k, j] * lambda[j] * phi[k][j];
    }
  }

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

  to_vector(gamma_std) ~ normal(0, 1);
  for(i in 1:U) phi[i] ~ cauchy(0, 1);
  lambda ~ cauchy(0, 1);


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
