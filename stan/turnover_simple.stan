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
      beta[t] ~ multi_normal(u[t] * gamma, quad_form_diag(Omega, tau));
    }
  }

  to_vector(gamma) ~ normal(0, 5);
  //gamma[, 1:D] can get more coherent priors
  
  
  for(n in 1:N) {
    sight[n, ] ~ bernoulli(pred[n, ]);
  }
}
generated quantities {
  int sight_tilde[N, T];  // observed presence
  
  for(n in 1:N) {
    for(t in 1:T) {
      sight_tilde[n, t] <- bernoulli_rng(pred[n, t]);
    }
  }
}
