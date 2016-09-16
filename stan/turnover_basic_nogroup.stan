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
  real<lower=0,upper=1> phi;
  corr_matrix[D] Omega;
  vector<lower=0>[D] tau;
  vector[D] gamma;
  //matrix[U, D] gamma;
  vector[D] beta[T-1];
}
transformed parameters {
  matrix[N, T-1] pred;
  
  for(t in 1:(T-1)) {
    for(n in 1:N) {
      pred[n, t] = inv_logit(x[n, ] * beta[t, ]);
    }
  } 
}
model {
  tau ~ cauchy(0, 1);
  Omega ~ lkj_corr(2);
  
  {
    matrix[D, D] Sigma_beta;
    Sigma_beta = quad_form_diag(Omega, tau);
    
    beta ~ multi_normal(gamma, Sigma_beta);
    //for(t in 1:T) {
    //  beta[t] ~ multi_normal(u[t] * gamma, Sigma_beta);
    //}
  }

  to_vector(gamma) ~ normal(0, 1);

  {
    int prod_state;

    for(n in 1:N) {
      prod_state = (1 - sight[n, 1]);
      sight[n, 1] ~ bernoulli(phi);
      for(t in 2:T) {
        prod_state = prod_state * (1 - sight[n, t - 1]);
        sight[n, t] ~ bernoulli(sight[n, t - 1] * pred[n, t - 1] + 
            prod_state * pred[n, t - 1]);
      }
    }
  }
}
