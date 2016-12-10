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
  real state_space_lp(int[] y, int I, real psi, real phi, real gamma, real p) {
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
          sl = bernoulli_lpmf(z[1] | psi);
          prod_term = 1 - z[1];

          // z as function of trait state
          {
            vector[S - 1] gg;
            for(j in 2:S) {
              prod_term = prod_term * (1 - z[j - 1]);
              gg[j - 1] = bernoulli_lpmf(z[j] | (z[j - 1] * gamma) + 
                  prod_term * psi);
            }
            sl = sl + sum(gg);
          }

          // finally, y as function of z and p
          {
            vector[S] hh;
            for(k in 1:S) {
              hh[k] = bernoulli_lpmf(y[k] | z[k] * p);
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

  int I[M];  // augmented? 1 = no, 0 = yes

  int sighta[M, T];  // observed presence

//  vector[N] mass;
}
parameters {
  real<lower=0,upper=1> psi;
  real<lower=0,upper=1> phi;
  real<lower=0,upper=1> gamma;
  real<lower=0,upper=1> p;

//  vector[N] mass_est;
}
transformed parameters {
//  vector[M] mass_full;
//
//  mass_full[1:N] = mass;
//  mass_full[(N + 1):M] = mass_est;
}
model {
//  mass ~ normal(mu_mass, sigma_mass);
//  mass_est ~ normal(mu_mass, sigma_mass);

  for(n in 1:M) {
    target += state_space_lp(sighta[n], I[n], psi, phi, gamma, p);
  }
}
