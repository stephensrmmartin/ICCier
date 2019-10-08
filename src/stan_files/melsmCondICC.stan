functions {
  vector ICC(vector mu_sd, vector sigma){
    int N = rows(sigma);
    vector[N] mu_var = mu_sd .* mu_sd;
    vector[N] sigma2 = sigma .* sigma;
    vector[N] ICC = mu_var ./ (mu_var + sigma2);
    return(ICC);
  }
  // doi: 10.1111/2041-210X.12225; but with predicted sigma_l instead.
  // TODO: Revisit this. May need to compute average x*S*x' for each /k/.
  // Requires finding all x entries for k, combining into one frame, then computing average.
  vector ICC_adjusted(int[] group, matrix x_loc_l1, matrix Omega, matrix mu_gamma_group_random_sd, vector sigma){
    int N = rows(sigma);
    int K = rows(mu_gamma_group_random_sd);
    int dim = rows(Omega);
    int Q_l1 = cols(x_loc_l1);
    matrix[dim,dim] cov[K];
    vector[N] numerator;
    vector[N] ICC_adjusted;

    for(k in 1:K){
      cov[k] = quad_form_diag(Omega,mu_gamma_group_random_sd[k]);
    }
    // Orig. author takes average numerator; we do not.
    // This is different: I use mean_ICC = E(Var / (Var + ErrVar))
    // They use mean_ICC = E(Var) / (E(Var) + ErrVar)
    // I use mean predicted ICC; they have ICC of mean values.
    for(n in 1:N){
      numerator[n] = x_loc_l1[n,] * cov[group[n]][1:Q_l1,1:Q_l1] * x_loc_l1[n,]';
    }
    ICC_adjusted = (numerator) ./ (numerator + (sigma .* sigma));
    return(ICC_adjusted);
  }
  vector ICC_adjusted_v2(int[] group, matrix x_loc_l1, matrix Omega, matrix mu_gamma_group_random_sd, vector sigma){
    int N = rows(sigma);
    int K = rows(mu_gamma_group_random_sd);
    int dim = rows(Omega);
    int Q_l1 = cols(x_loc_l1);
    int n_K[K] = rep_array(0,K);
    vector[N] numerator;
    vector[N] ICC_adjusted_out;
    for(n in 1:N){ // Count number of observations for each K.
      n_K[group[n]] += 1;
    }
    for(k in 1:K){
      int rows_K[n_K[k]]; // Rows indices for k
      int count = 1;
      real numerator_K;
      matrix[dim,dim] cov_K = quad_form_diag(Omega,mu_gamma_group_random_sd[k]);
      for(n in 1:N){
        if(group[n] == k){
          rows_K[count] = n;
          count += 1;
        }
      }
      //numerator_K = trace(x_loc_l1[rows_K,] * cov_K[1:Q_l1,1:Q_l1] * x_loc_l1[rows_K,]')/n_K[k];
      numerator_K = trace_quad_form(cov_K[1:Q_l1,1:Q_l1],x_loc_l1[rows_K,]')/n_K[k];
      numerator[rows_K] = rep_vector(numerator_K, n_K[k]);
    }
    ICC_adjusted_out = (numerator) ./ (numerator + sigma .* sigma);
    return(ICC_adjusted_out);

  }
}
data {
  int N;
  int K;
  int Q_l1;
  int Q_l2;
  int P_l1;
  int P_l2;
  int R_l2;
  int group[N];
  matrix[N,Q_l1] x_loc_l1;
  matrix[K,Q_l2] x_loc_l2;
  matrix[N,P_l1] x_sca_l1; // log(sigma) ~ x_sca_l1
  matrix[K,P_l2] x_sca_l2; // gammas ~ x_sca_l2 // For now, assume x_sca_l2 also predicts var(RE)
  matrix[K,R_l2] x_bet_l2;
  vector[N] y;
  int<lower=0,upper=1> adjust_icc;
}

parameters {
  // Location
  // real beta0; // Intercept-only location
  matrix[Q_l2, Q_l1] beta0;

  // Scale
  matrix[K, P_l1 + Q_l1] mu_gamma_group_random_z;
  cholesky_factor_corr[P_l1 + Q_l1] mu_gamma_group_random_cor_L;
  // vector[P_l1] mu_gamma_group_random_log_sd;
  matrix[R_l2, P_l1 + Q_l1] eta; // log(sd(u)) = x_sca_l2*eta
  matrix[P_l2, P_l1] gamma;
}

transformed parameters {
  matrix[K,P_l1 + Q_l1] mu_gamma_group_random_sd = exp(x_bet_l2*eta);
  matrix[K,P_l1 + Q_l1] mu_gamma_group_random;
  matrix[K,Q_l1] mu_group;
  matrix[K,P_l1] gamma_group;
  vector[N] yhat;
  vector[N] shat;

  for(k in 1:K){
    mu_gamma_group_random[k] = mu_gamma_group_random_z[k] * diag_pre_multiply(mu_gamma_group_random_sd[k],mu_gamma_group_random_cor_L)';
  }
  mu_group = x_loc_l2*beta0 + mu_gamma_group_random[,1:Q_l1];
  gamma_group = x_sca_l2*gamma + mu_gamma_group_random[,(Q_l1 + 1):];

  yhat = rows_dot_product(x_loc_l1, mu_group[group]);
  shat = exp(rows_dot_product(x_sca_l1,gamma_group[group]));
}

model {
  // Predictions
  // vector[N] yhat = mu_group[group];
  // vector[N] shat = exp(rows_dot_product(x_sca_l1,gamma_group[group]));

  // Priors
  to_vector(mu_gamma_group_random_z) ~ std_normal();
  mu_gamma_group_random_cor_L ~ lkj_corr_cholesky(1);
  to_vector(eta) ~ student_t(3,0,5);
  to_vector(gamma) ~ student_t(3,0,5);
  to_vector(beta0) ~ normal(0,10);

  // Likelihood
  y ~ normal(yhat,shat);
}

generated quantities {
//  vector[N] icc; = ICC(mu_gamma_group_random_sd[group,1], shat);
  corr_matrix[P_l1 + Q_l1] Omega = multiply_lower_tri_self_transpose(mu_gamma_group_random_cor_L);
  // vector[N] icc = adjust_icc ? ICC_adjusted(group, x_loc_l1, Omega, mu_gamma_group_random_sd, shat) : ICC(mu_gamma_group_random_sd[group,1], shat);
  vector[N] icc = adjust_icc ? ICC_adjusted_v2(group, x_loc_l1, Omega, mu_gamma_group_random_sd, shat) : ICC(mu_gamma_group_random_sd[group,1], shat);
  vector[N] log_lik;
  real icc_mean = mean(icc);
  real icc_sd = sd(icc);
  //vector ICC_adjusted(int[] group, matrix x_loc_l1, matrix Omega, matrix mu_gamma_group_random_sd, vector sigma){
 // vector[N] icc_adjusted = ICC_adjusted(group, x_loc_l1, Omega, mu_gamma_group_random_sd, shat);
  //real icc_adjusted_mean = mean(icc_adjusted);
  //real icc_adjusted_sd = sd(icc_adjusted);
  {
    for(n in 1:N){
      log_lik[n] = normal_lpdf(y[n] | yhat[n],shat[n]);
    }
  }
}
