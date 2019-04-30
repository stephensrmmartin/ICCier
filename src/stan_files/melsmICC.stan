functions {
  vector ICC(vector mu_sd, vector sigma){
    int N = rows(sigma);
    vector[N] mu_var = mu_sd .* mu_sd;
    vector[N] sigma2 = sigma .* sigma;
    return(mu_var ./ (mu_var + sigma2));
  }
}
data {
  int N;
  int K;
  int P_l1;
  int P_l2;
  int group[N];
  matrix[N,P_l1] x_sca_l1; // log(sigma) ~ x_sca_l1
  matrix[K,P_l2] x_sca_l2; // gammas ~ x_sca_l2 // For now, assume x_sca_l2 also predicts var(RE)
  vector[N] y;
}

parameters {
  // Location
  real beta0; // Intercept-only location

  // Scale
  matrix[K, P_l1 + 1] mu_gamma_group_random_z;
  cholesky_factor_corr[P_l1 + 1] mu_gamma_group_random_cor_L;
  // vector[P_l1] mu_gamma_group_random_log_sd;
  matrix[P_l2, P_l1 + 1] eta; // log(sd(u)) = x_sca_l2*eta
  matrix[P_l2, P_l1] gamma;
}

transformed parameters {
  matrix[K,P_l1 + 1] mu_gamma_group_random_sd = exp(x_sca_l2*eta);
  matrix[K,P_l1 + 1] mu_gamma_group_random;
  vector[K] mu_group;
  matrix[K,P_l1] gamma_group;
  vector[N] yhat;
  vector[N] shat;

  for(k in 1:K){
    mu_gamma_group_random[k] = mu_gamma_group_random_z[k] * diag_pre_multiply(mu_gamma_group_random_sd[k],mu_gamma_group_random_cor_L)';
  }
  mu_group = beta0 + mu_gamma_group_random[,1];
  gamma_group = x_sca_l2*gamma + mu_gamma_group_random[,2:];

  yhat = mu_group[group];
  shat = exp(rows_dot_product(x_sca_l1,gamma_group[group]));
}

model {
  // Predictions
  // vector[N] yhat = mu_group[group];
  // vector[N] shat = exp(rows_dot_product(x_sca_l1,gamma_group[group]));

  // Priors
  to_vector(mu_gamma_group_random_z) ~ std_normal();
  mu_gamma_group_random_cor_L ~ lkj_corr_cholesky(1);
  to_vector(eta) ~ std_normal();
  to_vector(gamma) ~ std_normal();
  beta0 ~ normal(0,10);

  // Likelihood
  y ~ normal(yhat,shat);
}

generated quantities {
  vector[N] icc = ICC(mu_gamma_group_random_sd[group,1], shat);
  vector[N] log_lik;
  corr_matrix[P_l1 + 1] Omega = multiply_lower_tri_self_transpose(mu_gamma_group_random_cor_L);
  {
    for(n in 1:N){
      log_lik[n] = normal_lpdf(y[n] | yhat[n],shat[n]);
    }
  }
}
