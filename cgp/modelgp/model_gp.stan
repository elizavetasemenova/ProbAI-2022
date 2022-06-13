data {
  int N_treatment;                            // total number of observations (treatment)
  int N_conc;                                 // number of unique concentrations
  //vector[N_conc] x;                           // vector of unique treatment concentrations
  real x[N_conc];
  vector[N_treatment] y_treatment;            // vector of observations (treatment)
  int x_index[N_treatment];                   // index of treatment dose for each observation
  
  int N_control;                             //  total number of control observations
  vector[N_control] y_control;               //  vector of control observations
}
parameters {
  real mu;                                  //  baseline (mean of conrols and treatments at zero dose)                          
  real<lower=0> sigma_rep;                  //  replicate-to-replicate variation   
  real<lower=0> sigma;                      //  curve uncertainty 
  real<lower=0> rho;                        //  GP-lengthscale
  real<lower=0> eta;                        //  GP-amplitude
  vector[N_conc] y;                         //  mean response, indexed by unique concentrations
  real<lower=1> nu;                         //  degrees of freedom
}
transformed parameters {
  vector[N_treatment] y_mean;               //  mean response, indexed by number of (possibly, repeated) concentrations
  
  for (i in 1:N_treatment) {
      y_mean[i] = y[x_index[i]];
  }
}
model {
  
  // Priors
  mu ~ normal(0, 0.1);
  sigma ~ inv_gamma(1, 2);
  sigma_rep ~ inv_gamma(1, 0.1);
  eta ~ normal(1, 1);
  //rho ~ gamma(50, 20);
  //rho ~ gamma(40, 20);
  rho ~ gamma(30, 20);
  //theta_raw ~ uniform(0, 1);
  nu ~ gamma(2,0.1);
  //g ~ gamma(10,1);
  
  // GP
  {
    // Construct covariance matrix
    matrix[N_conc, N_conc] K = cov_exp_quad(x, eta, rho) + diag_matrix(rep_vector(sigma^2, N_conc));
    matrix[N_conc, N_conc] L = cholesky_decompose(K);
    
    y ~ multi_normal_cholesky(rep_vector(mu, N_conc), L);
  }
  
   // likelihood - controls
    y_control ~ student_t(nu, mu, sqrt(sigma^2 + sigma_rep^2));
    // likelihood - treatment
    y_treatment ~ student_t(nu, y_mean, sigma_rep);
    
}

generated quantities{
   vector[N_treatment] y_pred; // Vector of observations (non-control)
   
   for (i in 1:N_treatment) {
     y_pred[i] = student_t_rng(3, y_mean[i], sigma_rep);
   }
   
}
