data {
  int N_treatment;                            // total number of observations (treatment)
  int N_conc;                                 // number of unique concentrations
  vector[N_conc] x;                           // vector of unique treatment concentrations
  vector[N_treatment] y_treatment;            // vector of observations (treatment)
  int x_index[N_treatment];                   // index of treatment dose for each observation
  
  int N_control;                             //  total number of control observations
  vector[N_control] y_control;               //  vector of control observations
}
transformed data {
  real min_x = min(x);                      //  minimal concentration
  real max_x = max(x);                      //  maximal concentration
}
parameters {
  real mu;                                  //  baseline (mean of conrols and treatments at zero dose)                          
  real<lower=0> sigma_rep;                  //  replicate-to-replicate variation   
  real<lower=0> sigma;                      //  curve uncertainty 
  
  real<lower=0> rho;                        //  GP-lengthscale
  real<lower=0, upper=1> theta_raw;         
  real<lower=0> eta;                        //  GP-amplitude
  
  vector[N_conc] y;                         //  mean response, indexed by unique concentrations
  
  real<lower=1> nu;                         //  degrees of freedom
  real<lower=0.1> g;                        //  changepoint steepness 
}
transformed parameters {
  real theta;                               //  PoD
  vector[N_treatment] y_mean;               //  mean response, indexed by number of (possibly, repeated) concentrations
  
  for (i in 1:N_treatment) {
      y_mean[i] = y[x_index[i]];
  }
  
  theta = min_x + theta_raw * (max_x - min_x);
  
}
model {
  
  // Priors
  mu ~ normal(0, 0.1);
  sigma ~ inv_gamma(1, 2);
  sigma_rep ~ inv_gamma(1, 0.1);
  eta ~ normal(1, 1);
  //rho ~ gamma(50, 20);
  rho ~ gamma(30, 20);
  theta_raw ~ uniform(0, 1);
  nu ~ gamma(2,0.1);
  g ~ gamma(10,1);
  
  // GP
  {
    matrix[N_conc, N_conc] K1;
    matrix[N_conc, N_conc] K2;
    matrix[N_conc, N_conc] K;
    matrix[N_conc, N_conc] L;
    real d1;
    real d2;
    
    // Construct covariance matrices
    for (i in 1:N_conc) {
      for (j in 1:N_conc) {
        d1 = x[i]-theta;
        d2 = x[j]-theta;
        K1[i,j] = 0; 
        K2[i,j] = eta^2* d1^2 * d2^2 * exp(-(d1-d2)^2 / (rho^2));
        
        K[i, j] = (1 - inv_logit(g*d1))*K1[i, j]*(1 - inv_logit(g*d2)) + inv_logit(g*d1)*K2[i, j]*inv_logit(g*d2);
    }
    {
      real diag = K[i, i];
      K[i, i] = diag + sigma^2;
    }
  }
  L = cholesky_decompose(K);
  
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
