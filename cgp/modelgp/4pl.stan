data {
  int N_treatment;                            // Total number of observations (treatment)
  int N_conc;                                 // Number of unique concentrations
  vector[N_conc] x;                           // Vector of unique treatment log10-concentrations
  vector[N_treatment] y_treatment;            // Vector of observations (treatment)
  int x_index[N_treatment];                   // Index of treatment dose for each observation
  
  int N_control;                             // Total number of control observations
  vector[N_control] y_control;               // Vector of control observations
}
transformed data {
}
parameters {
  real<lower=min(y_treatment), upper=max(y_control)> d;   // S0
  real<lower=min(y_treatment), upper=max(y_control)> a;   // Sinf = Dmax
  real b;                                                 // Hill coefficient
  real<lower=min(x), upper=max(x)> c;                     //log(DC50)
  real<lower=0> sigma_rep;          
  real<lower=0> sigma; 
  vector[N_conc] y;
}
transformed parameters {
  vector[N_treatment] y_mean;
  vector[N_conc] y_4pl;
  
  for (i in 1:N_conc){
    y_4pl[i] = d + (a-d)/(1 + exp(- b * (x[i] - c)));  
  }
  
  for (i in 1:N_treatment) {
      y_mean[i] = y[x_index[i]];
  }
}

model {
  
   // Priors
   
   sigma ~ inv_gamma(1, 2);
   sigma_rep ~ inv_gamma(1, 0.1);
   d ~ normal(mean(y_control), 0.1);
   a ~ normal(min(y_treatment), 0.1);
   b ~ normal(1, 5);
   
   {
     y ~ normal(y_4pl, sigma);
   }

   // likelihood
   y_control ~ normal(d, sqrt(sigma^2 + sigma_rep^2));
   y_treatment ~ normal(y_mean, sigma_rep);
    
}

generated quantities{
   vector[N_treatment] y_pred; // Vector of observations (non-control)
   
   for (i in 1:N_treatment) {
     y_pred[i] = normal_rng(y_mean[i], sigma_rep);
   }
}
