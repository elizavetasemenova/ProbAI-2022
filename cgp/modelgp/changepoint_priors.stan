data {
  int<lower=1> N;
  real x[N];
  
  real<lower=0> rho;
  real<lower=0> alpha;
  
  real theta;
  
  real g;
  real k2;
}

transformed data {
  
  matrix[N, N] L;
  real d1;
  real d2;
  
  {
    
    matrix[N, N] K;
    matrix[N, N] K1;
    matrix[N, N] K2 =  k2 * cov_exp_quad(x, alpha, rho);
    
     for (i in 1:N) {
      for (j in 1:N) {
        d1 = x[i]-theta;
        d2 = x[j]-theta;
        K1[i,j] = 0; 
        //K2[i,j] = 0; 
        //K2[i,j] = k2 * alpha^2* d1^2 * d2^2 * exp(-(d1-d2)^2 / (rho^2));
        
        K[i, j] = (1 - inv_logit(g*d1))*K1[i, j]*(1 - inv_logit(g*d2)) + inv_logit(g*d1)*K2[i, j]*inv_logit(g*d2);
    }
    {
      real diag = K[i, i];
      K[i, i] = diag + 1e-10;
    }
  }
  
  L = cholesky_decompose(K);
    
  }
 
}

parameters {}
model {}

generated quantities {
  vector[N] f = multi_normal_cholesky_rng(rep_vector(0, N), L);
}
