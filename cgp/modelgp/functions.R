#---------------------------------------------------------------------------------------------------------------
# kriging: CalcSigma
#---------------------------------------------------------------------------------------------------------------
CalcSigma <- function(x1, x2, N_s, eta_, rho_, sigma_){
  N1 <- length(x1)
  N2 <- length(x2)
  Sigma <- array(data=0, dim = c(N_s, N1, N2))
  sigma2 <- sigma_ * sigma_
  
  if (N1 == N2){
    for (i in 1:N1)
      Sigma[, i, i] = sigma2[i]
  }
  
  for (i in 1:N_s){
    for (j in 1:N1){
      for (k in 1:N2){
        
        Sigma[i, j, k] = Sigma[i, j, k] + eta_[i]^2 * exp(-(x1[j]-x2[k])^2/(rho_[i]^2))
        
      }
    }
  }
  return(Sigma)
}

#---------------------------------------------------------------------------------------------------------------
# kriging: CalcSigma_g
#---------------------------------------------------------------------------------------------------------------
CalcSigma_g <- function(x1, x2, N_s, eta_, rho_, sigma_, theta_, k_){
  N1 <- length(x1)
  N2 <- length(x2)
  Sigma <- array(data=0, dim = c(N_s, N1, N2))
  sigma2 <- sigma_ * sigma_
  
  if (N1 == N2){
    for (i in 1:N1)
      Sigma[, i, i] = sigma2[i]
  }
  
  for (i in 1:N_s){
    for (j in 1:N1){
      for (k in 1:N2){
        
        d1 = x1[j] - theta_[i]
        d2 = x2[k] - theta_[i]
        K2 = eta_[i]^2 * d1^2 * d2^2 * exp(-(d1-d2)^2/(rho_[i]^2))
        
        Sigma[i, j, k] = Sigma[i, j, k] +  inv_logit(k_[i]*d1)* K2 *inv_logit(k_[i]*d2)
        
      }
    }
  }
  return(Sigma)
}

#---------------------------------------------------------------------------------------------------------------
# pred y 4PL
#---------------------------------------------------------------------------------------------------------------
y_4pl <- function(x, N_s, d_, a_, b_, c_){
  N <- length(x)
  
  y_pred <- matrix(NA, nrow = N_s, ncol = N)
  
  for (i in 1:N_s){
    for (j in 1:N){
      
      y_pred[i, j] = d_[i] + (a_[i]-d_[i])/(1 + exp(- b_[i] * (x[j] - c_[i])));  
      
    }
  }
  return(y_pred)
}

#---------------------------------------------------------------------------------------------------------------
# inverse logit
#---------------------------------------------------------------------------------------------------------------
inv_logit <- function(x){
  1/ (1 + exp(-x))
}

#---------------------------------------------------------------------------------------------------------------
# mode
#---------------------------------------------------------------------------------------------------------------
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#---------------------------------------------------------------------------------------------------------------
# DC50
#---------------------------------------------------------------------------------------------------------------

DC50_func <- function(s0, gp, x, eps=0.1){
  DC50_1_y <- min(gp) + 0.5* (s0 - min(gp))
  diff <- abs(gp - DC50_1_y)
  DC50_1 <- min(x[ which(diff - min(diff) < eps)])
  return(list(DC50_1))
}


#---------------------------------------------------------------------------------------------------------------
# Plot GP prior realizations
#---------------------------------------------------------------------------------------------------------------

plot_gp_prior_realizations <- function(fit, xs, title) {
  samples <- extract(fit)
  I <- length(samples$f[,1])
  
  plot_idx <- seq(1, I, 50)
  N <- length(plot_idx)
  
  plot(1, type="n", 
       xlab=expression('log'[10]*'-concentration [M]'),
       ylab="response", 
       main=title,
       xlim=c(-8, -5), ylim=c(-3.5, 3.5),
       xaxt='n', yaxt='n',
       cex.lab = 1,
       cex.main=1)
  
  qs <- apply(samples$f, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
  
  col3 <- "#00998a" 
  
  polygon(c(xs,rev(xs)),c(qs[2, ],rev(qs[1,])), col=alpha(col3,0.2), border=NA)
  
  for (n in 1:N)
    lines(xs, samples$f[plot_idx[n],], col=alpha("coral3",0.2), lwd=2)
  
  legend("topleft", 
         legend = c("GP draws"),
         col = c("coral3"),
         pch = c(NA), 
         lty=c(1),
         lwd=c(1),
         bty = "n", 
         pt.cex = 1, 
         cex = 1, 
         text.col = "black",
         inset = c(0.05, 0.01),
         border = c(NA)
  )
  
}
