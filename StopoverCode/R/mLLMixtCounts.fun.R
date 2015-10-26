#' @export
#-- Function to create vectors and matrices of parameters --#

mLLMixtCounts.fun <- function(invect)
{
  if(p.m == "common"){ 
    N.in <- exp(invect[1:S]); 
    par.index <- S
    p.in <- matrix(1, nrow = S, ncol = TIME)
  } else { 
    N.in <- exp(invect[1:S])
    cvec.in <- invect[(S+1):(S+qp+1)]#c0.in <- invect[S+1]; 
    #c1.in <- invect[S+2]; 
    par.index <- S + qp + 1
    p.temp <- matrix(cvec.in[1], nrow = S, ncol = TIME)
    for(q in 1:qp) p.temp <- p.temp + cvec.in[q+1]*cov.p[,,q] #cov.p is now an array of SxTxqp dimensions
    p.in <- expo(p.temp)
  }
  
  if (M>1) {
    d0.in <- invect[(par.index + 1):(par.index + M - 1)]; 
    par.index <- par.index + M - 1
    d1.in <- NULL
    w.in <- matrix(0,S,M)
    if(w.m == "common") { for(i in 1:S) w.in[i,] <- lgp(d0.in)}
    if(w.m == "cov") { d1.in <- invect[par.index + 1]; for(i in 1:S) w.in[i,] <- lgp(d0.in + d1.in * cov.w[i]); par.index <- par.index + 1}
  } else {
    w.in <- matrix(1,S,M)
  }
  
  b0.in <- invect[(par.index + 1):(par.index + M)]; par.index <- par.index + M
  b1.in <- NULL
  if(mu.m == "cov") { 
    b1.in <- invect[(par.index + 1)]; 
    par.index <- par.index + 1
  }
  
  sigma.in <- exp(invect[par.index + 1]); 
  par.index <- par.index + 1
  if(sigma.m == "het") { 
    sigma.in <- c(sigma.in, exp(invect[(par.index + 1):(par.index + M - 1)])); 
    par.index <- par.index + M - 1
  }
  sigma.in <- matrix(sigma.in, S, M, byrow = TRUE)
  
  mu.in <- matrix(NA,S,M)
  betta.in <- matrix(0,S,TIME)
  
  if(mu.m == "common") {
    for(i in 1:S) {
      for(m in 1:M) {
        mu.in[i,m] <- exp(b0.in[m])
        betta.in[i,] <- betta.in[i,]+c(w.in[i,m]*c(pnorm(1,mean=mu.in[i,m],sd=sigma.in[i,m]),pnorm(2:(TIME-1),mean=mu.in[i,m],sd=sigma.in[i,m])-pnorm(1:(TIME-2),mean=mu.in[i,m],sd=sigma.in[i,m]),1-pnorm(TIME-1,mean=mu.in[i,m],sd=sigma.in[i,m])))
      }
    }
  }
  
  if(mu.m == "cov") {
    for(i in 1:S) {
      for(m in 1:M) {
        mu.in[i,m] <- exp(b0.in[m] + b1.in*cov.mu[i])
        betta.in[i,] <- betta.in[i,]+c(w.in[i,m]*c(pnorm(1,mean=mu.in[i,m],sd=sigma.in[i,m]),pnorm(2:(TIME-1),mean=mu.in[i,m],sd=sigma.in[i,m])-pnorm(1:(TIME-2),mean=mu.in[i,m],sd=sigma.in[i,m]),1-pnorm(TIME-1,mean=mu.in[i,m],sd=sigma.in[i,m])))
      }
    }
  }
  
  a0.in <- invect[par.index + 1]; par.index <- par.index + 1
  a1.in <- NULL
  if((phi.m == "logit.t")|(phi.m == "logit.a")) { 
    a1.in <- invect[par.index + 1]; par.index <- par.index + 1
  }
  a2.in <- NULL
  if((phi.m == "quadr.t")|(phi.m == "quadr.a")) {
    a1.in <- invect[par.index + 1];  
    a2.in <- invect[par.index + 2]; 
    par.index <- par.index + 2
  }
  
  phi.in <- array(0,c(TIME-1,TIME-1,S))
  
  if(phi.m == "const") {
    for(i in 1:S) {	
      phi.in[,,i] <- expo(a0.in)
      phi.in[,,i][lower.tri(phi.in[,,i])]<-0
    }
  }
  
  if(phi.m == "logit.t") {
    for(i in 1:S){	
      for(b in 1:(TIME-1)){
        phi.in[b,,i] <- expo(a0.in + a1.in*cov.phi[i,])
      }
      phi.in[,,i][lower.tri(phi.in[,,i])]<-0
    }
  }
  
  if(phi.m == "logit.a") {
    for(i in 1:S) {
      for(j in 1:(TIME-1)) {
        phi.in[j,(j:(TIME-1)),i] <- expo(a0.in + a1.in*cov.phi[i,])[1:(TIME-j)]
      }
      phi.in[,,i][lower.tri(phi.in[,,i])]<-0
    }
  }
  
  if(phi.m == "quadr.t") {
    for(i in 1:S) {	
      for(b in 1:(TIME-1)) {
        phi.in[b,,i] <- expo(a0.in + a1.in*cov.phi[i,] + a2.in*cov.phi[i,]^2)
      }
      phi.in[,,i][lower.tri(phi.in[,,i])]<-0
    }
  }
  
  if(phi.m == "quadr.a"){
    for(i in 1:S) {
      for(j in 1:(TIME-1)) {
        phi.in[j,(j:(TIME-1)),i] <- expo(a0.in + a1.in*cov.phi[i,] + a2.in*cov.phi[i,]^2)[1:(TIME-j)]
      }
      phi.in[,,i][lower.tri(phi.in[,,i])]<-0
    }
  }
  
  LL.MixtCounts.c(N.in, p.in, betta.in, phi.in)
}
