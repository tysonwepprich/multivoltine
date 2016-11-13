#' @export
#---- Function to fit the model ----#
mLLMixtCounts.fit <- function(p.m = "common", w.m = "common", mu.m = "common", sigma.m = "hom", phi.m = "const")
{
  ### Maximise log-likelihood ###
  this.fit <- NA
  this.fit <- try(optim(par=pars.start,fn=mLLMixtCounts.fun,method="L-BFGS-B",
                        hessian=F,control=list(maxit=10000,trace=1)), silent = TRUE) #maxit usually 10000, trace = 6
                        # hessian=T,control=list(maxit=10000,trace=1)), silent = TRUE) #maxit usually 10000, trace = 6
  # this.fit <- try(optim(par=pars.start,fn=mLLMixtCounts.fun,method="L-BFGS-B",
  #                       lower = c(0, -Inf, -Inf, -Inf, 0, -Inf, -Inf, logit(.2), -Inf, -Inf),
  #                       upper = c(10000000, Inf, Inf, Inf, log(31), Inf, Inf, logit(.8), Inf, Inf),
  #                       hessian=F,control=list(maxit=10000,trace=1)), silent = TRUE) #maxit usually 10000, trace = 6
  # 
  #  "N" = N.st, "cvec" = cvec.st,  "d0" = d0.st, "d1" = d1.st, 
  #"b0" = b0.st, "b1" = b1.st, "sigma" = sigma.st, "a0" = a0.st, 
  #"a1" = a1.st, "a2" = a2.st) 
  
  
  if(is.list(this.fit)) {
    outvect <- this.fit$par
    #VC <- solve(this.fit$hessian)
    
    #c0.out <- NULL; 
    #c1.out <- NULL
    cvec.out <- NULL 
    if(p.m == "common") { 
      N.out <- exp(outvect[1:S]); 
      par.index <- S
      p.out <- matrix(1, nrow = S, ncol = TIME)
    } else {
      N.out <- exp(outvect[1:S])
      cvec.out <- outvect[(S+1):(S+qp+1)]
      par.index <- S + qp + 1
      p.temp <- matrix(cvec.out[1], nrow = S, ncol = TIME)
      for(q in 1:qp) p.temp <- p.temp + cvec.out[q+1]*cov.p[,,q]
      p.out <- expo(p.temp)
      #c0.out <- outvect[S+1]; 
      #c1.out <- outvect[S+2]; 
      #par.index <- S + 2
      #p.out <- matrix(expo(c0.out + c1.out*cov.p), byrow = F, nrow = S, ncol = TIME)
    }
    d0.out <- NULL
    d1.out <- NULL
    if (M>1) {
      d0.out <- outvect[(par.index+1):(par.index+M-1)];
      par.index <- par.index + M - 1
      d1.out <- NULL
      w.out <- matrix(0,S,M)
      if(w.m == "common") { 
        for(i in 1:S) w.out[i,] <- lgp(d0.out)
      }
      if(w.m == "cov") { 
        d1.out <- outvect[par.index+1]; 
        for(i in 1:S) w.out[i,] <- lgp(d0.out + d1.out * cov.w[i]); 
        par.index <- par.index + 1
        if(M==2){
          for(i in 1:S){
            form <- sprintf("~x1 + x2*%f",cov.w[i])
          }
        }
      }
    } else {
      w.out <- matrix(1,S,M)
    }
    
    mu.out <- matrix(NA,S,M)
    
    if(mu.m == "common") {
      b0.out <- outvect[(par.index + 1):(par.index + M)]; 
      b1.out <- NULL
      for(m in 1:M) {
        mu.out[,m] <- exp(b0.out[m])
      }
      
      par.index <- par.index + M
    }
    
    if(mu.m == "cov") {
      b0.out <- outvect[(par.index + 1):(par.index + M)]; 
      b1.out <- outvect[(par.index + M + 1)]
      for(i in 1:S) {
        for(m in 1:M) {
          mu.out[i,m] <- exp(b0.out[m] + b1.out*cov.mu[i])
          form <- sprintf("~x1 + x2*%f",cov.mu[i])
        }
      }
      par.index <- par.index + M + 1
    }
    
    sigma.out <- exp(outvect[par.index + 1]); 
    par.index <- par.index + 1;
    if(sigma.m == "het") { 
      sigma.out <- c(sigma.out, exp(outvect[(par.index + 1):(par.index + M - 1)])); 
      par.index <- par.index + M - 1
    }
    
    sigma.out <- matrix(sigma.out, S, M, byrow = TRUE)
    
    betta.out <- matrix(0,S,TIME)
    
    for(i in 1:S) {
      for(m in 1:M) {
        betta.out[i,] <- betta.out[i,]+c(w.out[i,m]*c(pnorm(1,mean=mu.out[i,m],sd=sigma.out[i,m]),pnorm(2:(TIME-1),mean=mu.out[i,m],sd=sigma.out[i,m])-pnorm(1:(TIME-2),mean=mu.out[i,m],sd=sigma.out[i,m]),1-pnorm(TIME-1,mean=mu.out[i,m],sd=sigma.out[i,m])))
      }
    }
    
    phi.out <- array(0,c(TIME-1,TIME-1,S))
    
    if(phi.m == "const") {
      a0.out <- outvect[par.index + 1]; a1.out <- NULL; a2.out <- NULL
      for(i in 1:S) {	
        phi.out[,,i] <- expo(a0.out)
        phi.out[,,i][lower.tri(phi.out[,,i])]<-0
      }
    }
    
    if(phi.m == "logit.t") {
      a0.out <- outvect[par.index + 1]; 
      a1.out <- outvect[par.index + 2]; 
      a2.out <- NULL; par.index <- par.index + 2
      
      for(i in 1:S) {	
        for(j in 1:(TIME-1)) {
          form <- sprintf("~x1 + x2*%f",cov.phi[i,j])
          phi.out[,j,i] <- expo(a0.out + a1.out*cov.phi[i,j])
        }
        phi.out[,,i][lower.tri(phi.out[,,i])]<-0
      }
    }
    
    if(phi.m == "logit.a") {
      a0.out <- outvect[par.index + 1]; 
      a1.out <- outvect[par.index + 2]; 
      a2.out <- NULL; par.index <- par.index + 2
      
      for(i in 1:S) {
        for(b in 1:(TIME-1)){
          for(j in b:(TIME-1)){
            form <- sprintf("~x1 + x2*%f",cov.phi[i,j-b+1])
            phi.out[b,j,i] <- expo(a0.out + a1.out*cov.phi[i,j-b+1])
          }
        }
        phi.out[,,i][lower.tri(phi.out[,,i])]<-0
      }
    }
    
    if(phi.m == "quadr.t") {
      a0.out <- outvect[par.index + 1]; 
      a1.out <- outvect[par.index + 2]; 
      a2.out <- outvect[par.index + 3]; 
      par.index <- par.index + 3
      
      for(i in 1:S) {	
        for(j in 1:(TIME-1)) {
          form <- sprintf("~x1 + x2*%f + x3*%g",cov.phi[i,j],cov.phi[i,j]^2)
          phi.out[,j,i] <- expo(a0.out + a1.out*cov.phi[i,j] + a2.out*cov.phi[i,j]^2)
        }
        phi.out[,,i][lower.tri(phi.out[,,i])]<-0
      }
    }
    
    if(phi.m == "quadr.a") {
      a0.out <- outvect[par.index + 1]; 
      a1.out <- outvect[par.index + 2]; 
      a2.out <- outvect[par.index + 3]; 
      par.index <- par.index + 3
      
      for(i in 1:S) {
        for(b in 1:(TIME-1)) {
          for(j in b:(TIME-1)) {
            form <- sprintf("~x1 + x2*%f + x3*%g",cov.phi[i,j-b+1],cov.phi[i,j-b+1]^2)
            phi.out[b,j,i] <- expo(a0.out + a1.out*cov.phi[i,j-b+1] + a2.out*cov.phi[i,j-b+1]^2)
          }
        }
        phi.out[,,i][lower.tri(phi.out[,,i])]<-0
      }
    }
    
    list("ll.val"=-this.fit$value,"npar"=length(this.fit$par),"N.est"=N.out,"cest" = cvec.out,"p.est"=p.out,"d0.est"=d0.out,"d1.est"=d1.out,"b0.est"=b0.out,"b1.est"=b1.out,"w.est"=w.out,"mu.est"=mu.out,"sigma.est"=sigma.out,"betta.est"=betta.out,"inter.est"=a0.out,"slope.est"=a1.out,"slope2.est"=a2.out,"phi.est"=phi.out,"Hessian"=this.fit$hessian) 
    
  }else{
    list("ll.val"=NA)
  }
}