#' @export
#-------- Function to generate starting values -----------#
startVals <- function(p.m = "common", w.m = "common", mu.m = "common", sigma.m = "hom", phi.m = "const")
{
  
  N.st <- NULL
  cvec.st <- NULL #c0.st <- NULL; c1.st <- NULL
  
  if(p.m == "common"){ 
    N.st <- log(apply(counts,1,max,na.rm=TRUE)+50)	
  } else { 
    N.st <- log(apply(counts,1,max,na.rm=TRUE)+50)
    cvec.st <- rep(0, qp+1) #c0.st <- 0; c1.st <- 1
  }
  
  d0.st <- NULL; if(M>1) d0.st <- plg(rep(1/M, M))    #relative weights
  d1.st <- NULL
  if(w.m == "cov") d1.st <- 0 #slope for weight when covariate is available
  
  b0.st <- sort(log(sample(1:TIME,M,F)))   #mean arrival times
  b1.st <- NULL
  if(mu.m == "cov") b1.st <- 0 #slope for mu when covariate is available 
  
  sigma.st <- log(6)  #standard deviations of arrival times
  if(sigma.m == "het") sigma.st <- c(sigma.st, rep(log(6), M-1))
  
  a0.st <- logit(sample(seq(.3,.7,.1),1))	 #survival probabilities
  a1.st <- NULL; a2.st <- NULL
  if((phi.m == "logit.t")|(phi.m == "logit.a")) a1.st <- 0
  if((phi.m == "quadr.t")|(phi.m == "quadr.a")){a1.st <- 0; a2.st <- 0}
  
  list("N" = N.st, "cvec" = cvec.st,  "d0" = d0.st, "d1" = d1.st, "b0" = b0.st, "b1" = b1.st, "sigma" = sigma.st, "a0" = a0.st, "a1" = a1.st, "a2" = a2.st) 
}
