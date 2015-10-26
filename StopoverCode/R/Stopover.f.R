#' @export
#------------------- Function to calculate mean stopover at each site ---------------------------------#

Stopover.f <- function(model)
{
  
  Stopover <- rep(0, S)
  for(i in 1:S){
    for(b in 1:TIME){
      for(d in b:TIME){
        if((d>b)&(d<TIME)) Stopover[i] <- Stopover[i] + (d-b+1)*model$betta.est[i,b]*prod(model$phi.est[b,b:(d-1),i])*(1-model$phi.est[b,d,i])
        if((d>b)&(d==TIME)) Stopover[i] <- Stopover[i] + (d-b+1)*model$betta.est[i,b]*prod(model$phi.est[b,b:(d-1),i])
        if((d==b)&(d<TIME)) Stopover[i] <- Stopover[i] + (d-b+1)*model$betta.est[i,b]*(1-model$phi.est[b,d,i])
        if(b==TIME) Stopover[i] <- Stopover[i] + (d-b+1)*model$betta.est[i,b]
      }
    }
  }
  Stopover
  
}								
