#' @export
#------------------- Function to calculate the standardised deviance residuals of fitted model ---------------------------------#

StDevResid.f <- function(FittedVal)
{
  
  StDevResid <- matrix(0, S, T)
  
  for(i in 1:S){
    for(j in 2:T){
      if(counts[i,j] > -1) {if(counts[i,j]==0) counts[i,j] <- 0.0001; StDevResid[i,j] <- sign(counts[i,j]-FittedVal[i,j])*sqrt(2*(counts[i,j]*log(counts[i,j]/FittedVal[i,j])-(counts[i,j]-FittedVal[i,j])))
      }}
  }
  StDevResid
  
}	
