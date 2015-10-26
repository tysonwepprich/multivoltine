#' @export
#------------------- Function to calculate the residual deviance of fitted model ---------------------------------#

ResDev.f <- function(FittedVal)
{
  
  ResDev <- 0
  for(i in 1:S){
    for(j in 2:T){
      if(counts[i,j] > -1) {if(counts[i,j]==0) counts[i,j] <- 0.0001; ResDev <- ResDev + 2*(counts[i,j]*log(counts[i,j]/FittedVal[i,j])-(counts[i,j]-FittedVal[i,j]))
      }}
  }
  ResDev
  
}		