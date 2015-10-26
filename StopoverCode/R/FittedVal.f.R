#' @export
#------------------- Function to obtain fitted values ---------------------------------#

FittedVal.f <- function(model){
  S <- length(model$N.est)
  t <- dim(model$betta.est)[2]
  FittedVal <- matrix(0, S, t)
  for(i in 1:S){
    FittedVal[i,] <- model$N.est[i]*model$betta.est[i,]*model$p.est[i,]
    for(j in 2:t){
      for(b in 1:(j-1)){
        FittedVal[i,j] <- FittedVal[i,j] + model$N.est[i]*model$betta.est[i,b]*prod(model$phi.est[b,b:(j-1),i])*model$p.est[i,j]
      }
    }
  }
  FittedVal
}