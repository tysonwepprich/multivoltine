#' @export
#-------------------  LOGIT FUNCTION  ---------------------------

# If p is really close to 0 or 1, set logit to some max or min.

logit <- function(inval)  # Can be a single value, vector or matrix
{
  low.val  <- 0.00001
  high.val <- 0.99999
  inval[inval<low.val] <- low.val
  inval[inval>high.val] <- high.val
  # Output vector:
  log(inval/(1-inval))
}