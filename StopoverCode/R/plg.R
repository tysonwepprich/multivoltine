#' @export
#-------------------  PLG FUNCTION  ---------------------------

# Proportions to logits (one less element)
# Use for transforming beta to logit gamma.
# Slight increment in cumsum, to stop NaN if beta_j  = 0.

plg <- function(invec)
{
  # How many elements?
  nel <- length(invec)
  # Case of 2 els:
  if (nel==2) outvec <- invec[1]
  # Case of more than 2 els:
  if (nel>2)
  {
    # Cumulative sum:
    cs <- cumsum(invec)+0.000000001 # Stop NaNs
    # Output vector:
    outvec <- invec[1:(nel-1)]/(1-c(0,cs[1:(nel-2)]))
  }
  # return logit:
  logit(outvec)
}