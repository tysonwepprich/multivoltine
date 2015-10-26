#' @export
#-------------------  LGP FUNCTION  ---------------------------

# Logit to proportions (one extra element)
# Use for logit gamma to beta.
# Works for vector of length 1 or more.

lgp <- function(invec)
{
  # How many elements?
  nel <- length(invec)
  # Case of 1 el:
  if (nel==1) p.v <- expo(invec)
  # Case of more than 1 el:
  if (nel>1)
  {
    # Get inverse-logit vector:
    thisvec <- expo(invec)
    # Convert to proportions vector: need cum product of 1-thisvec
    cp <- cumprod(1-thisvec)
    # Return proportions vector:
    p.v <- thisvec*c(1,cp[1:(nel-1)])
  }
  c(p.v,1-sum(p.v))
}
