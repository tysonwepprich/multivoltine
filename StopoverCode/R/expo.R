#' @export
#-------------------  EXPO FUNCTION  ---------------------------

# Inverse of logit function, acts on a vector or matrix.

expo <- function(inval)
  1/(1+exp(-inval))

