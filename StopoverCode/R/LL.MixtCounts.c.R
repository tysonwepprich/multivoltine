#' @useDynLib StopoverCode MixtCountsLL
#' @export
### The likelihood ###

LL.MixtCounts.c <- function(N.val, p.mat, beta.mat, phi.arr)
{
  .C("MixtCountsLL",
     as.integer(counts),
     as.integer(S),
     as.integer(TIME),
     as.double(N.val),
     as.double(p.mat),
     as.double(beta.mat),
     as.double(phi.arr),
     LL = double (1),
     resultcode = integer(1),
     lambda = double (1)
  )$LL
}