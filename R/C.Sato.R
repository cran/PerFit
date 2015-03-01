########################################################################################
########################################################################################
# C (Sato, 1975):
########################################################################################
########################################################################################
C.Sato <- function(matrix, 
                   NA.method="NPModel", Save.MatImp=FALSE, 
                   IP=NULL, IRT.PModel="2PL", Ability=NULL, Ability.PModel="ML", mu=0, sigma=1)
{
  matrix      <- as.matrix(matrix)
  N           <- dim(matrix)[1]; I <- dim(matrix)[2]
  IP.NA       <- is.null(IP); Ability.NA  <- is.null(Ability)
  # Sanity check - Data matrix adequacy:
  Sanity.dma(matrix, N, I)
  # Dealing with missing values:
  res.NA <- MissingValues(matrix, NA.method, Save.MatImp, IP, IRT.PModel, Ability, Ability.PModel, mu, sigma)
  matrix <- res.NA[[1]]
  # Sanity check - Perfect response vectors:
  part.res  <- Sanity.prv(matrix, N, I)
  NC        <- part.res$NC
  all.0s    <- part.res$all.0s
  all.1s    <- part.res$all.1s
  matrix.sv <- matrix
  matrix    <- part.res$matrix.red
  # Compute PFS:
  pi             <- colMeans(matrix.sv)
  pi.ord         <- sort(pi, decreasing=TRUE)
  matrix.ord     <- matrix[, order(pi, decreasing=TRUE)]
  num            <- cov(t(matrix.ord), pi.ord)
  matrix.easiest <- (col(matrix) <= NC) * 1
  den            <- cov(t(matrix.easiest), pi.ord)
  res.red        <- as.vector(1 - num / den)
  # Compute final PFS vector:
  res <- final.PFS(res.red, all.0s, all.1s, N)
  # Export results:
  export.res.NP(matrix.sv, N, res, "C.Sato", part.res, Ncat=2, NA.method, 
                IRT.PModel, res.NA[[2]], Ability.PModel, res.NA[[3]], IP.NA, Ability.NA, res.NA[[4]])
}
