########################################################################################
########################################################################################
# U3 (van der Flier, 1980, 1982):
########################################################################################
########################################################################################
U3 <- function(matrix, 
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
  pi <- colMeans(matrix.sv); qi <- 1-pi
  # If there are answer options not chosen by any respondent then some entries in pi are 0 or 1.
  # Below all corresponding logs are set from Inf to 0.
  # (Reason: They carry no individual information regarding aberrant response behavior.):
  log.odds     <- log(pi/qi)
  log.odds[is.infinite(log.odds)] <- 0
  log.odds.ord <- sort(log.odds,decreasing=TRUE)
  # 
  sum.first.logodds <- cumsum(log.odds.ord)[NC]
  logodds.ordrev    <- sort(log.odds, decreasing=FALSE)
  sum.last.logodds  <- cumsum(logodds.ordrev)[NC]
  res.red           <- as.vector((sum.first.logodds - as.vector(matrix %*% log.odds)) / (sum.first.logodds - sum.last.logodds))
  # Compute final PFS vector:
  res <- final.PFS(res.red, all.0s, all.1s, N)
  # Export results:
  export.res.NP(matrix.sv, N, res, "U3", part.res, Ncat=2, NA.method, 
                IRT.PModel, res.NA[[2]], Ability.PModel, res.NA[[3]], IP.NA, Ability.NA, res.NA[[4]])
}
