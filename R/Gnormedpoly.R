########################################################################################
########################################################################################
# Polytomous items: Number of Guttman errors
# Gpoly reduces to G for 0-1 data (Ncat=2) (currently Gpoly needs G to be loaded)
# Gnormed.poly reduces to Gnormed for 0-1 data (Ncat=2) (currently Gnormed.poly needs G to be loaded)
########################################################################################
########################################################################################

Gnormed.poly <- function(matrix, Ncat,
                          NA.method="NPModel", Save.MatImp=FALSE, 
                          IP=NULL, IRT.PModel="GRM", Ability=NULL, Ability.PModel="EAP")
{
  matrix      <- as.matrix(matrix)
  N           <- dim(matrix)[1]; I <- dim(matrix)[2]; M <- Ncat-1
  IP.NA       <- is.null(IP); Ability.NA  <- is.null(Ability)
  # Sanity check - Data matrix adequacy:
  Sanity.dma.poly(matrix, N, I, M)
  # Dealing with missing values:
  res.NA <- MissingValues.poly(matrix, Ncat, NA.method, Save.MatImp, IP, IRT.PModel, Ability, Ability.PModel)
  matrix <- res.NA[[1]]
  # Perfect response vectors allowed (albeit uninformative).
  # Compute PFS:
  # Numerator:
  probs.ISD      <- matrix(NA, nrow=I, ncol=M)
  for (m in 1:M) {probs.ISD[,m] <- colMeans(matrix >= m)}
  f.scoresISD    <- function (x) {c(rep(1,x), rep(0,M-x))}
  matrix.ISD     <- matrix(unlist(lapply(t(matrix), f.scoresISD)), byrow=TRUE, nrow=N)
  probs.ISD.vect <- as.vector(t(probs.ISD))
  matrix.ISD.ord <- matrix.ISD[, order(probs.ISD.vect, decreasing=TRUE)]
  num            <- G(matrix.ISD.ord)$PFscores[,1]
  # Denominator: 
  NC        <- rowSums(matrix)
  ranks.ISD <- matrix(rank(I*M - probs.ISD.vect, ties.method="first"), nrow=I, byrow=TRUE)
  if (Ncat>2) {cumranks.ISD <- cbind(rep(0,I), t(apply(ranks.ISD, 1, cumsum)))}
  if (Ncat==2) {cumranks.ISD <- cbind(rep(0,I), ranks.ISD)}
  V        <- matrix(rep(cumranks.ISD[1,], Ncat), ncol=Ncat, byrow=FALSE)
  add.term <- matrix(rep(cumranks.ISD[2,], nrow(V)), ncol=Ncat, byrow=TRUE)
  V        <- V + add.term
  T        <- sapply(2:sum(dim(V)), function(x) {max(V[col(V)+row(V) == x])})
  for (i in 3:I) 
  {
    V        <- matrix(rep(T, Ncat), ncol=Ncat, byrow=FALSE)
    add.term <- matrix(rep(cumranks.ISD[i,], nrow(V)), ncol=Ncat, byrow=TRUE)
    V        <- V + add.term
    T        <- sapply(2:sum(dim(V)), function(x) {max(V[col(V)+row(V) == x])})
  }
  maxG       <- T - sapply(0:(I*M), function(x) {.5*x*(x+1)})
  maxG[c(1,length(maxG))] <- 1 # so that vectors (0,0,...,0) and (M,M,...,M), which have 0 Guttman errors, are divided by 1 instead of 0
  # 
  res <- num / maxG[NC+1]
  # Export results:
  export.res.NP(matrix, N, res, "Gnormed.poly", vector("list", 5) , Ncat=Ncat, NA.method, 
               IRT.PModel, res.NA[[2]], Ability.PModel, res.NA[[3]], IP.NA, Ability.NA, res.NA[[4]])

}