########################################################################################
########################################################################################
# Polytomous items: U3p
# U3poly reduces to U3 for 0-1 data (Ncat=2)
########################################################################################
########################################################################################
U3poly <- function(matrix, Ncat,
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
  NC        <- rowSums(matrix)
  probs.ISD <- matrix(NA,nrow=I,ncol=M)
  for (m in 1:M)
  {
    probs.ISD[,m] <- colMeans(matrix >= m)
  }
  f.scoresISD    <- function (x) {c(rep(1,x), rep(0,M-x))}
  matrix.ISD     <- matrix(unlist(lapply(t(matrix), f.scoresISD)), byrow=TRUE, nrow=N)
  probs.ISD.vect <- as.vector(t(probs.ISD))
  # If there are answer options not chosen by any respondent then some entries in 'probs.ISD.vect' are 0 and others are 1.
  # Below all corresponding logs are set from Inf to 0.
  # (Reason: They carry no information regarding aberrant response behavior).
  logits.ISD <- log(probs.ISD.vect/(1-probs.ISD.vect))
  logits.ISD[is.infinite(logits.ISD)] <- 0 # a vector
  W <- as.vector(matrix.ISD %*% logits.ISD)
  # 
  logits.ISD <- matrix(logits.ISD, nrow=I, byrow=TRUE) # a matrix
  if (Ncat>2)  {cumlogits.ISD <- cbind(rep(0,I), t(apply(logits.ISD, 1, cumsum)))}
  if (Ncat==2) {cumlogits.ISD <- cbind(rep(0,I), logits.ISD)}
  V        <- matrix(rep(cumlogits.ISD[1,], Ncat), ncol=Ncat, byrow=FALSE)
  add.term <- matrix(rep(cumlogits.ISD[2,], nrow(V)), ncol=Ncat, byrow=TRUE)
  V        <- V + add.term
  T        <- sapply(2:sum(dim(V)),function(x) {min(V[col(V)+row(V) == x])})
  for (i in 3:I) {
    V        <- matrix(rep(T,Ncat), ncol=Ncat, byrow=FALSE)
    add.term <- matrix(rep(cumlogits.ISD[i,], nrow(V)), ncol=Ncat, byrow=TRUE)
    V        <- V + add.term
    T        <- sapply(2:sum(dim(V)),function(x) {min(V[col(V)+row(V) == x])})
  }
  maxW <- sapply(0:(I*M), function(x) {sum(sort(as.vector(logits.ISD),decreasing=TRUE)[0:x])})
  minW <- T
  den  <- maxW - minW
  den[c(1,length(den))] <- 1 # so that vectors (0,0,...,0) and (M,M,...,M), which have 0 Guttman errors, are divided by 1 instead of 0
  #   
  res <- (maxW[NC+1] - W) / den[NC+1]
  # Export results:
  export.res.NP(matrix, N, res, "U3poly", vector("list", 5) , Ncat=Ncat, NA.method, 
               IRT.PModel, res.NA[[2]], Ability.PModel, res.NA[[3]], IP.NA, Ability.NA, res.NA[[4]])
}
