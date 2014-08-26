library(irtoys)

# Determine cutoff for class "PerFit":
cutoff <- function(x,method="Quantile",Qlvl=.05,Blvl=.05,Breps=1000,UDlvl=NA) { #x = an object from 'PerFit' class
  upp.PFS <- c("Cstar","C.Sato","U3","ZU3","G","Gnormed","Gpoly","Gnormed.poly","U3poly","D.KB");
  low.PFS <- c("r.pbis","NCI","Ht","A.KB","E.KB","lz","lzstar","lzpoly");
  if (method == "Quantile") {
    prop.flagged <- Qlvl;
    if (any(x$PFStatistic == upp.PFS)) {
      tail <- "upper";
      Qlvl.use <- 1-Qlvl}
    if (any(x$PFStatistic == low.PFS)) {
      tail <- "lower";
      Qlvl.use <- Qlvl}    
    cutoff.use <- quantile(x$PFscores,probs=Qlvl.use);
  }
  #
  if (method == "Bootstrap") {
    if (any(x$PFStatistic == upp.PFS)) {
      tail <- "upper";
      Blvl.use <- 1-Blvl}
    if (any(x$PFStatistic == low.PFS)) {
      tail <- "lower";
      Blvl.use <- Blvl}   
    Bvec <- c();
    for (i in 1:Breps){
      Bvec <- c(Bvec,quantile(sample(x$PFscores,size=length(x$PFscores),replace=TRUE),probs=Blvl.use))
    }
    cutoff.use <- median(Bvec);
    if (any(x$PFStatistic == upp.PFS)) {prop.flagged <- sum(x$PFscores >= cutoff.use) / length(x$PFscores)}
    if (any(x$PFStatistic == low.PFS)) {prop.flagged <- sum(x$PFscores <= cutoff.use) / length(x$PFscores)}
  }
  #
  if (method == "UserDefined") {
    cutoff.use <- UDlvl;
    if (any(x$PFStatistic == upp.PFS)) {
      prop.flagged <- sum(x$PFscores >= cutoff.use) / length(x$PFscores);
      tail <- "upper"}
    if (any(x$PFStatistic == low.PFS)) {
      prop.flagged <- sum(x$PFscores <= cutoff.use) / length(x$PFscores);
      tail <- "lower"}
  }
  #
  list(cutoff=as.numeric(cutoff.use),prop.flagged=prop.flagged,tail=tail)
}

# Define plot() function for class "PerFit":
plot.PerFit <- function (x,type="Density",both.scale=TRUE,cutoff=TRUE,method="Quantile",Qlvl=.05,Blvl=.05,Breps=1000,UDlvl=NA,
                         Xlabel=NA,Xcex=1.5,title=NA,Tcex=1.5,...) { #x = an object from 'PerFit' class
  upp.PFS <- c("Cstar","C.Sato","U3","ZU3","G","Gnormed","Gpoly","Gnormed.poly","U3poly","D.KB");
  low.PFS <- c("r.pbis","NCI","Ht","A.KB","E.KB","lz","lzstar","lzpoly");
  cutoff.res <- cutoff(x,method,Qlvl,Blvl,Breps,UDlvl);
  x.line <- cutoff.res$cutoff;
  perc.flagged <- round(100 * cutoff.res$prop.flagged,2);
  direction <- paste(", ",cutoff.res$tail," ",sep="");
  # Find correct scale for y-axis:
  ymax.hist <- max(hist(x$PFscores,plot=FALSE)$density);
  ymax.dens <- max(density(x$PFscores)$y);
  ymax <- switch(type,
                 Density=ymax.dens,
                 Histogram=ymax.hist,
                 Both=if (both.scale == TRUE) {max(ymax.dens,ymax.hist)} else {min(ymax.dens,ymax.hist)})
  par(mar=c(4,3.5,2,1)+.1,las=1);
  hist(x$PFscores,freq=FALSE,border="white",ann=FALSE,ylim=c(0,ymax));
  #
  if (cutoff == TRUE) {
    if (any(x$PFStatistic == upp.PFS)) {rect(x.line,0,par("usr")[2], par("usr")[4],col="lightpink1",border=NA)}
    if (any(x$PFStatistic == low.PFS)) {rect(par("usr")[1],0, x.line,par("usr")[4],col="lightpink1",border=NA)}
  }
  #
  if (type == "Histogram") {
    par(new=TRUE)
    hist(x$PFscores,freq=FALSE,col="lightblue",ann=FALSE,ylim=c(0,ymax))}
  #
  if (type == "Density") {
    points(density(x$PFscores),type="l",lwd=2,ann=FALSE,ylim=c(0,ymax))}
  #
  if (type == "Both") {
    par(new=TRUE)
    hist(x$PFscores,freq=FALSE,col="lightblue",ann=FALSE,ylim=c(0,ymax))
    points(density(x$PFscores),type="l",lwd=2,ann=FALSE,ylim=c(0,ymax))}
  box(col="black")
  #
  if (cutoff == FALSE) {
    tmp <- if (is.na(Xlabel)) {x$PFStatistic} else {Xlabel};
    mtext(side=1,text=tmp,line=2.5,col="black",cex=1.5,font=1)}
  #
  if (cutoff == TRUE) {
    abline(v=x.line,lwd=2)
    tmp <- if (is.na(Xlabel)) {paste(x$PFStatistic," (cutoff=",round(x.line,3),direction,perc.flagged,"%)",sep="")} else {Xlabel};
    mtext(side=1,text=tmp,line=2.5,col="black",cex=Xcex,font=1) 
  }
  #
  tmp <- if (is.na(title)) {"Distribution"} else {title};
  mtext(side=3,text=tmp,line=.5,col="black",cex=Tcex,font=2)
}

flagged.resp <- function(matrix,scores=TRUE,ord=TRUE,x,method="Quantile",Qlvl=.05,Blvl=.05,Breps=1000,UDlvl=NA) { # matrix = score matrix
  upp.PFS <- c("Cstar","C.Sato","U3","ZU3","G","Gnormed","Gpoly","Gnormed.poly","U3poly","D.KB");
  low.PFS <- c("r.pbis","NCI","Ht","A.KB","E.KB","lz","lzstar","lzpoly");
  if (any(x$PFStatistic == upp.PFS)) {flagged.subs <- which(x$PFscores >= cutoff(x,method,Qlvl,Blvl,Breps,UDlvl)$cutoff)};
  if (any(x$PFStatistic == low.PFS)) {flagged.subs <- which(x$PFscores <= cutoff(x,method,Qlvl,Blvl,Breps,UDlvl)$cutoff)};
  Ps <- round(apply(matrix,2,mean),3);
  # Not ordered by pvalue:
  if (ord == F) {
    flagged.scores <- matrix[flagged.subs,];
    colnames(flagged.scores) <- paste("It",1:dim(matrix)[2],sep="");
  }
  # Ordered by pvalue:
  if (ord == T) {
    matrix.ord <- matrix[,order(Ps,decreasing=TRUE)]; # ordered from easy to difficult
    flagged.scores <- matrix.ord[flagged.subs,];
    colnames(flagged.scores) <- paste("It",order(Ps,decreasing=TRUE),sep="");
    Ps <- sort(Ps,decreasing=TRUE);
  }
  #rownames(flagged.scores) <- NULL; #paste("Resp. ",flagged.subs,sep="");
  flagged.scores <- as.matrix(flagged.scores);
  rownames(flagged.scores) <- NULL;
  res <- if (scores == FALSE) {
    cbind(FlaggedID=flagged.subs,PFscores=x$PFscores[flagged.subs])} else {
      list(Scores=cbind(FlaggedID=flagged.subs,flagged.scores,PFscores=x$PFscores[flagged.subs]),MeanItemValue=Ps)}
  res
}

########################################################################################
########################################################################################
# C* (Harnisch & Linn, 1981):
########################################################################################
########################################################################################
Cstar <- function(matrix){
  matrix <- as.matrix(matrix);
  N <- dim(matrix)[1]; I <- dim(matrix)[2];
  NC <- apply(matrix,1,sum);
  uniqueNC <- sort(unique(NC));
  if (min(uniqueNC)==0 || max(uniqueNC)==I){
    stop('Remove all-0s and/or all-1s response vectors from data.')};
  pi <- apply(matrix,2,mean);
  pi.ord <- sort(pi,decreasing=TRUE);
  res <- rep(NA,N);
  for (i in 1:length(uniqueNC)){
    Respond.i <- NC==uniqueNC[i];
    matrix.i <- matrix(matrix[Respond.i,],nrow=sum(Respond.i),byrow=FALSE);
    num <- as.vector(sum(pi.ord[1:uniqueNC[i]]) - matrix.i %*% pi);
    den <- sum(pi.ord[1:uniqueNC[i]]) - sum(pi.ord[(I-uniqueNC[i]+1):I]);
    res[Respond.i] <- num/den;
  }
  res <- list(PFscores=round(res,4),PFStatistic="Cstar");
  class(res) <- "PerFit";
  res
}

########################################################################################
########################################################################################
# C (Sato, 1975):
########################################################################################
########################################################################################
C.Sato <- function(matrix){
  matrix <- as.matrix(matrix);
  N <- dim(matrix)[1]; I <- dim(matrix)[2];
  NC <- apply(matrix,1,sum);
  uniqueNC <- sort(unique(NC));
  if (min(uniqueNC)==0 || max(uniqueNC)==I){
    stop('Remove all-0s and/or all-1s response vectors from data.')};
  pi <- apply(matrix,2,mean);
  pi.ord <- sort(pi,decreasing=TRUE);
  res <- rep(NA,N);
  for (i in 1:length(uniqueNC)){
    Respond.i <- NC==uniqueNC[i];
    matrix.i <- matrix(matrix[Respond.i,],nrow=sum(Respond.i),byrow=FALSE);
    matrix.iORD <- matrix(matrix.i[,order(pi,decreasing=TRUE)],nrow=sum(Respond.i),byrow=FALSE);
    res[Respond.i] <- 1-cov(t(matrix.iORD),pi.ord) / cov(c(rep(1,uniqueNC[i]),rep(0,I-uniqueNC[i])),pi.ord)
  }
  res <- list(PFscores=round(res,4),PFStatistic="C.Sato");
  class(res) <- "PerFit";
  res
}

########################################################################################
########################################################################################
# r.pbis, personal point-biserial correlation (Brennan, 1980, cited in Harhisch & Linn, 1981):
########################################################################################
########################################################################################
r.pbis <- function(matrix){
  N <- dim(matrix)[1]; I <- dim(matrix)[2];
  NC <- apply(matrix,1,sum);
  uniqueNC <- sort(unique(NC));
  if (min(uniqueNC)==0 || max(uniqueNC)==I){
    stop('Remove all-0s and/or all-1s response vectors from data.')};
  pi <- apply(matrix,2,mean);
  ri <- list(PFscores=as.vector(round(apply(matrix,1,function(vect){cor(vect,pi)}),4)),PFStatistic="r.pbis");
  class(ri) <- "PerFit";
  ri
}

########################################################################################
########################################################################################
# U3 (van der Flier, 1980, 1982):
########################################################################################
########################################################################################
U3 <- function(matrix){
  matrix <- as.matrix(matrix);
  N <- dim(matrix)[1]; I <- dim(matrix)[2];
  NC <- apply(matrix,1,sum);
  uniqueNC <- sort(unique(NC));
  if (min(uniqueNC)==0 || max(uniqueNC)==I){
    stop('Remove all-0s and/or all-1s response vectors from data.')};
  pi <- apply(matrix,2,mean); qi <- 1-pi;
  # If there are answer options not chosen by any respondent then some entries in pi are 0 or 1.
  # Below all corresponding logs are set from Inf to 0 
  # (reason: They carry no individual information regarding aberrant response behavior):
  log.odds <- log(pi/qi); log.odds <- sapply(log.odds,function(x){if (is.infinite(x)) 0 else x});
  log.odds.ord <- sort(log.odds,decreasing=TRUE);
  res <- rep(NA,N);
  for (i in 1:length(uniqueNC)){
    Respond.i <- NC==uniqueNC[i];
    matrix.i <- matrix(matrix[Respond.i,],nrow=sum(Respond.i),byrow=FALSE);
    num.U3 <- sum(log.odds.ord[1:uniqueNC[i]]) - matrix.i %*% log.odds;
    den.U3 <- sum(log.odds.ord[1:uniqueNC[i]]) - sum(log.odds.ord[(I-uniqueNC[i]+1):I]);
    res[Respond.i] <- num.U3/den.U3;  
  }
  res <- list(PFscores=round(res,4),PFStatistic="U3");
  class(res) <- "PerFit";
  res
}

########################################################################################
########################################################################################
# ZU3 (van der Flier, 1980, 1982):
########################################################################################
########################################################################################
ZU3 <- function(matrix){
  matrix <- as.matrix(matrix);
  N <- dim(matrix)[1]; I <- dim(matrix)[2];
  NC <- apply(matrix,1,sum);
  uniqueNC <- sort(unique(NC));
  if (min(uniqueNC)==0 || max(uniqueNC)==I){
    stop('Remove all-0s and/or all-1s response vectors from data.')};
  pi <- apply(matrix,2,mean); qi <- 1-pi;
  # If there are answer options not chosen by any respondent then some entries in pi are 0 or 1.
  # Below all corresponding logs are set from Inf to 0 (reason: They carry no information regarding aberrant response behavior):
  log.odds <- log(pi/qi); log.odds <- sapply(log.odds,function(x){if (is.infinite(x)) 0 else x});
  log.odds.ord <- sort(log.odds,decreasing=TRUE);
  res <- rep(NA,N);
  for (i in 1:length(uniqueNC)){
    Respond.i <- NC==uniqueNC[i];
    matrix.i <- matrix(matrix[Respond.i,],nrow=sum(Respond.i),byrow=FALSE);
    alpha <- sum(pi*log.odds) + sum(pi*qi*log.odds)*(uniqueNC[i]-sum(pi)) / sum(pi*qi);
    beta <- sum(pi*qi*(log.odds)^2) - (sum(pi*qi*log.odds))^2/sum(pi*qi);
    exp.val <- (sum(log.odds.ord[1:uniqueNC[i]]) - alpha) / (sum(log.odds.ord[1:uniqueNC[i]]) - sum(log.odds.ord[(I-uniqueNC[i]+1):I]));
    var.val <- beta / ((sum(log.odds.ord[1:uniqueNC[i]]) - sum(log.odds.ord[(I-uniqueNC[i]+1):I]))^2);
    num.U3 <- sum(log.odds.ord[1:uniqueNC[i]]) - matrix.i %*% log.odds;
    den.U3 <- sum(log.odds.ord[1:uniqueNC[i]]) - sum(log.odds.ord[(I-uniqueNC[i]+1):I]);
    res[Respond.i] <- ((num.U3/den.U3) - exp.val) / sqrt(var.val);
  }
  res <- list(PFscores=round(res,4),PFStatistic="ZU3");
  class(res) <- "PerFit";
  res
}

########################################################################################
########################################################################################
# G (van der Flier, 1977; Meijer, 1994):
########################################################################################
########################################################################################
G <- function(matrix){
  N <- dim(matrix)[1]; I <- dim(matrix)[2];
  NC <- apply(matrix,1,sum);
  uniqueNC <- sort(unique(NC));
  pi <- apply(matrix,2,mean);
  matrix.ord <- matrix[,order(pi,decreasing=TRUE)];
  per.row <- function(vect){
    NC <- sum(vect);
    sum <- 0;
    for (i in 1:(I-1)){sum <- sum + sum(diff(vect,lag=i)==1)}
    sum;
  }
  res <- apply(matrix.ord,1,per.row);
  res <- list(PFscores=as.vector(round(res,4)),PFStatistic="G");
  class(res) <- "PerFit";
  res
}

########################################################################################
########################################################################################
# Gnormed (van der Flier, 1977; Meijer, 1994):
# This statistic is perfectly linearly related to NCI (Tatsuoka & Tatsuoaka, 1982, 1983)
# NCI = 1-2Gnormed
########################################################################################
########################################################################################
Gnormed <- function(matrix){
  N <- dim(matrix)[1]; I <- dim(matrix)[2];
  NC <- apply(matrix,1,sum);
  uniqueNC <- sort(unique(NC));
  pi <- apply(matrix,2,mean);
  matrix.ord <- matrix[,order(pi,decreasing=TRUE)];
  per.row <- function(vect){
    NC <- sum(vect);
    if (NC %% I != 0) {
      sum <- 0;
      for (i in 1:(I-1)){sum <- sum + sum(diff(vect,lag=i)==1)}
      sum/(NC*(I-NC))} else {0} # all-0s or all-1s vector -> Gnormed=0
  }
  res <- apply(matrix.ord,1,per.row);
  res <- list(PFscores=as.vector(round(res,4)),PFStatistic="Gnormed");
  class(res) <- "PerFit";
  res
}

########################################################################################
########################################################################################
# NCI (Tatsuoka & Tatsuoaka, 1982, 1983):
# This statistic is perfectly linearly related to Gnormed (van der Flier, 1977; Meijer, 1994)
# NCI = 1-2Gnormed
########################################################################################
########################################################################################
NCI <- function(matrix){
  N <- dim(matrix)[1]; I <- dim(matrix)[2];
  NC <- apply(matrix,1,sum);
  uniqueNC <- sort(unique(NC));
  pi <- apply(matrix,2,mean);
  matrix.ord <- matrix[,order(pi,decreasing=TRUE)];
  per.row <- function(vect){
    NC <- sum(vect);
    if (NC %% I != 0) {
      sum <- 0;
      for (i in 1:(I-1)){sum <- sum + sum(diff(vect,lag=i)==1)}
      1 - 2*sum/(NC*(I-NC))} else {1} # all-0s or all-1s vector -> NCI=1
  }
  res <- apply(matrix.ord,1,per.row);
  res <- list(PFscores=as.vector(round(res,4)),PFStatistic="NCI");
  class(res) <- "PerFit";
  res
}

########################################################################################
########################################################################################
# Ht (Sijtsma, 1986)
########################################################################################
########################################################################################
Ht <- function(matrix){
  N <- dim(matrix)[1]; I <- dim(matrix)[2];
  NC <- apply(matrix,1,sum);
  uniqueNC <- sort(unique(NC));
  if (min(uniqueNC)==0 || max(uniqueNC)==I){
    stop('Error: Remove all-0s and/or all-1s response vectors from data.')};
  res <- c();
  singlePs <- apply(matrix,1,mean);
  tot.score <- apply(matrix,2,sum);
  num <- apply(matrix,1,function(vect){cov(vect,tot.score - vect)*(I-1)/I});
  df <- data.frame(1:N,singlePs,num);
  df.ord <- df[order(df[,2]),];
  singlePs.ord <- df.ord[[2]];
  pos <- which(diff(singlePs.ord,lag=1)>0);  
  less <- sapply(c(pos,N),function(x){sum(singlePs.ord[1:(x-1)])}); if (pos[1]==1){less[1] <- 0};
  less <- rep(less,c(pos[1],diff(c(pos,N),lag=1))) * (1-singlePs.ord);
  more <- sapply(pos,function(x){sum(1-singlePs.ord[(x+1):N])}); more <- c(more,0);
  more <- rep(more,c(pos[1],diff(c(pos,N),lag=1))) * singlePs.ord;
  den <- less + more;
  df.ord <- data.frame(df.ord,den);
  df <- df.ord[order(df.ord[,1]),];
  res <- df$num / df$den;
  res <- list(PFscores=round(res,4),PFStatistic="Ht");
  class(res) <- "PerFit";
  res
}

########################################################################################
########################################################################################
# A, D, E (Kane & Brennan, 1980)
########################################################################################
########################################################################################
A.KB <- function(matrix){
  N <- dim(matrix)[1]; I <- dim(matrix)[2];
  NC <- apply(matrix,1,sum);
  uniqueNC <- sort(unique(NC));
  if (min(uniqueNC)==0 || max(uniqueNC)==I){
    stop('Error: Remove all-0s and/or all-1s response vectors from data before proceeding.')};
  pi <- apply(matrix,2,mean);
  a <- apply(matrix,1,function(vect){vect %*% pi});
  res <- list(PFscores=as.vector(a),PFStatistic="A.KB");
  class(res) <- "PerFit";
  res
}

D.KB <- function(matrix){
  N <- dim(matrix)[1]; I <- dim(matrix)[2];
  NC <- apply(matrix,1,sum);
  uniqueNC <- sort(unique(NC));
  if (min(uniqueNC)==0 || max(uniqueNC)==I){
    stop('Error: Remove all-0s and/or all-1s response vectors from data before proceeding.')};
  pi <- apply(matrix,2,mean);
  pi.ord <- sort(pi,decreasing=TRUE);
  a <- apply(matrix,1,function(vect){vect %*% pi});
  a.max <- apply(matrix,1,function(vect){sum(pi.ord[1:sum(vect)])});
  res <- list(PFscores=as.vector(a.max-a),PFStatistic="D.KB");
  class(res) <- "PerFit";
  res
}

E.KB <- function(matrix){
  N <- dim(matrix)[1]; I <- dim(matrix)[2];
  NC <- apply(matrix,1,sum);
  uniqueNC <- sort(unique(NC));
  if (min(uniqueNC)==0 || max(uniqueNC)==I){
    stop('Error: Remove all-0s and/or all-1s response vectors from data before proceeding.')};
  pi <- apply(matrix,2,mean);
  pi.ord <- sort(pi,decreasing=TRUE);
  a <- apply(matrix,1,function(vect){vect %*% pi});
  a.max <- apply(matrix,1,function(vect){sum(pi.ord[1:sum(vect)])});
  res <- list(PFscores=as.vector(a/a.max),PFStatistic="E.KB");
  class(res) <- "PerFit";
  res
}

########################################################################################
########################################################################################
# Polytomous items: Number of Guttman errors
# Gpoly reduces to G for 0-1 data (Ncat=2) (currently Gpoly needs G to be loaded)
# Gnormed.poly reduces to Gnormed for 0-1 data (Ncat=2) (currently Gnormed.poly needs G to be loaded)
########################################################################################
########################################################################################

Gpoly <- function(matrix,Ncat) {
  N <- dim(matrix)[1]; I <- dim(matrix)[2]; M <- Ncat-1;
  probs.ISD <- matrix(NA,nrow=I,ncol=M);
  for (m in 1:M) {probs.ISD[,m] <- apply(matrix,2,function(vec) {sum(vec >= m)/N})};
  f.scoresISD <- function (x) {if (x == 0) {rep(0,M)} else {c(rep(1,x),rep(0,M-x))}};
  matrix.ISD <- matrix(unlist(lapply(t(matrix),f.scoresISD)),byrow=TRUE,nrow=N);
  probs.ISD.vect <- as.vector(t(probs.ISD));
  matrix.ISD.ord <- matrix.ISD[,order(probs.ISD.vect,decreasing=TRUE)];
  res <- G(matrix.ISD.ord)$PFscores;
  res <- list(PFscores=round(res,4),PFStatistic="Gpoly");
  class(res) <- "PerFit";
  res
}

Gnormed.poly <- function(matrix,Ncat) {
  N <- dim(matrix)[1]; I <- dim(matrix)[2]; M <- Ncat-1;
  NC <- apply(matrix,1,sum);
  probs.ISD <- matrix(NA,nrow=I,ncol=M);
  for (m in 1:M) {probs.ISD[,m] <- apply(matrix,2,function(vec) {sum(vec >= m)/N})};
  probs.ISD.vect <- as.vector(t(probs.ISD));
  ranks.ISD <- matrix(rank(I*M - probs.ISD.vect,ties.method="first"),nrow=I,byrow=TRUE);
  if (Ncat>2) {cumranks.ISD <- cbind(rep(0,I),t(apply(ranks.ISD,1,cumsum)))};
  if (Ncat==2) {cumranks.ISD <- cbind(rep(0,I),ranks.ISD)};
  V <- matrix(rep(cumranks.ISD[1,],Ncat),ncol=Ncat,byrow=FALSE);
  add.term <- matrix(rep(cumranks.ISD[2,],nrow(V)),ncol=Ncat,byrow=TRUE);
  V <- V + add.term;
  T <- sapply(2:sum(dim(V)),function(x) {max(V[col(V)+row(V) == x])});
  for (i in 3:I) {
    V <- matrix(rep(T,Ncat),ncol=Ncat,byrow=FALSE);
    add.term <- matrix(rep(cumranks.ISD[i,],nrow(V)),ncol=Ncat,byrow=TRUE);
    V <- V + add.term;
    T <- sapply(2:sum(dim(V)),function(x) {max(V[col(V)+row(V) == x])});
  }
  maxG <- T - sapply(0:(I*M),function(x) {.5*x*(x+1)});
  maxG[c(1,length(maxG))] <- 1; # so that vectors (0,0,...,0) and (M,M,...,M), which have 0 Guttman errors, are divided by 1 instead of 0
  res <- Gpoly(matrix,Ncat)$PFscores / maxG[NC+1];
  res <- list(PFscores=round(res,4),PFStatistic="Gnormed.poly");
  class(res) <- "PerFit";
  res
}

########################################################################################
########################################################################################
# Polytomous items: U3p
# U3poly reduces to U3 for 0-1 data (Ncat=2)
########################################################################################
########################################################################################
U3poly.W <- function(matrix,Ncat) {
  N <- dim(matrix)[1]; I <- dim(matrix)[2]; M <- Ncat-1;
  probs.ISD <- matrix(NA,nrow=I,ncol=M);
  for (m in 1:M) {probs.ISD[,m] <- apply(matrix,2,function(vec) {sum(vec >= m)/N})};
  f.scoresISD <- function (x) {if (x == 0) {rep(0,M)} else {c(rep(1,x),rep(0,M-x))}};
  matrix.ISD <- matrix(unlist(lapply(t(matrix),f.scoresISD)),byrow=TRUE,nrow=N);
  probs.ISD.vect <- as.vector(t(probs.ISD));
  # If there are answer options not chosen by any respondent then some entries in 'probs.ISD.vect' are 0 and others are 1.
  # Below all corresponding logs are set from Inf to 0 (reason: They carry no information regarding aberrant response behavior):
  logits.ISD <- log(probs.ISD.vect/(1-probs.ISD.vect));
  logits.ISD <- sapply(logits.ISD,function(x){if (is.infinite(x)) 0 else x}); # a vector
  as.vector(matrix.ISD %*% logits.ISD)
}

U3poly <- function(matrix,Ncat) {
  N <- dim(matrix)[1]; I <- dim(matrix)[2]; M <- Ncat-1;
  NC <- apply(matrix,1,sum);
  probs.ISD <- matrix(NA,nrow=I,ncol=M);
  for (m in 1:M) {probs.ISD[,m] <- apply(matrix,2,function(vec) {sum(vec >= m)/N})};
  probs.ISD.vect <- as.vector(t(probs.ISD));
  # If there are answer options not chosen by any respondent then some entries in 'probs.ISD.vect' are 0 and others are 1.
  # Below all corresponding logs are set from Inf to 0 (reason: They carry no information regarding aberrant response behavior):
  logits.ISD <- log(probs.ISD.vect/(1-probs.ISD.vect));
  logits.ISD <- matrix(sapply(logits.ISD,function(x){if (is.infinite(x)) 0 else x}),nrow=I,byrow=TRUE) # a matrix
  if (Ncat>2)  {cumlogits.ISD <- cbind(rep(0,I),t(apply(logits.ISD,1,cumsum)))};
  if (Ncat==2) {cumlogits.ISD <- cbind(rep(0,I),logits.ISD)};
  V <- matrix(rep(cumlogits.ISD[1,],Ncat),ncol=Ncat,byrow=FALSE);
  add.term <- matrix(rep(cumlogits.ISD[2,],nrow(V)),ncol=Ncat,byrow=TRUE);
  V <- V + add.term;
  T <- sapply(2:sum(dim(V)),function(x) {min(V[col(V)+row(V) == x])});
  for (i in 3:I) {
    V <- matrix(rep(T,Ncat),ncol=Ncat,byrow=FALSE);
    add.term <- matrix(rep(cumlogits.ISD[i,],nrow(V)),ncol=Ncat,byrow=TRUE);
    V <- V + add.term;
    T <- sapply(2:sum(dim(V)),function(x) {min(V[col(V)+row(V) == x])});
  }
  #maxW <- sapply(0:(I*M),function(x) {sum(sort(log(probs.ISD.vect/(1-probs.ISD.vect)),decreasing=TRUE)[0:x])});
  maxW <- sapply(0:(I*M),function(x) {sum(sort(as.vector(logits.ISD),decreasing=TRUE)[0:x])});
  minW <- T;
  den <- maxW - minW;
  den[c(1,length(den))] <- 1; # so that vectors (0,0,...,0) and (M,M,...,M), which have 0 Guttman errors, are divided by 1 instead of 0
  res <- (maxW[NC+1] - U3poly.W(matrix,Ncat)) / den[NC+1];
  res <- list(PFscores=round(res,4),PFStatistic="U3poly");
  class(res) <- "PerFit";
  res
}

########################################################################################
########################################################################################
# lz, lzstar
########################################################################################
########################################################################################

lz <- function (matrix,ip=NA,model="2PL",ability=NA,method="ML",mu=0,sigma=1) {
  # Estimate item parameters if not provided:
  if (is.na(ip)[1]) {ip <- est(matrix,model,engine="ltm",rasch=TRUE,nqp=20)$est}
  # Estimate ability parameters if not provided (using 'irtoys'):
  if (is.na(ability)[1]) {
    ability <- switch(method,
                      ML=mlebme(matrix,ip,mu,sigma,method="ML")[,1],
                      BM=mlebme(matrix,ip,mu,sigma,method="BM")[,1],
                      WL=wle(matrix,ip)[,1])
  }
  #  
  N <- dim(matrix)[1]; I <- dim(matrix)[2];
  A <- ip[,1]; B <- ip[,2]; C <- ip[,3];
  P <- t(sapply(1:N,function (x) {C+(1-C) / (1+exp(-A*(ability[x] - B)))})); Q <- 1-P;
  l0 <- apply(matrix*log(P) + (1-matrix)*log(Q),1,sum);
  El0 <- apply(P*log(P) + Q*log(Q),1,sum);
  Vl0 <- apply(P*Q*(log(P/Q))^2,1,sum);
  res <- list(PFscores=as.vector(round((l0 - El0) / sqrt(Vl0),4)),PFStatistic="lz");
  class(res) <- "PerFit";
  res
}

lzstar <- function (matrix,ip=NA,model="2PL",method="ML",mu=0,sigma=1) {
  # Estimate item parameters if not provided:
  if (is.na(ip)[1]) {ip <- est(matrix,model,engine="ltm",rasch=TRUE)$est}
  # Estimate ability parameters using 'irtoys':
  ability <- switch(method,
                    ML=mlebme(matrix,ip,mu,sigma,method="ML")[,1],
                    BM=mlebme(matrix,ip,mu,sigma,method="BM")[,1],
                    WL=wle(matrix,ip)[,1]);  
  #  
  N <- dim(matrix)[1]; I <- dim(matrix)[2];
  A <- ip[,1]; B <- ip[,2]; C <- ip[,3];
  P <- t(sapply(1:N,function (x) {C+(1-C) / (1+exp(-A*(ability[x] - B)))})); Q <- 1-P;
  d1P <- t(sapply(1:N,function (x) {(1-C)*A*exp(A*(ability[x] - B)) / (1+exp(A*(ability[x] - B)))^2}));
  d2P <- t(sapply(1:N,function (x) {(1-C)*(A^2)*exp(A*(ability[x] - B))*(1-exp(A*(ability[x] - B))) / (1+exp(A*(ability[x] - B)))^3}));
  ri <- d1P/(P*Q);
  r0 <- switch(method,
               ML=0,
               BM=(mu-ability) / (sigma^2),
               WL=apply((d1P*d2P)/(P*Q),1,sum) / (2*apply((d1P^2)/(P*Q),1,sum)));
  wi <- log(P/Q);
  Wn <- apply((matrix - P)*wi,1,sum);
  sigma2n <- apply((wi^2)*P*Q,1,sum) / I;
  cn <- apply(d1P*wi,1,sum) / apply(d1P*ri,1,sum);
  wi.tilde <- wi - diag(cn) %*% ri;
  tau2n <- apply((wi.tilde^2)*P*Q,1,sum) / I;
  EWn <- -cn * r0;
  VWn <- I * tau2n;
  res <- list(PFscores=as.vector(round((Wn - EWn) / sqrt(VWn),4)),PFStatistic="lzstar");
  class(res) <- "PerFit";
  res
}

########################################################################################
########################################################################################
# lzp
########################################################################################
########################################################################################
lzpoly <- function (matrix,Ncat,ip=NA,model="GRM",ability=NA,method="EAP") {
  if (sum(is.na(ip)[1] + is.na(ability)[1]) == 1){
    stop('Provide either no estimated parameters (i.e., "ip=NA" and "ability=NA") or all 
         estimated parameters (both "ip" and "ability").')};
  #
  N <- dim(matrix)[1]; I <- dim(matrix)[2]; M <- Ncat-1;
  matrix2 <- data.frame(apply(matrix,2,as.factor)); # eliminates item levels with no answers
  # Estimate item parameters if not provided:
  if (is.na(ip)[1]) {
    ip <- switch(model,
                 PCM=gpcm(matrix2,constraint="rasch",IRT.param=TRUE),
                 GPCM=gpcm(matrix2,constraint="gpcm",IRT.param=TRUE),
                 GRM=grm(matrix2,constrained=FALSE,IRT.param=TRUE));
    ip.coef <- coef(ip);
  } else {ip.coef <- ip};  
  if (is.list(ip.coef)) { # NOT all answer categories of all items were used
    abs.freqs <- apply(matrix,2,table);
    abs.freqs <- lapply(abs.freqs,function (vect) as.numeric(names(vect)));
    tmp <- matrix(NA,nrow=I,ncol=Ncat);
    for (i in 1:I) {
      tmp[i,abs.freqs[[i]][-length(abs.freqs[[i]])]+1] <- ip.coef[[i]][-length(abs.freqs[[i]])];
      tmp[i,Ncat] <- ip.coef[[i]][length(ip.coef[[i]])]
    }
    ip.coef <- tmp;
  }
  # Estimate ability parameters if not provided (using 'ltm'):
  if (is.na(ability)[1]) {
    ability <- factor.scores(ip,resp.patterns=matrix2,method=method);
    ability <- ability$score.dat[,ncol(ability$score.dat)-1]; 
  }
  #  
  if (model == "GRM") {
    # P.ISRF is N x I*M:
    P.ISRF <- t(sapply(ability,function(x){as.vector(t(1/(1+exp(-ip.coef[,ncol(ip.coef)]*(x-ip.coef[,-ncol(ip.coef)])))))}));
    # Fix for datasets with non-chosen answer options (the NA entries):
    #   1s for NAs in the first item steps
    #   0s for NAs in the last item steps
    #   entry (x+1) for NAs in entry x
    if (sum(is.na(ip.coef)) > 0) {
      first.cols <- (which(is.na(ip.coef[,1]))-1)*M+1; P.ISRF[,first.cols][is.na(P.ISRF[,first.cols])] <- 1;
      last.cols <- which(is.na(ip.coef[,M]))*M; P.ISRF[,last.cols][is.na(P.ISRF[,last.cols])] <- 0;
      middle.cols <- sort(which(is.na(t(cbind(rep(0,I),ip.coef[,-c(1,M,Ncat)],rep(0,I))))),decreasing=TRUE);
      for (i in 1:length(middle.cols)){P.ISRF[,middle.cols] <- P.ISRF[,middle.cols+1]}
    }
    P.CRF <- matrix(,nrow=N,ncol=I*Ncat);
    for (i in 1:I) {
      P.ISRF.item <- cbind(rep(1,N),P.ISRF[,((i-1)*M+1):(i*M)],rep(0,N));
      P.CRF[,((i-1)*Ncat+1):(i*Ncat)] <- t(apply(-P.ISRF.item,1,diff))
    }
  }
  if (model == "PCM" | model == "GPCM") {
    lin.it <- t(sapply(ability,function(x){as.vector(t(ip.coef[,ncol(ip.coef)]*(x-ip.coef[,-ncol(ip.coef)])))}));
    lin.it[,is.na(as.vector(t(ip.coef[,-ncol(ip.coef)])))] <- 0; # NAs -> 0 to eliminate these terms from the sums
    P.CRF <- matrix(,nrow=N,ncol=I*Ncat);
    for (i in 1:I) {
      num <- t(exp(apply(cbind(rep(0,N),lin.it[,((i-1)*M+1):(i*M)]),1,cumsum)));
      P.CRF[,((i-1)*Ncat+1):(i*Ncat)] <- num / apply(num,1,sum)
    }  
  }
  f.scores <- function (x) {vec<- rep(0,Ncat); vec[x+1] <- 1; vec};
  matrix.01 <- matrix(unlist(lapply(t(matrix),f.scores)),byrow=TRUE,nrow=N);
  # If there are answer options not chosen by any respondent then some entries in 'P.CRF' might be 0.
  # Below all corresponding logs are set from Inf to 0 (reason: They carry no information regarding aberrant response behavior):
  log.P.CRF <- log(P.CRF);
  log.P.CRF <- matrix(sapply(log.P.CRF,function(x){if (is.infinite(x)) 0 else x}),nrow=N,byrow=FALSE) # a matrix
  #
  l0p <- apply(matrix.01 * log.P.CRF,1,sum);
  El0p <- apply(P.CRF * log.P.CRF,1,sum);
  V.row <- function(vect) {
    tot <- 0;
    for (i in 1:I) {
      vect.part <- vect[((i-1)*Ncat+1):(Ncat*i)];
      log.vect.part <- log(vect.part);
      # Convert Inf's to 0 (reason: They carry no information regarding aberrant response behavior):
      log.vect.part <- sapply(log.vect.part,function(x){if (is.infinite(x)) 0 else x})
      part1 <- (as.matrix(vect.part) %*% vect.part);
      part2 <- matrix(rep(log.vect.part,Ncat),nrow=Ncat);
      part3 <- matrix(rep(log.vect.part,Ncat),nrow=Ncat) - matrix(rep(log.vect.part,Ncat),nrow=Ncat,byrow=TRUE);
      #part3 <- log(as.matrix(vect.part) %*% (1/vect.part));
      tot <- tot + sum(part1 * part2 * part3);
    }
    tot
  }
  Vl0p <- apply(P.CRF,1,V.row);  
  res <- list(PFscores=round((l0p - El0p) / sqrt(Vl0p),4),PFStatistic="lzpoly");
  class(res) <- "PerFit";
  res
  }















