PRF <- function(matrix,respID,h=.09,N.FPts=101) {
  matrix <- as.matrix(matrix);
  Diffs <- 1 - apply(matrix,2,mean);
  focal.pts <- seq(0,1,length.out=N.FPts);
  GaussKernel <- function(x) {dnorm(x,mean=0,sd=1)};
  KernArg <- expand.grid(focal.pts,Diffs);
  KernArg <- matrix(KernArg[,1] - KernArg[,2],nrow=length(focal.pts),byrow=F) / h;
  weights <- GaussKernel(KernArg) / apply(GaussKernel(KernArg),1,sum);
  list(PRFdiffs=focal.pts,PRFest=as.vector(weights %*% matrix[respID,]))
}

PRF.VarBands <- function (matrix,respID,h=.09,N.FPts=101,alpha=.05) {
  matrix <- as.matrix(matrix);
  focal.pts <- seq(0,1,length.out=N.FPts);
  I <- dim(matrix)[2];
  PRFscores <- as.vector(PRF(matrix,respID,h,N.FPts)$PRFest);
  # Jackknife estimate of the SE:
  PRF.SEmat <- matrix(rep(NA,length(focal.pts) * I),nrow=length(focal.pts));
  for (it in 1:I){PRF.SEmat[,it] <- PRF(matrix[,-it],respID,h,N.FPts)$PRFest};
  PRF.SE <- sqrt( ((I-1)/I) * apply((PRF.SEmat - apply(PRF.SEmat,1,mean))^2,1,sum) );
  crit.val <- qnorm(1-alpha,mean=0,sd=1);
  list(PRF.VarBandsLow=PRFscores-crit.val*PRF.SE,PRF.VarBandsHigh=PRFscores+crit.val*PRF.SE)
}

PRFplot <- function (matrix,respID,h=.09,N.FPts=101,VarBands=FALSE,VarBands.area=FALSE,alpha=.05,
                     Xlabel=NA,Xcex=1.5,Ylabel=NA,Ycex=1.5,title=NA,Tcex=1.5) {
  matrix <- as.matrix(matrix);
  res1 <- PRF(matrix,respID,h,N.FPts);
  res2 <- PRF.VarBands(matrix,respID,h,N.FPts,alpha);
  par(mar=c(4,4,2,1)+.1,las=1);
  plot(res1$PRFdiffs,res1$PRFest,lwd=2,type="l",axes=F,ann=F,frame.plot=T,xlim=c(0,1),ylim=c(0,1));
  tmpx <- if (is.na(Xlabel)) {"Item difficulty"} else {Xlabel};
  axis(1,at=seq(0,1,by=.2)); mtext(side=1,text=tmpx,line=2.5,col="black",cex=Xcex,font=1);
  tmpy <- if (is.na(Ylabel)) {"Probability correct answer"} else {Ylabel};
  axis(2,at=seq(0,1,by=.2)); mtext(side=2,text=tmpy,line=2.8,col="black",cex=Ycex,font=1,las=3);
  if (VarBands == TRUE & VarBands.area == TRUE) {
    polygon(c(res1$PRFdiffs,rev(res1$PRFdiffs)), c(res2$PRF.VarBandsHigh, rev(res2$PRF.VarBandsLow)),col = "lightpink1",border=NA);
    points(res1$PRFdiffs,res1$PRFest,lwd=2,type="l",ann=F);
    points(res1$PRFdiffs,res2$PRF.VarBandsLow,type="l",ann=F,lty=2,lwd=1.5);
    points(res1$PRFdiffs,res2$PRF.VarBandsHigh,type="l",ann=F,lty=2,lwd=1.5);
  }
  if (VarBands == TRUE & VarBands.area == FALSE) {
    points(res1$PRFdiffs,res2$PRF.VarBandsLow,type="l",ann=F,lty=2,lwd=1.5);
    points(res1$PRFdiffs,res2$PRF.VarBandsHigh,type="l",ann=F,lty=2,lwd=1.5);
  }
  tmp <- if (is.na(title)) {paste("PRF (respID # ",respID,")",sep="")} else {title};
  mtext(side=3,text=tmp,line=.5,col="black",cex=Tcex,font=2);
}




