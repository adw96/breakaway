betta_pic <- function(y, se, x=1:length(y), ylimu=NA, myy=NA, mymain=NA, mycol=rep("black", length(y)), labs=NA, mypch=rep(16, length(y)), myxlim=c(0.8*min(x, na.rm=T), 1.2*max(x, na.rm=T))) {
  n <- length(y)
  ylimu <- ifelse(is.na(ylimu), max(y+2*se), ylimu)
  par(xpd=NA)
  plot(0, 0, type="n", xlim=myxlim, ylim=c(0, ylimu), xlab="", bty="n", ylab=myy, main=mymain)
  for (i in 1:n) {
    if(!is.na(y[i]) & !is.na(x[i])) {
      points(x[i], y[i], pch=mypch[i], col=mycol[i])
      lines(c(x[i], x[i]), c(max(0, y[i]-1.96*se[i]), y[i]+1.96*se[i]), col=mycol[i])
    }
  }
  if(!is.na(labs)) axis(1, at=1:length(y), labels=labs, las=2, cex=0.8)
}
