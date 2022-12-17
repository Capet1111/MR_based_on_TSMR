############defination of forrestplot function############
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
  #import the data
  rt <- read.table(coxFile, header=T, sep=",", check.names=F)
  exposure <- rt[,1]
  method <- rt[,2]
  or <- sprintf("%.3f",rt$"OR")
  orLow  <- sprintf("%.3f",rt$"OR.95L")
  orHigh <- sprintf("%.3f",rt$"OR.95H")
  odds.ratio <- paste0(or,"(",orLow,"-",orHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #export the image
  pdf(file=forestFile, width=6.6, height=4.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #plot the clinical information
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  
  text.cex=0.8
  text(0,n:1,exposure,adj=0,cex=text.cex)
  text(0,n:1,method,adj=0,cex=text.cex)
  
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex)
  text(1.5-0.5*0.2,n+1,'p-value',cex=text.cex,font=2,adj=1)
  
  text(3.1,n:1,odds.ratio,adj=1,cex=text.cex)
  text(3.1,n+1,'Odds ratio',cex=text.cex,font=2,adj=1)
  
  #plot forrest
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Odds ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
  dev.off()
}

