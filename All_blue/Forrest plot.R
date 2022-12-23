#set plot information
File="mr_OR.csv"
#set arguments
forestFile = "test.pdf"
Width = 10
Height = 12

arrowCol = "#000000"
boxCol = c("#d90429","#06d6a0")
xlim_right = 3.5


##import the data---------------------------------------------------------------
rt <- na.omit(read.table(File, header=T, sep=",", check.names=F, fill = T))
dat = rt
dat$odds.ratio = sprintf("%s (%s-%s)",
                         sprintf("%.3f",rt$"OR"),
                         sprintf("%.3f",rt$"OR.95L"),
                         sprintf("%.3f",rt$"OR.95H"))
dat$pval = ifelse(rt$pval < 0.001, "<0.001", sprintf("%.3f", rt$pval))
dat$arrow_OR.95H =as.numeric(ifelse(dat$OR.95H >= xlim_right, xlim_right, sprintf("%.3f", dat$OR.95H)))

boxCol = rep(boxCol,times = length(levels(as.factor(dat$exposure))),each = length(levels(as.factor(dat$method))))


##Create a PDF file-------------------------------------------------------------
pdf(file=forestFile, width=Width, height = Height)

#set y axis limitation based on row number
n <- nrow(dat); nRow <- n+1 ; ylim <- c(1,nRow)

#seperate the picture
layout(matrix(c(1,2), nc = 2),width=c(4,2) )


##plot the clinical information-------------------------------------------------
xlim = c(0,12); text.cex=1; par(mar=c(4,1,1,0))#set marginal
#plot blank picture
plot(1, xlim=xlim, ylim=ylim, type="n", axes=F, xlab="", ylab="")
#text clinical information
#[1]
x1 = 0
text(x=x1, y=n:1,dat$exposure, adj=0, cex=text.cex)
text(x=x1, y=n+1,'Exposure', cex=text.cex, font=2,adj=0)

#[2]
x1 = x1 + max(nchar(dat$exposure))/8 +1.5
text(x=x1, y=n:1,dat$method, adj=0, cex=text.cex)
text(x=x1, y=n+1,'Method', cex=text.cex, font=2,adj=0)

#[3]
x1 = x1 + max(nchar(dat$method))/8 +1.5
text(x=x1, y=n:1, dat$pval,adj=0,cex=text.cex)
text(x=x1, y=n+1, 'P value',cex=text.cex,font=2,adj=0)

#[4]
x1 = x1 + max(nchar(dat$pval))/8 +1.5
text(x=x1, y=n:1,dat$"odds.ratio",adj=0,cex=text.cex)
text(x=x1, y=n+1,'Odds ratio',cex=text.cex,font=2,adj=0)
#[5]


##plot forest-------------------------------------------------------------------
par(mar=c(4,0.5,1,2), mgp=c(2,0.5,0))
xlim_forest = c(0,xlim_right)


#plot blank picture
plot(1,xlim=xlim_forest,ylim=ylim, type="n",axes=F, ylab="", xaxs="i", xlab="Odds ratio")

#draw arrows
arrows(dat$OR.95L, n:1, 
       dat$arrow_OR.95H, n:1,   #coordinates of arrow
       angle=90, code=3, length=0.03, col=arrowCol, lwd=2)#shape of arrow
#draw points
points(x=dat$OR, y=n:1, cex=1.5,
       pch = 15, #shape of point
       col = boxCol #color of point
)

#draw none_effect line 
abline(v=1,col="black", lty=2, lwd=2); axis(1)
#-------------------------------------------------------------------------------
dev.off()





