#!/usr/bin/R
#
#
#   $Id: graph.R,v 1.2 2011/02/09 09:34:59 artigue Exp $
#
#
#    USAGE: R-dev --vanilla --slave --args <$i> < $SRC/graph.R 
#
#    OUTPUT: graph_<i>.png, sdps.dat, <i>.pvalue
#
#    NEED: conservation.dat  <zscore.dat>
#


args <- commandArgs(TRUE)


i <- grep("--args", commandArgs()) 
n <- as.numeric(commandArgs()[i <- i + 1]) 

conservation.dat<-"conservation.dat"
zscore.dat<-paste(n,".zscore",collapse = NULL,sep="")
graph<-paste("graph_",n,".png", collapse = NULL,sep="")


# Read data
conservation <-read.table(conservation.dat)
zscore <-read.table(zscore.dat)


# Deduce number of amino-acids in the pocket
param1<-length(conservation[,1])

pvalue <--pnorm(zscore[,5],lower.tail=T,log.p=T)
to <-sort(conservation[,3]*0.01,index.return = TRUE)

filename<-paste(n,".pvalue",sep="")
cat("",file=filename,append=FALSE)

filesdp<-"sdps.dat"

for(kk in 1:param1){
tmp <- pvalue[kk]
cat(kk, tmp, "\n",file=filename,append=TRUE)
if(tmp < 1e-2){cat(kk,"\t",n,"\n",file=filesdp, append=TRUE) }
}


# Graphics
options(bitmapType='cairo')
bitmap(graph,type="png16m", height = 1500, width = 86*param1, res = 72, units = "px", pointsize = 32)
#png(file=graph,bg="white", height = 1500, units = "px",width=86*param1,pointsize = 32,res=NA)

plot(-pnorm(zscore[,5],lower.tail=T,log.p=T), type="h", xlab="AA Position in the Pocket", ylab="P-value",main="ASMC",col = "blue", lwd=40, axes=TRUE,log="y", ylim=c(1e-10,1.1))

abline(h=1e-4,col="red",lty=3,lwd=4)
par(new=TRUE)

plot(conservation[,3], type="h",col = "green", lwd=15,ylim=c(0,100),axes=FALSE,ylab="",xlab="")

plot_colors <-c("green", "blue")

legend("topleft", c("Conserv.","P-value"), cex=0.8, bty="n", lty=1, lwd=5, col=plot_colors)


dev.off()

q()
