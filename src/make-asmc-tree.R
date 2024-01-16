#! /home/eelisee/miniconda3/envs/asmc/bin/R
#
#    $Id: make-asmc-tree.R,v 1.7 2011/02/09 14:57:36 artigue Exp $	
#
#
#    USAGE: R-dev --vanilla --slave --args $REP < $SRC/make-asmc-tree.R 
#
#    OUTPUT: graph.png
#
#    ARGUMENTS: Number of clusters 
#
#    NEED: ImageMagick
#
#
#
 
library(pixmap)
library(ape)

args <- commandArgs(TRUE)
i <- grep("--args", commandArgs())
file <- as.character(commandArgs()[i <- i + 1])
data.tree<-read.nexus(file)


leafs<-data.tree$tip.label
nbclus<-length(data.tree$tip.label)

bitmap("asmc-tree.png",type="png16m",height = 3600, width = 2400, res = 72, units = "px",)

# space occupancy for layout:
left<-rep(1, nbclus)
middle<-rep(2, nbclus)
right<- 3:(nbclus+2)

m<-matrix(c(left,middle,right),nbclus, 3, byrow = FALSE)
layout(m,widths=c(2,1,2))

## Big summary logo on left
par(mar = c(0.1, 0.1, 0.1, 0.1))

big<- read.pnm("allseq.ppm")

plot(big)

par(mar = c(0.1, 0.1, 1.2, 0.1))
## Dendrogram in middle

#if (nbclus > 1)
 plot(data.tree, edge.width =2)


par(mar = c(0.01, 0.01, 0.01, 0.01))

## Sublogos on the right side
for(i in nbclus:1)
{
clus<-leafs[i]
fileppm<-paste(clus,".ppm",sep="")
small<- read.pnm(fileppm)
plot(small)
par(mar = c(0.1, 0.1, 0.1, 0.1))
}

dev.off()

q()
