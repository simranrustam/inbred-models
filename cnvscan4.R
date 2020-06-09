library(cobs)
library(quantreg)
library(parallel)
library(gplots)
library(limma)
library(grid)
library(gridExtra)

interpolate<-function(dat,curve){

interpolate_points<-function(row,dat,curve){
    MY_X=dat[row,1]
    MY_Y=dat[row,2]
    VAL1=tail(which(curve[,1]<MY_X),1)
    VAL2=VAL1+1
    X <- curve[c(VAL1,VAL2),1]
    Y <- curve[c(VAL1,VAL2),2]
    INTERP_Y=approx(X,Y,xout=MY_X)$y
    INTERP_Y
}
res<-unlist(mclapply(seq(nrow(dat)),interpolate_points,dat=dat,curve=curve,mc.cores=8))
res
}

################################################
# prepare the analysis
################################################

# figure out what the project is
IDX=grep("RP",rev(unlist(strsplit(getwd(),"/")))[1:2])
PROJ=rev(unlist(strsplit(getwd(),"/")))[IDX]

TG=paste(PROJ,".genes_fmt.tsv",sep="")
TE6=paste(PROJ,".1e6_fmt.tsv",sep="")
TE5=paste(PROJ,".1e5_fmt.tsv",sep="")
TE4=paste(PROJ,".1e4_fmt.tsv",sep="")
TE3=paste(PROJ,".1e3_fmt.tsv",sep="")

CORG=paste("plots/",PROJ,".genes_cor.pdf",sep="")
CORE6=paste("plots/",PROJ,".1e6_cor.pdf",sep="")
CORE5=paste("plots/",PROJ,".1e5_cor.pdf",sep="")
CORE4=paste("plots/",PROJ,".1e4_cor.pdf",sep="")
CORE3=paste("plots/",PROJ,".1e3_cor.pdf",sep="")

REGG=paste(PROJ,".genes_regions.tsv",sep="")
REGE6=paste(PROJ,".1e6_regions.tsv",sep="")
REGE5=paste(PROJ,".1e5_regions.tsv",sep="")
REGE4=paste(PROJ,".1e4_regions.tsv",sep="")
REGE3=paste(PROJ,".1e3_regions.tsv",sep="")

CVG=paste("plots/",PROJ,".genes_cv.pdf",sep="")
CVE6=paste("plots/",PROJ,".1e6_cv.pdf",sep="")
CVE5=paste("plots/",PROJ,".1e5_cv.pdf",sep="")
CVE4=paste("plots/",PROJ,".1e4_cv.pdf",sep="")
CVE3=paste("plots/",PROJ,".1e3_cv.pdf",sep="")

HMG=paste("plots/",PROJ,".genes_hm.pdf",sep="")
HME6=paste("plots/",PROJ,".1e6_hm.pdf",sep="")
HME5=paste("plots/",PROJ,".1e5_hm.pdf",sep="")
HME4=paste("plots/",PROJ,".1e4_hm.pdf",sep="")
HME3=paste("plots/",PROJ,".1e3_hm.pdf",sep="")

HCG=paste(PROJ,".genes_hc.txt",sep="")
HCE6=paste(PROJ,".1e6_hc.txt",sep="")
HCE5=paste(PROJ,".1e5_hc.txt",sep="")
HCE4=paste(PROJ,".1e4_hc.txt",sep="")
HCE3=paste(PROJ,".1e3_hc.txt",sep="")

dir.create("plots")

################################################
# genes
################################################

x<-read.table(TG,header=T,row.names=1)
x<-x[which(rowSums(x)>=10),]
x<-sweep(x, 2, colSums(x), FUN="/")*1000000

pdf(CORG)
plotMDS(x,cex=0.8,main="MDS plot")
hist(as.numeric(cor(x)),main="Pearson correlation")
hist(as.numeric(dist(t(x))),main="Eucledian distance")
mycor<-cbind(summary(as.numeric(cor(x)))  ,  summary(as.numeric(dist(t(x)))) )
colnames(mycor)<-c("Pearson","Euclidean")
mycor<-signif(mycor,4)
grid.newpage()
grid.table(mycor)
dev.off()

mysd<-apply(x,1,sd)
mean<-apply(x,1,mean)
y<-data.frame(log10(mean),mysd/mean)
colnames(y)=c("logMean","cv")
Rbs.9 <- cobs(y$logMean,y$cv, nknots=20,constraint="none",tau=0.99)
Rbs.median <- cobs(y$logMean,y$cv,nknots=20,constraint="none",tau=0.5)
pred<-data.frame(predict(Rbs.9))

y$interpolated<-interpolate(y,pred)
y$diff=y$cv-y$interpolated
yy <- y[which(y$diff>0.05),]
yy <- yy[order(-yy$diff),]
yy <- head(yy,100)

write.table(yy,file=REGG)

pdf(CVG)
plot(y$logMean,y$cv,pch=18,cex=0.5,xlab="log10(mean)",ylab="CV")
lines(predict(Rbs.9), col = "red", lwd = 1.0)
lines(predict(Rbs.median), col = "blue", lwd = 1.0)
points(yy$logMean,yy$cv)
text(yy$logMean,yy$cv+0.02,labels=rownames(yy),cex=0.5)
dev.off()

pdf(HMG)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
zz<-x[which(rownames(x) %in% rownames(yy)),]
heatmap.2(as.matrix(zz),margin=c(8, 22),cexRow=0.3,trace="none",
  cexCol=0.8,col=my_palette,scale="row")

heatmap.2(cor(t(zz)),trace="none",scale="none",margins=c(12,12),
  cexRow=0.3, cexCol=0.3)

hr <- hclust(as.dist(1-cor(t(zz), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- rev(redgreen(75))

colfunc <- colorRampPalette(c("blue", "white", "red"))
heatmap.2(as.matrix(zz), main="Hierarchical Cluster", Rowv=as.dendrogram(hr), Colv=NA,
 dendrogram="row", scale="row", col = colfunc(15), density.info="none", trace="none",
 RowSideColors= myClusterSideBar, margin=c(8, 22) )
dev.off()

write.table(mycl,file=HCG,quote=F,sep="\t")

################################################
# 1Mbp bins
################################################
x<-read.table(TE6,header=T,row.names=1)
x<-x[which(rowSums(x)>=10),]
x<-sweep(x, 2, colSums(x), FUN="/")*1000000

pdf(CORE6)
plotMDS(x,cex=0.8,main="MDS plot")
hist(as.numeric(cor(x)),main="Pearson correlation")
hist(as.numeric(dist(t(x))),main="Eucledian distance")
mycor<-cbind(summary(as.numeric(cor(x)))  ,  summary(as.numeric(dist(t(x)))) )
colnames(mycor)<-c("Pearson","Euclidean")
mycor<-signif(mycor,4)
grid.newpage()
grid.table(mycor)
dev.off()

mysd<-apply(x,1,sd)
mean<-apply(x,1,mean)
y<-data.frame(log10(mean),mysd/mean)
colnames(y)=c("logMean","cv")
Rbs.9 <- cobs(y$logMean,y$cv, nknots=20,constraint="none",tau=0.99)
Rbs.median <- cobs(y$logMean,y$cv,nknots=20,constraint="none",tau=0.5)
pred<-data.frame(predict(Rbs.9))

y$interpolated<-interpolate(y,pred)
y$diff=y$cv-y$interpolated
yy <- y[which(y$diff>0.05),]
yy <- yy[order(-yy$diff),]
yy <- head(yy,200)

write.table(yy,file=REGE6)

pdf(CVE6)
plot(y$logMean,y$cv,pch=18,cex=0.5,xlab="log10(mean)",ylab="CV")
lines(predict(Rbs.9), col = "red", lwd = 1.0)
lines(predict(Rbs.median), col = "blue", lwd = 1.0)
points(yy$logMean,yy$cv)
#text(yy$logMean,yy$cv+0.02,labels=rownames(yy),cex=0.8)
dev.off()

#pdf(HME6)
#my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
#zz<-x[which(rownames(x) %in% rownames(yy)),]
#heatmap.2(as.matrix(zz),margin=c(8, 22),cexRow=0.65,trace="none",
#  cexCol=0.8,col=my_palette,scale="row")

#heatmap.2(cor(t(zz)),trace="none",scale="none",margins=c(12,12),
#  cexRow=0.8, cexCol=0.8)
#dev.off()

################################################
# 100 kbp bins
################################################

x<-read.table(TE5,header=T,row.names=1)
x<-x[which(rowSums(x)>=10),]
x<-sweep(x, 2, colSums(x), FUN="/")*1000000

pdf(CORE5)
plotMDS(x,cex=0.8,main="MDS plot")
hist(as.numeric(cor(x)),main="Pearson correlation")
hist(as.numeric(dist(t(x))),main="Eucledian distance")
mycor<-cbind(summary(as.numeric(cor(x)))  ,  summary(as.numeric(dist(t(x)))) )
colnames(mycor)<-c("Pearson","Euclidean")
mycor<-signif(mycor,4)
grid.newpage()
grid.table(mycor)
dev.off()

mysd<-apply(x,1,sd)
mean<-apply(x,1,mean)
y<-data.frame(log10(mean),mysd/mean)
colnames(y)=c("logMean","cv")
Rbs.9 <- cobs(y$logMean,y$cv, nknots=20,constraint="none",tau=0.99)
Rbs.median <- cobs(y$logMean,y$cv,nknots=20,constraint="none",tau=0.5)
pred<-data.frame(predict(Rbs.9))

y$interpolated<-interpolate(y,pred)
y$diff=y$cv-y$interpolated
yy <- y[which(y$diff>0.05),]
yy <- yy[order(-yy$diff),]
yy <- head(yy,200)

write.table(yy,file=REGE5)

pdf(CVE5)
plot(y$logMean,y$cv,pch=18,cex=0.5,xlab="log10(mean)",ylab="CV")
lines(predict(Rbs.9), col = "red", lwd = 1.0)
lines(predict(Rbs.median), col = "blue", lwd = 1.0)
points(yy$logMean,yy$cv)
text(yy$logMean,yy$cv+0.02,labels=rownames(yy),cex=0.8)
dev.off()

pdf(HME5)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
zz<-x[which(rownames(x) %in% rownames(yy)),]
heatmap.2(as.matrix(zz),margin=c(8, 22),cexRow=0.65,trace="none",
  cexCol=0.8,col=my_palette,scale="row")

heatmap.2(cor(t(zz)),trace="none",scale="none",margins=c(12,12),
  cexRow=0.8, cexCol=0.8)
dev.off()

################################################
# 10 kbp bins
################################################

x<-read.table(TE4,header=T,row.names=1)
x<-x[which(rowSums(x)>=10),]
x<-sweep(x, 2, colSums(x), FUN="/")*1000000

pdf(CORE4)
plotMDS(x,cex=0.8,main="MDS plot")
hist(as.numeric(cor(x)),main="Pearson correlation")
hist(as.numeric(dist(t(x))),main="Eucledian distance")
mycor<-cbind(summary(as.numeric(cor(x)))  ,  summary(as.numeric(dist(t(x)))) )
colnames(mycor)<-c("Pearson","Euclidean")
mycor<-signif(mycor,4)
grid.newpage()
grid.table(mycor)
dev.off()

mysd<-apply(x,1,sd)
mean<-apply(x,1,mean)
y<-data.frame(log10(mean),mysd/mean)
colnames(y)=c("logMean","cv")
Rbs.9 <- cobs(y$logMean,y$cv, nknots=20,constraint="none",tau=0.99)
Rbs.median <- cobs(y$logMean,y$cv,nknots=20,constraint="none",tau=0.5)
pred<-data.frame(predict(Rbs.9))

y$interpolated<-interpolate(y,pred)
y$diff=y$cv-y$interpolated
yy <- y[which(y$diff>0.05),]
yy <- yy[order(-yy$diff),]
yy <- head(yy,200)

write.table(yy,file=REGE4)

pdf(CVE4)
plot(y$logMean,y$cv,pch=18,cex=0.5,xlab="log10(mean)",ylab="CV")
lines(predict(Rbs.9), col = "red", lwd = 1.0)
lines(predict(Rbs.median), col = "blue", lwd = 1.0)
points(yy$logMean,yy$cv)
text(yy$logMean,yy$cv+0.02,labels=rownames(yy),cex=0.8)
dev.off()

pdf(HME4)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
zz<-x[which(rownames(x) %in% rownames(yy)),]
heatmap.2(as.matrix(zz),margin=c(8, 22),cexRow=0.65,trace="none",
  cexCol=0.8,col=my_palette,scale="row")

heatmap.2(cor(t(zz)),trace="none",scale="none",margins=c(12,12),
  cexRow=0.8, cexCol=0.8)
dev.off()

hr <- hclust(as.dist(1-cor(t(zz), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- rev(redgreen(75))

colfunc <- colorRampPalette(c("blue", "white", "red"))
heatmap.2(as.matrix(zz), main="Hierarchical Cluster", Rowv=as.dendrogram(hr), Colv=NA,
 dendrogram="row", scale="row", col = colfunc(15), density.info="none", trace="none",
 RowSideColors= myClusterSideBar, margin=c(8, 22) )
dev.off()

write.table(mycl,file=HCE4,quote=F,sep="\t")

################################################
# 1 kbp bins
################################################

# skip this
if(FALSE){
x<-read.table(TE3,header=T,row.names=1)
x<-x[which(rowSums(x)>=10),]
x<-sweep(x, 2, colSums(x), FUN="/")*1000000

pdf(CORE3)
plotMDS(x,cex=0.8,main="MDS plot")
hist(as.numeric(cor(x)),main="Pearson correlation")
hist(as.numeric(dist(t(x))),main="Eucledian distance")
mycor<-cbind(summary(as.numeric(cor(x)))  ,  summary(as.numeric(dist(t(x)))) )
colnames(mycor)<-c("Pearson","Euclidean")
mycor<-signif(mycor,4)
grid.newpage()
grid.table(mycor)
dev.off()

mysd<-apply(x,1,sd)
mean<-apply(x,1,mean)
y<-data.frame(log10(mean),mysd/mean)
colnames(y)=c("logMean","cv")
Rbs.9 <- cobs(y$logMean,y$cv, nknots=20,constraint="none",tau=0.99)
Rbs.median <- cobs(y$logMean,y$cv,nknots=20,constraint="none",tau=0.5)
pred<-data.frame(predict(Rbs.9))

y$interpolated<-interpolate(y,pred)
y$diff=y$cv-y$interpolated
yy <- y[which(y$diff>0.2),]
yy <- yy[order(-yy$diff),]
yy <- head(yy,200)

write.table(yy,file=REGE3)

png(CVE3,height=960,width=960)
plot(y$logMean,y$cv,pch=18,cex=0.5,xlab="log10(mean)",ylab="CV")
lines(predict(Rbs.9), col = "red", lwd = 1.0)
lines(predict(Rbs.median), col = "blue", lwd = 1.0)
points(yy$logMean,yy$cv)
text(yy$logMean,yy$cv+0.02,labels=rownames(yy),cex=0.8)
dev.off()

pdf(HME3)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
zz<-x[which(rownames(x) %in% rownames(yy)),]
heatmap.2(as.matrix(zz),margin=c(8, 22),cexRow=0.65,trace="none",
  cexCol=0.8,col=my_palette,scale="row")

heatmap.2(cor(t(zz)),trace="none",scale="none",margins=c(12,12),
  cexRow=0.8, cexCol=0.8)
dev.off()


hr <- hclust(as.dist(1-cor(t(zz), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- rev(redgreen(75))

colfunc <- colorRampPalette(c("blue", "white", "red"))
heatmap.2(as.matrix(zz), main="Hierarchical Cluster", Rowv=as.dendrogram(hr), Colv=NA,
 dendrogram="row", scale="row", col = colfunc(15), density.info="none", trace="none",
 RowSideColors= myClusterSideBar, margin=c(8, 22) )
dev.off()

write.table(mycl,file=HCE3,quote=F,sep="\t")
}
