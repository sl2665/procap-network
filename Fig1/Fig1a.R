source('~/Sandbox/procap/rscript/heatmap.R')

# Read table
pro = read.table('readcount/procap.ambr.norm.txt',stringsAsFactors=F)
pro.var = read.table('readcount/procap.var.norm.txt',stringsAsFactors=F)

# Locate targets
bcl2 = c(60790579, 60986613)
bcl2 = c(bcl2[1]*1.1-bcl2[2]*0.1,bcl2[2]*1.1-bcl2[1]*0.1)
slfn5 = c(33568000, 33572000)

bcl2.range = pro[,1]=="chr18" & pro[,2]>bcl2[1] & pro[,2]<bcl2[2]
slfn5.range = pro[,1]=="chr17" & pro[,2]>slfn5[1] & pro[,2]<slfn5[2]

slfn5.mat = pro[slfn5.range,-(1:2)]
bcl2.mat = pro[bcl2.range,-(1:2)]

# Make a multipanel scatterplot for slfn5
pdf("Fig6A.pdf",width=3.6,height=3.6)
par(mfrow=c(5,5))
par(cex=0.6)
par(mar=c(0.1,0.1,0.1,0.1))
par(oma=c(3,4,2,1))
par(tcl=-0.25)
par(mgp=c(1.5,.2,0))
for(i in 1:5) for(j in 1:5) {
	plot(log10(unlist(slfn5.mat[i,])+1),log10(unlist(slfn5.mat[j,])+1),pch=16,cex=0.5,xlab=NULL,ylab=NULL,axes=F)
	box(lwd=0.5)
}
dev.off()

# Make a correlation heatmap for bcl2
pdf("Fig6B.pdf",width=3.6,height=4)
par(mar=c(3,2,2,1),mgp=c(1.5,0.2,0))
plot(0,0,type='n',xlim=bcl2,ylim=bcl2,axes=F,xaxs='i',yaxs='i',xlab="",ylab="")
# Lines from TREs to heatmap
n = nrow(bcl2.mat)
pos = pro[bcl2.range,2]
heatmappos = (1:n-0.5)*(bcl2[2]-bcl2[1])/n + bcl2[1]
range = bcl2[2] - bcl2[1]
for(i in 1:n) lines(c(pos[i],heatmappos[i]),c(bcl2[2]+0.05*range,bcl2[2]),lwd=0.5,xpd=T)
for(i in 1:n) lines(c(bcl2[1]-0.05*range,bcl2[1]),bcl2[1]+bcl2[2]-c(pos[i],heatmappos[i]),lwd=0.5,xpd=T)
# Heatmap
cor.mat = cor(t(bcl2.mat))
for(i in 1:n) cor.mat[i,i] = 0
cor.mat = matrix(rep(cor.mat,each=20),ncol=ncol(cor.mat))
cor.mat = matrix(apply(cor.mat,2,rep,20),nrow=nrow(cor.mat))
heatmap.1col(cor.mat,xlim=bcl2,ylim=bcl2,zlim=c(-1,1))
heatmap.1col(matrix(-100:100/100,nrow=1),xlim=c(bcl2[2]-range*0.3,bcl2[2]),ylim=c(bcl2[1]-range*0.08,bcl2[1]-range*0.04),zlim=c(-1,1))
dev.off()
