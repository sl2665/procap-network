source("~/Documents/data/rscript/heatmap.R")

## Make arbitrary matrix
#mat = matrix(runif(60),ncol=2)-0.5
n = nrow(mat)
#rownames(mat) = Intersect$X1

# Value range
zlim = c(-15,15)

# Scalebar labels
at = c(-10, 0, 10)

# Hierachical clustering
cl = hclust(dist(mat), method="ward.D")
dn = as.dendrogram(cl)
ch = max(cl$height)

# Open pdf file
width=2
height=6
pdf("Fig9.pdf",width=width,height=height)
par(mar=c(3,3,3,1),mgp=c(1.5,0.2,0))

# Ratio between dendrogram and heatmap widths
rw = 5
# Ratio between scalebar and heatmap heights
rh = 20

# Plot dendrogram
par(yaxs='i',lwd=0.5)
plot(dn, yaxt='n', horiz=T, leaflab="none", xlim=c(-rw*ch,ch), ylim=c(-n/rh*2+0.5,n+0.5))

# Plot heatmap
vmat = mat[unlist(dn)[n:1],]
if(nrow(vmat)<height*100) {
	rep=floor(height*100/nrow(vmat))
	vmat=matrix(rep(vmat,each=rep),nrow=nrow(vmat)*rep)
}
if(ncol(vmat)<width*100) {
	rep=floor(width*100/ncol(vmat))
	vmat=t(matrix(rep(t(vmat),each=rep),nrow=ncol(vmat)*rep))
}

heatmap.1col(vmat, xlim=c(-rw*ch,0), ylim=c(0.5,n+0.5), zlim=zlim)

# Plot scalebar
heatmap.1col(matrix(0:100/100,nrow=1),xlim=c(-rw*ch,0),ylim=c(-n/rh*2+0.5,-n/rh+0.5))

# Plot scalebar axis
axis(1,lwd=0,tck=-0.02,lwd.ticks=0.5,at=((at-zlim[1])/(zlim[2]-zlim[1])-1)*rw*ch,labels=at,cex.axis=0.6)

# Print sample names
text(-rw*ch*1.05,1:n,labels=rownames(mat)[unlist(dn)],xpd=T,adj=1,cex=0.4)

dev.off()
