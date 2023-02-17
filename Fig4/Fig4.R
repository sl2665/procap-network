# Read strans specific normalzied reads
pl = read.table("Table.4d.procap.ambr.norm.pl.txt",header = 1)
mn = read.table("Table.4e.procap.ambr.norm.mn.txt",header = 1)

# Distance measures between adjascent tTREs 
l = nrow(pl) - 1
dist = pl$pos[-1] - pl$pos[-(l+1)]

# Correlation coefficients between adjascent tTREs
cor.div = sapply(1:l, function(i) cor(unlist(mn[i,-(1:2)]), unlist(pl[i+1,-(1:2)])))
cor.conv = sapply(1:l, function(i) cor(unlist(pl[i,-(1:2)]), unlist(mn[i+1,-(1:2)])))
cor.pl = sapply(1:l, function(i) cor(unlist(pl[i,-(1:2)]), unlist(pl[i+1,-(1:2)])))
cor.mn = sapply(1:l, function(i) cor(unlist(mn[i,-(1:2)]), unlist(mn[i+1,-(1:2)])))

# Distance between skipping tTREs
l.s = nrow(pl) - 2
dist.s = pl$pos[-(1:2)] - pl$pos[-((l.s+1):(l.s+2))]

# Correlation coefficients between skipping tTREs
cor.div.s = sapply(1:l.s, function(i) cor(unlist(mn[i,-(1:2)]), unlist(pl[i+2,-(1:2)])))
cor.conv.s = sapply(1:l.s, function(i) cor(unlist(pl[i,-(1:2)]), unlist(mn[i+2,-(1:2)])))
cor.pl.s = sapply(1:l.s, function(i) cor(unlist(pl[i,-(1:2)]), unlist(pl[i+2,-(1:2)])))
cor.mn.s = sapply(1:l.s, function(i) cor(unlist(mn[i,-(1:2)]), unlist(mn[i+2,-(1:2)])))

# Distance categories
d.o = order(dist)
d.o.s = order(dist.s)
d.o = d.o[dist[d.o]>250]
d.o.s = d.o.s[dist.s[d.o.s]>250]

d.s = split(d.o, ceiling(seq_along(d.o)/1000))
d.s.s = split(d.o.s, ceiling(seq_along(d.o.s)/1000))

# quantiles
q = c(0.05, 0.25, 0.5, 0.75, 0.95)
cor.div.q = t(matrix(unlist(lapply(d.s, function(x) quantile(cor.div[x], q, na.rm=T))), nrow=5))
cor.conv.q = t(matrix(unlist(lapply(d.s, function(x) quantile(cor.conv[x], q, na.rm=T))), nrow=5))
cor.ss.q = t(matrix(unlist(lapply(d.s, function(x) quantile(c(cor.pl[x],cor.mn[x]), q, na.rm=T))), nrow=5))
dist.q = t(sapply(d.s, function(x) c(dist[x[1]], dist[x[length(x)]])))

cor.div.q.s = t(matrix(unlist(lapply(d.s.s, function(x) quantile(cor.div.s[x], q, na.rm=T))), nrow=5))
cor.conv.q.s = t(matrix(unlist(lapply(d.s.s, function(x) quantile(cor.conv.s[x], q, na.rm=T))), nrow=5))
cor.ss.q.s = t(matrix(unlist(lapply(d.s.s, function(x) quantile(c(cor.pl.s[x],cor.mn.s[x]), q, na.rm=T))), nrow=5))
dist.q.s = t(sapply(d.s.s, function(x) c(dist.s[x[1]], dist.s[x[length(x)]])))

plot.q = function(dist.q, cor, xdev=0) {
	lines( unlist(t(dist.q))+xdev, rep(cor, each=2) )
}


#pdf("FigS7a.pdf")
#plot(0,0, xlim=c(500,10000), ylim=c(0.1,0.4), xaxs='i', type='n')
#plot.q(dist.q, cor.conv.q[,3], -250)
#plot.q(dist.q, cor.div.q[,3], 250)
#plot.q(dist.q, cor.ss.q[,3])
#dev.off()

source('~/Sandbox/procap/rscript/scatterplot.R')

custom.scatter.open(file="FigS7a.pdf",xlab="Distance (bp)",ylab="Correlation coefficients",height=3,width=4)
d.co = dist > 500 & dist < 10000
cor.ss = c(cor.pl[d.co],cor.mn[d.co])
custom.scatter.plot(rep(dist[d.co],2), cor.ss, xlim=c(500,10000), color='jet', pval=F, height=3, width=4, lwd=100/37)
loess.ss = loess.smooth(rep(dist[d.co],2), cor.ss, span=0.1, eval=200)
lines(loess.ss$x, loess.ss$y, lwd=2*100/37, col='white')
lines(loess.ss$x, loess.ss$y, lwd=100/37, col='black')
custom.scatter.label(xlab="Distance (bp)", ylab="Correlation coefficients")
custom.scatter.close()

custom.scatter.open(file="FigS7b.pdf",xlab="Distance (bp)",ylab="Correlation coefficients",height=3,width=4)
custom.scatter.plot(dist[d.co], cor.conv[d.co], xlim=c(500,10000), color='jet', pval=F, height=3, width=4, lwd=100/37)
loess.conv = loess.smooth(dist[d.co], cor.conv[d.co], span=0.1, eval=200)
lines(loess.conv$x, loess.conv$y, lwd=2*100/37/3*4, col='white')
lines(loess.conv$x, loess.conv$y, lwd=100/37/3*4, col='black')
custom.scatter.label(xlab="", ylab="")
custom.scatter.close()

custom.scatter.open(file="FigS7c.pdf",xlab="Distance (bp)",ylab="Correlation coefficients",height=3,width=4)
custom.scatter.plot(dist[d.co], cor.div[d.co], xlim=c(500,10000), color='jet', pval=F, height=3, width=4, lwd=100/37)
loess.div = loess.smooth(dist[d.co], cor.div[d.co], span=0.1, eval=200)
lines(loess.div$x, loess.div$y, lwd=2*100/37/3*4, col='white')
lines(loess.div$x, loess.div$y, lwd=100/37/3*4, col='black')
custom.scatter.label(xlab="", ylab="")
custom.scatter.close()

pdf("FigS7d.pdf",height=1.4, width=2)
par(mar=c(1,1,0.5,0.5), mgp=c(0.5,0.2,0))
plot(0,xlim=c(500,5000), ylim=c(0.1,0.4), axes=F, xlab='', ylab='', xaxs='i')
lines(loess.ss$x, loess.ss$y, lwd=2, col='#7094CD')
lines(loess.conv$x, loess.conv$y, lwd=2, col='#80C080')
lines(loess.div$x, loess.div$y, lwd=2, col='#CC72AE')
legend("topright",legend=c("Sense","Convergent","Divergent"),col=c('#7094CD','#80C080','#CC72AE'),lwd=2, bty='n',cex=7/12)
axis(1,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5, mgp=c(0.5,-0.1,0))
axis(2,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5)
box()
dev.off()

pdf("FigS7e.pdf",height=1.4, width=2)
par(mar=c(1,1,0.5,0.5), mgp=c(0.5,0.2,0))
plot(0,xlim=c(500,5000), ylim=c(0.1,0.4), axes=F, xlab='', ylab='', xaxs='i')
lines(loess.ss$x-250, loess.ss$y, lwd=2, col='#7094CD')
lines(loess.conv$x-500, loess.conv$y, lwd=2, col='#80C080')
lines(loess.div$x, loess.div$y, lwd=2, col='#CC72AE')
legend("topright",legend=c("Sense","Convergent","Divergent"),col=c('#7094CD','#80C080','#CC72AE'),lwd=2, bty='n',cex=7/12)
axis(1,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5, mgp=c(0.5,-0.1,0))
axis(2,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5)
box()
dev.off()

# Loess for skipping
d.co.s = dist.s > 500 & dist.s < 10000
cor.ss.s = c(cor.pl.s[d.co.s],cor.mn.s[d.co.s])
loess.ss.s = loess.smooth(rep(dist.s[d.co.s],2), cor.ss.s, span=0.1, eval=200)
loess.conv.s = loess.smooth(dist.s[d.co.s], cor.conv.s[d.co.s], span=0.1, eval=200)
loess.div.s = loess.smooth(dist.s[d.co.s], cor.div.s[d.co.s], span=0.1, eval=200)

pdf("FigS7f.pdf",height=1.4, width=1.6)
par(mar=c(1,1,0.5,0.5), mgp=c(0.5,0.2,0))
plot(0,xlim=c(500,5000), ylim=c(0.1,0.4), axes=F, xlab='', ylab='', xaxs='i')
lines(loess.ss$x, loess.ss$y, lwd=2, col='#7094CD')
lines(loess.ss.s$x, loess.ss.s$y, lwd=2, col='#7094CD', lty=5)
legend("topright",legend=c("Neighboring","Interlacing"),col=c('#7094CD','#7094CD'),lwd=2, bty='n',cex=7/12, lty=c(1,2))
axis(1,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5, mgp=c(0.5,-0.1,0))
axis(2,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5)
box()
dev.off()

pdf("FigS7g.pdf",height=1.4, width=1.6)
par(mar=c(1,1,0.5,0.5), mgp=c(0.5,0.2,0))
plot(0,xlim=c(500,5000), ylim=c(0.1,0.4), axes=F, xlab='', ylab='', xaxs='i')
lines(loess.conv$x, loess.conv$y, lwd=2, col='#80C080')
lines(loess.conv.s$x, loess.conv.s$y, lwd=2, col='#80C080', lty=5)
legend("topright",legend=c("Neighboring","Interlacing"),col=c('#80C080','#80C080'),lwd=2, bty='n',cex=7/12, lty=c(1,2))
axis(1,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5, mgp=c(0.5,-0.1,0))
axis(2,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5)
box()
dev.off()

pdf("FigS7h.pdf",height=1.4, width=1.6)
par(mar=c(1,1,0.5,0.5), mgp=c(0.5,0.2,0))
plot(0,xlim=c(500,5000), ylim=c(0.1,0.4), axes=F, xlab='', ylab='', xaxs='i')
lines(loess.div$x, loess.div$y, lwd=2, col='#CC72AE')
lines(loess.div.s$x, loess.div.s$y, lwd=2, col='#CC72AE', lty=5)
legend("topright",legend=c("Neighboring","Interlacing"),col=c('#CC72AE','#CC72AE'),lwd=2, bty='n',cex=7/12, lty=c(1,2))
axis(1,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5, mgp=c(0.5,-0.1,0))
axis(2,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5)
box()
dev.off()

