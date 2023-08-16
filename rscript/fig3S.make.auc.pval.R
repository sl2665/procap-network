source('rscript/scatterdist.R')

# Load tables
if(F) {
pro = read.table('readcount/procap.var.norm.txt',stringsAsFactors=F)
pro.near = read.table('window/procap.window.all-all.1M.txt',stringsAsFactors=F)
}


# Function to generate correlation coefficient plot for RNA-seq vs nTSS 
# returns the quantile plots of the c.c.'s
cor.pro = function(window.list, table1, table2, file = "proseq.cor.pdf", asym=F) {
	pro.scd = scatter.dist(table1, table2, window.list)

	x = unlist(pro.scd$dist)
	y = unlist(pro.scd$corr)

	x = x[is.finite(x)]
	y = y[is.finite(y)]
	if(asym) x=abs(x)

	cdc = quantile.curve(x,y,n=1000,p.val=c(0.05,0.5,0.95),spar=0)
	return(list(cdc=cdc$r, x=x, y=y, scd=pro.scd))
}

# Function to reselect cor.coefs from pre calculated ccs
recor.pro = function(window.list, window.ref, ref.cor) {
	corr = sapply(1:length(window.list), function(i) {
		if(!is.null(window.list[[i]])) 
			return(ref.cor$scd$corr[[i]][match(window.list[[i]],window.ref[[i]])])
		else
			return(NULL)
	})
	dist = sapply(1:length(window.list), function(i) {
		if(!is.null(window.list[[i]])) 
			return(ref.cor$scd$dist[[i]][match(window.list[[i]],window.ref[[i]])])
		else
			return(NULL)
	})
	x = unlist(dist)
	y = unlist(corr)
	x = x[is.finite(x)]
	y = y[is.finite(y)]

	cdc = quantile.curve(x,y,n=1000,p.val=c(0.05,0.5,0.95),spar=0)
	return(list(cdc=cdc$r, x=x, y=y, scd=list(dist=dist, corr=corr)))
}

# Function to make window list
make.window.list = function(window.table) {
	return( apply(window.table, 1, function(x) {
		if(x[3] <= x[1]) return(NULL)
		else return( (x[1]+1):x[3] ) } ))
}

# Function to exclude pairs crossing the sites from the window list 
exclude.window = function(window, pos1, pos2, pos, n.cross=1) {
	conv.pos = function(p) {
		chr = substring(p[,1],4)
		chr[chr=="X"] = "23"
		chr[chr=="Y"] = "24"
		chr = as.numeric(chr)
		chr[is.na(chr)] = 0
		return(chr*1000000000 + p[,2])
	}
	p1 = conv.pos(pos1)
	p2 = conv.pos(pos2)
	p = conv.pos(pos)
	p = sort(p)
	iv = findInterval(p1, p)
	n = length(window)
	n.pos = length(p)
	return( sapply(1:n, function(i) {
		dn.ptr = iv[i] + n.cross
		up.ptr = iv[i] - n.cross + 1
		if(dn.ptr > n.pos) dn.ptr = n.pos
		if(up.ptr < 1) up.ptr = 1
		dn.pos = p[dn.ptr]
		up.pos = p[up.ptr]
		rp2 = p2[window[[i]]]
		return(window[[i]][rp2 > up.pos & rp2 < dn.pos])
	}) )
}

overlap.window = function(window, pos1, pos2, pos, n.over=1, range=500) {
	conv.pos = function(p) {
		chr = substring(p[,1],4)
		chr[chr=="X"] = "23"
		chr[chr=="Y"] = "24"
		chr = as.numeric(chr)
		chr[is.na(chr)] = 0
		return(chr*1000000000 + p[,2])
	}
	p1 = conv.pos(pos1)
	p2 = conv.pos(pos2)
	p = conv.pos(pos)
	p = sort(p)
	iv1 = findInterval(p1, p)
	iv2 = findInterval(p2, p)
	n.p1 = length(p1)
	n.p2 = length(p2)
	n = length(p)
	t1 = c()
	t2 = c()
	for(i in 1:n.p1) {
		up.ptr = iv1[i]
		dn.ptr = iv1[i] + 1
		if(up.ptr < 1) up.ptr = 1
		if(dn.ptr > n) dn.ptr = n
		if(abs(p1[i]-p[up.ptr]) < range | abs(p[dn.ptr]-p1[i]) < range) t1[i] = T
		else t1[i] = F
	}

	for(i in 1:n.p2) {
		up.ptr = iv2[i]
		dn.ptr = iv2[i] + 1
		if(up.ptr < 1) up.ptr = 1
		if(dn.ptr > n) dn.ptr = n
		if(abs(p2[i]-p[up.ptr]) < range | abs(p[dn.ptr]-p2[i]) < range) t2[i] = T
		else t2[i] = F
	}

	return( sapply(1:n.p1, function(i) {
		rp1 = t1[i]
		rp2 = t2[window[[i]]]
		return(window[[i]][rp1+rp2 >= n.over])
	} ) )
}

# Function to subtract window lists
diff.window = function(window1, window2) {
	n = length(window1)
	return( sapply(1:n, function(i) {
		return(setdiff(window1[[i]], window2[[i]]))
	}) )
}

# Function to generate a correlation decay plot for ctcf(or other intersecting sites) 
cdc.intersect = function(window, pos, ref.cor, cross=5) {
	window.in = exclude.window(window, pro[,1:2], pro[,1:2], pos, n.cross=cross)
	window.ex = diff.window(window, window.in)

	pro.in.cor = recor.pro(window.in, window, ref.cor)
	pro.ex.cor = recor.pro(window.ex, window, ref.cor)

	cdc.in = list()
	cdc.ex = list()
	for(i in 1:1) {
		cdc.in[[i]] = pro.in.cor$cdc[,c(1,4)]	
		cdc.ex[[i]] = pro.ex.cor$cdc[,c(1,4)]
	}
	return(list(cdc.in = cdc.in, cdc.ex = cdc.ex))
}

cdc.overlap = function(window, pos, ref.cor) {
	window.o1 = overlap.window(window, pro[,1:2], pro[,1:2], pos, n.over=1)
	window.o2 = overlap.window(window, pro[,1:2], pro[,1:2], pos, n.over=2)
	window.o0 = diff.window(window, window.o1)
	
	pro.o0.cor = recor.pro(window.o0, window, ref.cor)
	pro.o1.cor = recor.pro(window.o1, window, ref.cor)
	pro.o2.cor = recor.pro(window.o2, window, ref.cor)

	cdc.list = list(pro.o0.cor$cdc, pro.o1.cor$cdc, pro.o2.cor$cdc)
	return(cdc.list)
}

auc = function(cdc.element, xlim=c(0,200)) {
	n = nrow(cdc.element)
	x = c(cdc.element[seq(1,n-1,by=2),1],cdc.element[n,1]) / 1000
	y = cdc.element[seq(1,n-1,by=2),2]
	ulim = max(which(x <= xlim[2]))
	llim = min(which(x >= xlim[1]))
	auc = sum((x[(llim+1):ulim] - x[llim:(ulim-1)]) * y[llim:(ulim-1)])
	if(llim>1) auc = auc + (x[llim]-xlim[1]) * y[llim-1]
	if(ulim<n/2) auc = auc + (xlim[2]-x[ulim]) * y[ulim]
	return(as.numeric(auc))
}

auc.cdc = function(cdc.list, cdc.in=1, cdc.ex=1, xlim=c(0,200)) {
	c.in = cdc.list$cdc.in[[cdc.in]]
	c.ex = cdc.list$cdc.ex[[cdc.ex]]
	return(c(auc(c.in,xlim=xlim), auc(c.ex,xlim=xlim)))
}

auc.cdc.over = function(cdc.list, xlim=c(0,200)) {
	return(c(auc(cdc.list[[1]][,c(1,4)],xlim=xlim), auc(cdc.list[[2]][,c(1,4)],xlim=xlim),
		auc(cdc.list[[3]][,c(1,4)],xlim=xlim)))
}

# End of function definition 
####################################
if(F) {
# Calculate correlation scatterplot along the distance
window = make.window.list(pro.near)
pro.cor  = cor.pro(window, pro, pro)

# Load random controls

rand.control = c("05k","10k","15k","20k","25k","30k","40k","50k","60k","70k")
rand.control.count = c(5,10,15,20,25,30,40,50,60,70)
auc.res2 = list()

for(i in 3:10) {
	auc.res2[[rand.control[i]]] = matrix(ncol=5,nrow=100)
	colnames(auc.res2[[rand.control[i]]]) = c("i1","i2","o0","o1","o2")
	for(j in 1:100) {
		tflist = read.table(paste("AUC_pval_estimation/rand",rand.control[i],"/shuf.",rand.control[i],".",formatC(j,width=3,format="d",flag="0"),".txt",sep=""), stringsAsFactors=F)
		#cdc.list = cdc.intersect(window, tflist, pro.cor, cross=2)
		cdc.ol.list = cdc.overlap(window, tflist, pro.cor)
		#auc.val = c(auc.cdc(cdc.list),auc.cdc.over(cdc.ol.list))
		auc.val = c(auc.cdc.over(cdc.ol.list))
		print(auc.val)
		auc.res2[[rand.control[i]]][j,3:5] = auc.val
	}
}
}


pdf("pdf/Fig3/Fig3S.ref.dAUC.pdf",width=4,height=3)
par(mar=c(3,3,1,1),mgp=c(1.5,0.4,0))
plot(0,0,xlim=c(0,80),ylim=c(-50,50),axes=FALSE,xlab='Number of random TF binding sites (k)',ylab=expression(Delta*AUC),xaxs='i',yaxs='i',type='n')
for(i in 1:10) points(rand.control.count[i]+runif(100,-1,1),auc.res2[[i]][,1]-auc.res2[[i]][,2],pch=16,cex=0.3)
auc.mean = unlist(lapply(auc.res2, function(x) mean(x[,1]-x[,2])))
auc.sd = unlist(lapply(auc.res2, function(x) sd(x[,1]-x[,2])))
lines(rand.control.count,auc.mean)
lines(rand.control.count,auc.mean-auc.sd,lty=2)
lines(rand.control.count,auc.mean+auc.sd,lty=2)
axis(1,lwd=0,lwd.ticks=0.5,tck=-0.03,cex.axis=0.8)
axis(2,lwd=0,lwd.ticks=0.5,tck=-0.03,cex.axis=0.8)
legend("topleft",legend=c('mean',expression(mean%+-%sd)),lwd=1,pch=32,lty=c(1,2),bty='n',ncol=1,cex=0.8)
#grid(lty=1,lwd=0.5)
box(lwd=0.5)
dev.off()

pdf("pdf/Fig3/Fig3S.ref.dAUC.ol.pdf",width=4,height=3)
par(mar=c(3,3,1,1),mgp=c(1.5,0.4,0))
plot(0,0,xlim=c(0,80),ylim=c(-10,10),axes=FALSE,xlab='Number of random TF binding sites (k)',ylab=expression(Delta*AUC),xaxs='i',yaxs='i',type='n')
for(i in 1:10) points(rand.control.count[i]+runif(100,-1,1),auc.res2[[i]][,4]-auc.res2[[i]][,3],pch=16,cex=0.3)
auc.mean = unlist(lapply(auc.res2, function(x) mean(x[,4]-x[,3])))
auc.sd = unlist(lapply(auc.res2, function(x) sd(x[,4]-x[,3])))
lines(rand.control.count,auc.mean)
lines(rand.control.count,auc.mean-auc.sd,lty=2)
lines(rand.control.count,auc.mean+auc.sd,lty=2)
axis(1,lwd=0,lwd.ticks=0.5,tck=-0.03,cex.axis=0.8)
axis(2,lwd=0,lwd.ticks=0.5,tck=-0.03,cex.axis=0.8)
legend("topleft",legend=c('mean',expression(mean%+-%sd)),lwd=1,pch=32,lty=c(1,2),bty='n',ncol=1,cex=0.8)
#grid(lty=1,lwd=0.5)
box(lwd=0.5)
dev.off()


