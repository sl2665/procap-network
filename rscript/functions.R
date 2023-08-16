source('rscript/general/scatterplot.R')
source('rscript/general/boxplot.R')
source('rscript/general/plot.R')

# Function to generate smooth quantile curves from a scatterplot distribution along the x
quantile.curve = function(x, y, xlim=NULL, n=NULL, n.interval=400, p.val=c(.05,.25,.5,.75,.95), spar=NULL) {
	if(is.null(xlim)) {
		xlim[1] = min(x)
		xlim[2] = max(x)
	}
	if(is.null(n)) x.interval = seq(xlim[1],xlim[2],length.out=n.interval+1)
	else if(length(x) <= n) {
		x.interval = c(0,1)
		n.interval = 1
	} else {
		x.interval = quantile(x,c(seq(0,1,by = n/length(x)),1))
		n.interval = length(x.interval)-1
	}
	r = t( sapply(1:n.interval, function(i) {
		return( quantile(y[x>=x.interval[i]&x<x.interval[i+1]],prob=p.val,na.rm=T))}))

	xs = (x.interval[1:n.interval]+x.interval[1:n.interval+1])/2
	rs = r
	if(n.interval>1) {
	  for(i in n.interval:1) if(any(is.na(r[i,]))) {
	    rs = rs[-i, ]
	    xs = xs[-i]
	  }
	}
	  
	if(n.interval>4) s = cbind( xs, apply(rs,2,function(i){return(smooth.spline(xs,i,spar=spar)$y)}))
	else s = r
	r = cbind( c(x.interval[1], rep(x.interval[-c(1,n.interval+1)],each=2), x.interval[n.interval+1]),
		matrix( rep(r, each=2), ncol=length(p.val) ))
	return(list(s=s, r=r))
}

# Function to calculate p-value from Pearson correlation coefficient
cor2pvalue = function(r, n) {
  r[!is.finite(r)] = 0
  t <- (r*sqrt(n-2))/sqrt(1-r^2)
  p <- 2*(1 - pt(abs(t),(n-2)))
  return(p)
}

## Function to gnerate the lists of correlation coefficients and the distace between two tables
scatter.dist = function(table1, table2, window.list, list1=NULL, list2=NULL) {
	if(is.null(list1)) {
	  col.id1 = col.id2 = 1:(ncol(table1)-2)
  } else {
		common.id = intersect(list1, list2)
		col.id1 = match(common.id, list1)
		col.id2 = match(common.id, list2)
  }
  n = length(col.id1)
  corrfunc = function(x) {
    i = as.numeric(window.list[[x]])
    if(is.null(i)) return(NULL)	
    return(cor(t(table1[x,col.id1+2]), t(table2[i,col.id2+2])))
  }
  corr = sapply(1:nrow(table1), corrfunc)
  pval = sapply(1:nrow(table1), function(x) cor2pvalue(corr[[x]], n))
	dist = sapply(1:nrow(table1), function(x) {
		i = as.numeric(window.list[[x]])
		if(is.null(i)) return(NULL)
		return( table2[i,2] - table1[x,2])
	})
	return(list(dist=dist, corr=corr, pval = pval))	
}

# Function to plot scatterplot of correlation coefficients along the distance x
plot.scatter.dist = function(x, y, cdc=NULL, file="test.pdf", xlim=NULL, ylim=c(-1,1), xlab="Distance (kb)",
                             ylab="Correlation coefficients", p.val=c(0.05,0.5,0.95), boxplot=F,
                             col = NULL, alpha = 1, grid = F) {
	if(!boxplot) {
		custom.scatter.open(file=file,xlab=xlab,ylab=ylab,height=3,width=4)
		custom.scatter.plot(x/1000,y,xlim=xlim/1000,ylim=ylim,height=3,width=4,col=col,pval=F, alpha = alpha, lwd = 1.25, border.width =  1.75)
		if(is.null(cdc)) cdc = quantile.curve(x,y,xlim=xlim,n=1000,p.val=p.val,spar=spar)
		for(i in 1:length(p.val)) {
			lines(cdc[,1]/1000, cdc[,i+1],lwd=3,col='white')
			lines(cdc[,1]/1000, cdc[,i+1],lwd=1.25,col='black')
		}
		custom.scatter.label(xlab=xlab,ylab=ylab)
		custom.scatter.close()
	} else {
		ndata = split(y,cut(x/1000,xlim/1000))
		n = length(xlim)-1
		custom.boxplot(ndata[1:n], ylim=ylim, ylab=ylab, xlab=xlab, filename=file,
		               labels=paste(xlim[1:n]/1000,xlim[1:n+1]/1000,sep="-"), width=0.9, grid=grid, lwd=1,
		               col = col, alpha = alpha, border.width = 1.5)
	}	
}

plot.cdc.fold = function(x, y, file="test.pdf", xlim=NULL, ylim=c(-1,1), xlab="Distance (kb)",
                             ylab="Correlation coefficients", p.val=c(0.05,0.5,0.95),
                             col = NULL, grid = F, lwd = 1) {
  custom.scatter.open(file=file,xlab=xlab,ylab=ylab,height=3,width=4)
  cdc.up = quantile.curve(-x[x<0], y[x<0], xlim = c(0, -xlim[1]), n = 1000, p.val = p.val, spar = 0.1)
  cdc.dn = quantile.curve(x[x>0], y[x>0], xlim = c(0, xlim[2]), n = 1000, p.val = p.val, spar = 0.1)

  custom.scatter.plot(x,y,xlim=c(0, xlim[2]/1000),ylim=ylim,height=3,width=4,col=col,pval=F,noimage=T)
  
  for(i in 1:length(p.val)) {
    lines(cdc.up$r[,1]/1000, cdc.up$r[,i+1], lwd=1, col=col[1])
    lines(cdc.dn$r[,1]/1000, cdc.dn$r[,i+1], lwd=1, col=col[2])
  }
  par(xpd=TRUE)
  custom.scatter.label(xlab=xlab,ylab=ylab)
  custom.scatter.close()
}

# Function to plot percentage of significant correlation coefficients along the distance x
# FDR of correlations
# New function to plot percentage of correlations with FDR less than 0.1 or 0.05, and number displayed on top
source("rscript/general/plot.R")
plot.fdr.dist = function(pcor, file="pdf/test.pdf",
                         xlim=c(0,1000*2^(0:10)), ylim = c(0, 1),
                         fdr = 0.1, main = "All tTRE", width = 1,
                         grid = F, cols = NULL, alpha = 1) {
  dt = data.frame(d = unlist(pcor$scd$dist), r = unlist(pcor$scd$corr), p = unlist(pcor$scd$pval))
  dt$fdr = p.adjust(dt$p, method = "fdr")
  dt$fdr[dt$r < 0] = 1
  ndata = split(dt$fdr,cut(dt$d/1000,xlim/1000))
  f = unlist(lapply(ndata, function(x) sum(x<fdr)/length(x))) * 100
  n = unlist(lapply(ndata, function(x) length(x)))
  nxlim = length(xlim)-1
  labels = paste(xlim[1:nxlim]/1000,xlim[1:nxlim+1]/1000,sep="-")
  custom.barplot(f, filename = file, ylim=c(0, 100), ylab="Percentage", xlab="Distance (kb)",
                 labels=labels, width=1, grid=grid, lwd=1.5, cols = cols, main = main, alpha = alpha)
}


# Function to generate correlation coefficient plot for PRO-cap vs nTSS 
# returns the quantile plots of the c.c.'s
cor.pro = function(window.list, table1, table2, asym=F) {
	pro.scd = scatter.dist(table1, table2, window.list)

	x = unlist(pro.scd$dist)
	y = unlist(pro.scd$corr)

	y[!is.finite(y)]=0
	if(asym) x=abs(x)

	cdc = quantile.curve(x,y,n=1000,p.val=c(0.05,0.5,0.95),spar=0)
	return(list(cdc=cdc$r, x=x, y=y, scd=pro.scd))
}

# Function to generate correlation coefficient plot for RNA-seq vs nTSS 
# returns the quantile plots of the c.c.'s
cor.rna_pro = function(window.list, pro.table) {
	rna.scd = scatter.dist(rna, pro.table, window.list, rna.list, pro.list)
	rna.scd$dist = assign.direction(rna.scd$dist)

	x = unlist(rna.scd$dist)
	y = unlist(rna.scd$corr)

	y[!is.finite(y)]=0
	

	cdc = quantile.curve(x,y,n=1000,p.val=c(0.05,0.5,0.95),spar=0)
	return(list(cdc=cdc$r, x=x, y=y, scd=rna.scd))
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
	y[!is.finite(y)]=0
	

	cdc = quantile.curve(x,y,n=1000,p.val=c(0.05,0.5,0.95),spar=0)
	return(list(cdc=cdc$r, x=x, y=y, scd=list(dist=dist, corr=corr)))
}


# Function to generate the plot only
plot.cor = function(cor.list, file = "proseq.cor.pdf", xlim=c(0,800000), xlog=F) {
	plot.scatter.dist(cor.list$x , cor.list$y, cor.list$cdc, file = file, xlim = xlim, xlog = xlog)
}

# Function to make window list
make.window.list = function(window.table) {
	return( apply(window.table, 1, function(x) {
		if(x[3] <= x[1]) return(NULL)
		else return( (x[1]+1):x[3] ) } ))
}

make.window.asym.list = function(window.table) {
	return( apply(window.table, 1, function(x) x[2]:x[3] ))
}

assign.direction = function(dist.list) {
	return( sapply(1:length(dist.list), function(i) {
		if(rna.pos[i,6]=="+") return(dist.list[[i]])
		else return(-dist.list[[i]])
	}) )
}

make.rna.window.list = function(window.table) {
	return( apply(window.table, 1, function(x) {
		if(x[3] < x[2]) return(NULL)
		else return( x[2]:x[3] )
	}) )
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
cdc.intersect = function(window, pos, ref.cor, max.cross=5) {
	window.in = list()
	window.ex = list()
	for(i in 1:max.cross) {
		window.in[[i]] = exclude.window(window, pro[,1:2], pro[,1:2], pos, n.cross=i)
		window.ex[[i]] = diff.window(window, window.in[[i]])
	}
	pro.in.cor = list()
	pro.ex.cor = list()
	for(i in 1:max.cross) {
		pro.in.cor[[i]] = recor.pro(window.in[[i]], window, ref.cor)
		pro.ex.cor[[i]] = recor.pro(window.ex[[i]], window, ref.cor)
	}
	cdc.in = list()
	cdc.ex = list()
	cdc.im = list()
	cdc.em = list()
	for(i in 1:max.cross) {
		cdc.in[[i]] = pro.in.cor[[i]]$cdc[,c(1,4)]	
		cdc.ex[[i]] = pro.ex.cor[[i]]$cdc[,c(1,4)]
		cdc.im[[i]] = pro.in.cor[[i]]$cdc[,c(1,3)]	
		cdc.em[[i]] = pro.ex.cor[[i]]$cdc[,c(1,3)]
	}
	return(list(cdc.in = cdc.in, cdc.ex = cdc.ex,
	            cdc.im = cdc.im, cdc.em = cdc.em))
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

cal.AUC = function(x, y, max.dist = 200000) {
  n=sum(x < max.dist)
  return(sum(y[2:n] * (x[2:n] - x[1:(n-1)])) / x[n])
}

plot.cdc = function(cdc.list, cdc.in=1, cdc.ex=2, file="test.pdf", factor=NULL,
                    xlim=c(0,200), xlog=F, ylim = c(0.3, 1), do.med = F) {
	pdf(file=file, width=3.6, height=2.4, pointsize=8)
	par(mar=c(2.5,2.5,1,1), mgp=c(1.5,0.3,0))
	ylab = "Correlation coefficients (95 percentile)"
	if(do.med) ylab = "Correlation coefficients (median)" 
	if(xlog) plot(0, 0, xlim=log10(xlim+1), ylim=ylim, axes=F, xlab="Distance (kb)",
		ylab=ylab, xaxs='i', yaxs='i', type='n')
	else plot(0, 0, xlim=xlim, ylim=ylim, axes=F, xlab="Distance (kb)",
		ylab=ylab, xaxs='i', yaxs='i', type='n')
	c.in = cdc.list$cdc.in[[cdc.in]]
	c.ex = cdc.list$cdc.ex[[cdc.ex]]
	if(do.med) {
	  c.in = cdc.list$cdc.im[[cdc.in]]
	  c.ex = cdc.list$cdc.em[[cdc.ex]]
	}
	AUC.in = cal.AUC(c.in[,1], c.in[,2])
	AUC.ex = cal.AUC(c.ex[,1], c.ex[,2])

	if(xlog) {
		lines(log10(c.in[,1]/1000+1), c.in[,2], col="#084FA0")
		lines(log10(c.ex[,1]/1000+1), c.ex[,2], col="#E09810")
	}
	else {
		lines(c.in[,1]/1000, c.in[,2], col="#084FA0")
		lines(c.ex[,1]/1000, c.ex[,2], col="#E09810")
	}
	axis(1,lwd=0,lwd.ticks=0.5, tck=-0.03, cex.axis=0.75)
	axis(2,lwd=0,lwd.ticks=0.5, tck=-0.03, cex.axis=0.75)
	legend("topright",
	       legend = c(as.expression(bquote("" <= .(cdc.in-1))),as.expression(bquote("" >= .(cdc.ex)))),
	       col=c('#084FA0','#E09810'),
		lwd=1, bty='n', pch=32, ncol=1, title=paste(factor,"intersection"))
	box(lwd=0.5)
	dev.off()
	return(c(AUC.in, AUC.ex))
}

plot.cdc.over = function(cdc.list, file="test.pdf", factor=NULL, xlim=c(0,200), xlog=F, do.med=F, ylim = c(0.3,1)) {
	pdf(file=file, width=3.6,height=2.4, pointsize=8)
	par(mar=c(2.5,2.5,1,1), mgp=c(1.5,0.3,0))
	ylab = "Correlation coefficients (95 percentile)"
	if(do.med) ylab = "Correlation coefficients (median)" 
	if(xlog) plot(0, 0, xlim=log10(xlim+1), ylim=ylim, axes=F, xlab="Distance (kb)",
		ylab=ylab, xaxs='i', yaxs='i', type='n')
	else plot(0, 0, xlim=xlim, ylim=ylim, axes=F, xlab="Distance (kb)",
		ylab=ylab, xaxs='i', yaxs='i', type='n')
	cc=4
	if(do.med) cc = 3
	if(xlog) {
		lines(log10(cdc.list[[1]][,1]/1000+1), cdc.list[[1]][,cc], col="#084FA0")
		lines(log10(cdc.list[[2]][,1]/1000+1), cdc.list[[2]][,cc], col="#80C080")
		lines(log10(cdc.list[[3]][,1]/1000+1), cdc.list[[3]][,cc], col="#E09810")
	} else {
		lines(cdc.list[[1]][,1]/1000, cdc.list[[1]][,cc], col="#084FA0")
		lines(cdc.list[[2]][,1]/1000, cdc.list[[2]][,cc], col="#80C080")
		lines(cdc.list[[3]][,1]/1000, cdc.list[[3]][,cc], col="#E09810")
	}
	AUC.in = cal.AUC(cdc.list[[1]][,1], cdc.list[[1]][,cc])
	AUC.ex = cal.AUC(cdc.list[[3]][,1], cdc.list[[3]][,cc])
	axis(1,lwd=0,lwd.ticks=0.5, tck=-0.03, cex.axis=0.75)
	axis(2,lwd=0,lwd.ticks=0.5, tck=-0.03, cex.axis=0.75)
	legend("topright", legend=c(0,1,2), col=c('#084FA0','#80C080','#E09810'),
		lwd=1, bty='n', pch=32, ncol=1, title=paste(factor,"overlap"))
	box(lwd=0.5)
	dev.off()
	return(c(AUC.in, AUC.ex))
}

# Function to take a list and match the distribution by sampling elements from multiple sets and return the id of matched sets
match.list = function(x, split=100) {
  x = lapply(x, log10)
  minx = min(unlist(x))
  maxx = max(unlist(x))
  range = seq(minx,maxx,length.out=split+1)
  xs = lapply(x, function(a) split(seq_along(a), cut(a,range,include.lowest=T)))
  len = matrix(unlist(lapply(xs,function(a) lapply(a,length))),ncol=length(x))
  mlen = apply(len,1,min)
  id = list()
  for(i in 1:length(x)) id[[i]] = sapply(1:split,function(a) sample(xs[[i]][[a]],mlen[a]))
  return(lapply(id, unlist, use.names=F))
}
