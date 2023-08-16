source('rscript/scatterplot.R')

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
	}
	else {
		x.interval = quantile(x,c(seq(0,1,by = n/length(x)),1))
		n.interval = length(x.interval)-1
	}
	r = t( sapply(1:n.interval, function(i) {
		return( quantile(y[x>=x.interval[i]&x<x.interval[i+1]],prob=p.val,na.rm=T))}))
	if(n.interval>1) for(i in 2:n.interval) if(any(is.na(r[i,]))) r[i,]=r[i-1,]
	xs = (x.interval[1:n.interval]+x.interval[1:n.interval+1])/2
	if(n.interval>4) s = cbind( xs, apply(r,2,function(i){return(smooth.spline(xs,i,spar=spar)$y)}))
	else s = r
	r = cbind( c(x.interval[1], rep(x.interval[-c(1,n.interval+1)],each=2), x.interval[n.interval+1]),
		matrix( rep(r, each=2), ncol=length(p.val) ))
	return(list(s=s, r=r))
}

# Function to gnerate the lists of correlation coefficients and the distace between two tables
scatter.dist = function(table1, table2, window.list, list1=NULL, list2=NULL) {
	if(is.null(list1)) col.id1 = col.id2 = 1:(ncol(table1)-2)
	else {
		common.id = intersect(list1, list2)
		col.id1 = match(common.id, list1)
		col.id2 = match(common.id, list2)
	}
	corr = sapply(1:nrow(table1), function(x) {
		i = as.numeric(window.list[[x]])
		if(is.null(i)) return(NULL)	
		return( cor( t(table1[x,col.id1+2]), t(table2[i,col.id2+2])) )
	})
	dist = sapply(1:nrow(table1), function(x) {
		i = as.numeric(window.list[[x]])
		if(is.null(i)) return(NULL)	
		return( table2[i,2] - table1[x,2])
	})
	return(list(dist=dist, corr=corr))	
}

# Function to plot scatterplot of correlation coefficients along the distance x
plot.scatter.dist = function(x, y, file="test.pdf", xlim=NULL, p.val=c(0.05,0.5,0.95), spar=NULL) {
	custom.scatter.open(file=file,xlab="Distance (kb)",ylab="Correlation coefficient",height=3,width=4)
	custom.scatter.plot(x/1000,y,xlim=xlim/1000,ylim=c(-1,1),height=3,width=4)
	cdc = quantile.curve(x,y,xlim=xlim,n=1000,p.val=p.val,spar=spar)
	if(spar==0) cdc = cdc$r
	else cdc = cdc$s
	for(i in 1:length(p.val)) {
		lines(cdc[,1]/1000, cdc[,i+1],lwd=3,col='white')
		lines(cdc[,1]/1000, cdc[,i+1],lwd=1,col='black')
	}
	custom.scatter.label(xlab="Distance (kb)",ylab="Correlation coefficient")
	custom.scatter.close()
	return(cdc)
}

replot.scatter.dist = function(x, y, cdc, file="test.pdf", xlim=NULL, xlog=F) {
	custom.scatter.open(file=file,xlab="Distance (kb)",ylab="Correlation coefficient",height=3,width=4)
	if(xlog) custom.scatter.plot(log10(x/1000+1),y,xlim=log10(xlim/1000+1),ylim=c(-1,1),height=3,width=4)
	else custom.scatter.plot(x/1000,y,xlim=xlim/1000,ylim=c(-1,1),height=3,width=4)
	if(xlog) for(i in 2:ncol(cdc)) {
		lines(log10(cdc[,1]/1000+1), cdc[,i],lwd=3,col='white')
		lines(log10(cdc[,1]/1000+1), cdc[,i],lwd=1,col='black')
	}
	else for(i in 2:ncol(cdc)) {
		lines(cdc[,1]/1000, cdc[,i],lwd=3,col='white')
		lines(cdc[,1]/1000, cdc[,i],lwd=1,col='black')
	}
	custom.scatter.label(xlab="Distance (kb)",ylab="Correlation coefficient")
	custom.scatter.close()
}

