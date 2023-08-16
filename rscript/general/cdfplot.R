custom.cdf = function(data, file='test.pdf', legend=NULL, legend.title=NULL, xlab=NULL, ylab="Cumulative fraction", main=NULL, xlim=NULL, ylim=NULL, col=NULL, lwd=2, height=2, width=3.236, add=F, step=NULL, text=NULL)
{
    n = length(data)
    ndata = lapply(data, function(x) return(sort(x[is.finite(x)])))
    
    mar=c(3,3,3,1)
    if(is.null(xlab)) mar[1]=1.5
    if(is.null(ylab)) mar[2]=1.5
    if(is.null(main)) mar[3]=1.5
    
    if(!add) pdf(file,width=0.2*(mar[2]+mar[4])+width,height= 0.2*(mar[1]+mar[3])+height)
    par(mar=mar, mgp=c(1.5,0.3,0))
    
    if(is.null(col)) {
		if(n==2) col = colorRampPalette( c('#201040', '#302070', '#104090', '#0060b0', '#00a0a0', '#80c080', '#d0c010', '#f0c010', '#d07010', '#b02030', '#701020'))(5)[c(2,4)]
		else col = colorRampPalette( c('#201040', '#302070', '#104090', '#0060b0', '#00a0a0', '#80c080', '#d0c010', '#f0c010', '#d07010', '#b02030', '#701020'))(n+2)[1:n+1]
	}

    if(is.null(xlim)) xlim = quantile(unlist(ndata),c(0.005,0.995))
    xrange = xlim[2] - xlim[1]
    if(is.null(ylim)) ylim = c(0, 1)
    plot(0,0,xlim=xlim,ylim=ylim,axes=FALSE,xlab=xlab,ylab=ylab,main=main,type='n')

    for(i in 1:n) {
        N = length(ndata[[i]])
        if(N>1) {
			if(is.null(step))
				lines(c(xlim[1]-xrange,rep(ndata[[i]],each=2),xlim[2]+xrange),c(0,0,rep(1:(N-1),each=2),N,N)/N, lwd=lwd, col=col[i])
			else {
				j = floor(seq(1,N,length.out=step))
				lines(c(xlim[1]-xrange,ndata[[i]][j],xlim[2]+xrange),c(0,(j-0.5)/N,1), lwd=lwd, col=col[i])
			}
		}
    }

    axis(1,tck=-0.03,cex.axis=0.8,lwd=0,lwd.tick=0.5)
    axis(2,tck=-0.03,cex.axis=0.8,lwd=0,lwd.tick=0.5)
    #grid(lty=1,lwd=0.5,col='#a0a0a0')
    box(lwd=0.5)
    if(!is.null(legend)) legend("bottomright", legend=legend, col=col, lwd=lwd, pch=32, bty='n', ncol=1, cex=0.8, title=legend.title)
	if(!is.null(text)) text(par("usr")[1], par("usr")[4], text, cex=0.8, adj=c(-0.1,1.3))
    if(!add) dev.off()
}

custom.cdfD = function(x, y, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, labels=NULL, file='test.pdf', n=4, q=NULL, height=2, width=3.236, col=NULL, main=NULL, lwd=2, add=F)
{
    l=is.finite(x)&is.finite(y)
    xs=order(x[l])
    ys=(y[l])[xs]
    if(!is.null(q)) ndata = split(ys,cut(seq_along(ys),q*length(ys)))
    else ndata=split(ys,ceiling(seq_along(ys)/(length(ys)/n)))

    p = bquote(p==.(format(ks.test(ndata[[1]],ndata[[length(ndata)]])$p.value,digits=2)))
    custom.cdf(ndata, xlim=xlim, ylim=ylim, ylab='Cumulative fraction', xlab=ylab, col=col, main=main,
		file=file, height=height, width=width, legend=labels, legend.title=xlab, lwd=lwd, add=add, text=p)
}
