custom.meta=function(data,filename='test.pdf',xlim=NULL,ylim=NULL,xlab='',ylab='',col=NULL,v=NULL,h=NULL,lwd=2)
{
	pdf(filename,height=3,width=4)
	par(mar=c(3,3,1,1),mgp=c(1.5,0.3,0))
	
	custom.metapanel(data,xlim,ylim,col,v,h,lwd)
	
	axis(1,lwd=0,lwd.ticks=1,tck=-0.02,cex.axis=0.6)
	axis(2,lwd=0,lwd.ticks=1,tck=-0.02,cex.axis=0.6)
	dev.off()
}

custom.metapanel=function(data,xlim=NULL,ylim=NULL,col=NULL,v=NULL,h=NULL,lwd=1.5)
{
	nc=ncol(data)
	nr=nrow(data)
	if(is.null(xlim)) xlim=c(data[1,1],data[nr,1])
	if(is.null(ylim)) ylim=c(min(data[,2:nc]),max(data[,2:nc]))
	plot(0,type='n',xlim=xlim,ylim=ylim,xaxt='n',yaxt='n',bty='n')
	if(!is.null(v)) abline(v=v,col='#d0d0d0')
	if(!is.null(h)) abline(h=h,col='#d0d0d0')

	if(is.null(col)) {
		if(nc==3) col=c('#b2182b','#2166ac')
		else col=rep('#000000',nc-1)
	}
	for(i in 2:nc) lines(data[1:nr,1],data[1:nr,i],col=col[i-1],lwd=lwd)
	box(lwd=0.5)
}

