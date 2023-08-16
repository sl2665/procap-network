jetc=colorRampPalette(c('#201040','#302070','#104090','#0060b0','#00a0a0','#80c080','#d0c010','#f0c010','#d07010','#b02030','#701020'))

lines.histo=function(h,col="black",lwd=1,lty=1)
{
	x=rep(h$breaks,each=2)
	x=x[-c(1,length(x))]
	y=rep(h$counts,each=2)
	lines(x,y,col=col,lwd=lwd,lty=lty)
}

lines.histo.d=function(h,col="black",lwd=1,lty=1)
{
	x=rep(h$breaks,each=2)
	x=x[-c(1,length(x))]
	y=rep(h$density,each=2)
	lines(x,y,col=col,lwd=lwd,lty=lty)
}
