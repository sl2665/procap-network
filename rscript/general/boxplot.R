custom.boxplot<-function(data, labels=NULL, ylim=NULL, ylab=NULL, xlab=NULL, main=NULL, filename='test.pdf',
                         cols=NULL, alpha = 1, paired = F, border.width = 1,
                         width=1, lwd=0.75, med.lwd = 2, outlier=F, staple=T, close=T, grid=F, xaxt=T, yaxt=T, at=NULL, boxwd=0.75)
{
    n=length(data)
    mar=c(3,3,4,1)
    if(is.null(labels)) mar[1]=1
    if(is.null(xlab)) mar[3]=mar[3]-1.5
    if(is.null(ylab)) mar[2]=1.5
    if(is.null(main)) mar[3]=mar[3]-1.5

    if(dev.cur()==1) {
		pdf(filename,width=0.2*(mar[2]+mar[4])+0.25*n*width,height= 2+0.2*(mar[1]+mar[3]))
		par(mar=mar)
	}
    par(mgp=c(1.5,0.3,0))
    
	ndata=list()
    for(i in 1:n) {ndata[[i]]=data[[i]][is.finite(data[[i]])]}
   if(paired) n = n/2 
	 if(is.null(cols)) cols=rep('#ffffff',n)
    else if(cols[1]=='br') cols=colorRampPalette(c('#4393c3','#92c5de','#d1e5f0','#f7f7f7','#fddbc7','#f4a582','#d6604d'))(n)
    else if(cols[1]=='g') cols=colorRampPalette(c('#ffffff','#808080'))(n)
    else if(cols[1]=='rb') cols=colorRampPalette(c('#4393c3','#92c5de','#d1e5f0','#f7f7f7','#fddbc7','#f4a582','#d6604d'))(n)[n:1]
    else if(cols[1]=='p') cols=colorRampPalette(c('#352187','#0f5cdd','#1481d6','#06a4ca','#2eb7a4','#87bf77','#d1bb59','#fec832','#f9fb0e'))(n)
    else if(cols[1]=='jet') cols=colorRampPalette(c('#201040','#302070','#104090','#0060b0','#00a0a0','#80c080','#d0c010','#f0c010','#d07010','#b02030','#701020'))(n+2)[1:n+1]
    else cols=colorRampPalette(cols)(n)
	if(paired) {
	  cols = c(rbind(rep("#ffffff", n), cols))
	  n = n*2
	}
	cols = alpha(cols, alpha)
	if(is.null(ylim)) ylim=c(min(unlist(ndata)),max(unlist(ndata)))
	plot(0,0,xlim=c(0.5,n+0.5),ylim=ylim,type='n',xaxt='n',yaxt='n',xlab=NA,ylab=ylab,bty='n')
	if(grid) grid(nx=NA,ny=NULL,lwd=0.5,lty=1,col="#C0C0C0")
    boxplot(ndata,ylim=ylim,xaxt='n',pch=16,cex=0.5,main=main,ylab=ylab,cex.lab=0.8,cex.axis=0.8,col=cols,
		yaxt='n',staplewex=0.5*staple,boxwex=boxwd,boxlwd=lwd,whisklty=1,whisklwd=lwd,
		staplelwd=lwd,medlwd=med.lwd,outpch=16,outcex=0.2,outcol='black',outbg='black',
		outline=outlier,frame=F,add=T, at=at)
    if(xaxt) axis(1,at=1:(n+1)-0.5,lab=NA,lwd=0,lwd.ticks=0.5,cex.axis=0.8,tck=0.03)
    if(yaxt) axis(2,cex.axis=0.7,lwd=0,lwd.ticks=0.5,tck=-0.03)
    if(is.list(labels)) {
      n.lab = length(labels)
      for(i in 1:n.lab) text(labels=labels[[i]],x=i+0.5,y=ylim[1]-(ylim[2]-ylim[1])*0.1,srt=45,xpd=TRUE,cex=0.8,pos=2)
    } else {
      text(labels=labels,x=1:n+0.5,y=ylim[1]-(ylim[2]-ylim[1])*0.1,srt=45,xpd=TRUE,cex=0.8,pos=2)
    }
	box(lwd=border.width)
    mtext(xlab,3,0.5,cex=0.8)
    if(close) dev.off()
}

custom.boxplotD<-function(x, y, ylim=NULL, labels=c('Lowest','Lower','Higher','Highest'), xlab=NULL, ylab=NULL, filename='test.pdf', n=4, q=NULL, width=1, cols=NULL, outlier=F, main=NULL)
{
    l=is.finite(x)&is.finite(y)
    xs=order(x[l])
    ys=(y[l])[xs]
    if(!is.null(q)) ndata = split(ys,cut(seq_along(ys),q*length(ys)))
    else ndata=split(ys,ceiling(seq_along(ys)/(length(ys)/n)))
    
    custom.boxplot(ndata, ylim=ylim, ylab=ylab, xlab=xlab, cols=cols, outlier=outlier, main=main, filename=filename, width=width, labels=labels)
}


