custom.barplot<-function(data, labels=NULL, ylim=NULL, ylab=NULL, xlab=NULL,
                         main=NULL, filename='test.pdf', cols=NULL, width=1, lwd=1,
                         close=T, grid=F, xaxt=T, yaxt=T, at=NULL, boxwd=0.75,
                         alpha = 1)
{
  n=length(data)
  mar=c(3,3,3,1)
  if(is.null(labels)) mar[1]=1
  if(is.null(xlab)) mar[3]=mar[3]-1.5
  if(is.null(ylab)) mar[2]=1.5
  if(is.null(main)) mar[3]=mar[3]-1.5
  
  if(dev.cur()==1) {
    pdf(filename,width=0.2*(mar[2]+mar[4])+0.25*n*width,height= 2+0.2*(mar[1]+mar[3]))
    par(mar=mar)
  }
  par(mgp=c(2,0.3,0))
  
  if(is.null(cols)) cols=rep('#ffffff',n)
  else if(cols[1]=='br') cols=colorRampPalette(c('#4393c3','#92c5de','#d1e5f0','#f7f7f7','#fddbc7','#f4a582','#d6604d'))(n)
  else if(cols[1]=='g') cols=colorRampPalette(c('#ffffff','#808080'))(n)
  else if(cols[1]=='rb') cols=colorRampPalette(c('#4393c3','#92c5de','#d1e5f0','#f7f7f7','#fddbc7','#f4a582','#d6604d'))(n)[n:1]
  else if(cols[1]=='p') cols=colorRampPalette(c('#352187','#0f5cdd','#1481d6','#06a4ca','#2eb7a4','#87bf77','#d1bb59','#fec832','#f9fb0e'))(n)
  else if(cols[1]=='jet') cols=colorRampPalette(c('#201040','#302070','#104090','#0060b0','#00a0a0','#80c080','#d0c010','#f0c010','#d07010','#b02030','#701020'))(n+2)[1:n+1]
  else cols=colorRampPalette(cols)(n)
  cols = alpha(cols, alpha)
  if(is.null(ylim)) ylim=c(min(unlist(data)),max(unlist(data)))
  plot(0,0,xlim=c(1-0.4,n+0.4),ylim=ylim,type='n',xaxt='n',yaxt='n',xlab=xlab,ylab=ylab,bty='n',main = main)
  if(grid) grid(nx=NA,ny=NULL,lwd=0.5,lty=1,col="#C0C0C0")
  n=length(data)
  rect((1:n)-0.4, 0, (1:n)+0.4, data, col = cols, cex = 0.5)
  if(xaxt) axis(1,at=1:n+0.5,lab=NA,lwd=0,lwd.ticks=0.5,cex.axis=0.8,tck=0.03)
  if(yaxt) axis(2,cex.axis=0.7,lwd=0,lwd.ticks=0.5,tck=-0.03)
  text(labels=labels,x=1:n+0.5,y=rep(ylim[1]-(ylim[2]-ylim[1])*0.075,n),srt=45,xpd=TRUE,cex=0.75,pos=2)
  box(lwd=0.5)
  #mtext(xlab,3,0.5,cex=0.8)
  if(close) dev.off()
}

plot_holder = function() {
  pdf("test.pdf",width=4,height=3)
  par(mar=c(3,3,1,1),mgp=c(1.5,0.4,0))
  custom.color=colorRampPalette(c('#201040','#302070','#104090','#0060b0','#00a0a0','#80c080','#d0c010','#f0c010','#d07010','#b02030','#701020'))
  
  xlim=NULL
  ylim=NULL
  xlab=''
  ylab=''
  
  legend=c('','','','')
  col=custom.color(4)
  
  plot(0,0,xlim=xlim,ylim=ylim,axes=FALSE,xlab=xlab,ylab=ylab,xaxs='i',yaxs='i',type='n')
  axis(1,lwd=0,lwd.ticks=0.5,tck=-0.03,cex.axis=0.8)
  axis(2,lwd=0,lwd.ticks=0.5,tck=-0.03,cex.axis=0.8)
  grid(lty=1,lwd=0.5)
  
  
  
  
  
  
  
  
  box(lwd=0.5)
  legend("topright",legend=legend,col=col,lwd=1,pch=16,bty='n',ncol=2,cex=0.8)
  dev.off()
}
