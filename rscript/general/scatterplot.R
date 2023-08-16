custom.scatter = function(x, y, filename='test.pdf', xlim=NULL, ylim=NULL, main=NULL, xlab=NULL, ylab=NULL, scalelab=NULL, width=3, height=3, res=300, rad=8, scalewidth=0.2, diag=FALSE, h=NULL, v=NULL, color=NULL)
{
    custom.scatter.open(filename,width,height,res,main,xlab,ylab,scalelab,scalewidth)
    custom.scatter.plot(x,y,xlim,ylim,width,height,res,rad,scalewidth,diag,h,v,color)
    custom.scatter.label(main,xlab,ylab,scalelab,scalewidth)
    custom.scatter.close()
}

custom.scatter.open = function(filename='test.pdf', width=3, height=3, res=300, main=NULL, xlab=NULL, ylab=NULL, scalelab=NULL, scalewidth=0.2)
{
    is.main = 0+!is.null(main)
    is.xlab = 0+!is.null(xlab)
    is.ylab = 0+!is.null(ylab)
    is.scalelab = 0+!is.null(scalelab)

    pdf(filename, width=width+scalewidth*2+0.2*(4+is.ylab+is.scalelab), height=height+0.2*(3+is.xlab+2*is.main))
    par(mar=c(2+is.xlab,2+is.ylab,1+2*is.main,2+2*scalewidth/0.2+is.scalelab),mgp=c(1.5,0.3,0))
}

custom.scatter.plot = function(x, y, xlim=NULL, ylim=NULL, width=3, height=3, res=300, rad=8,
                               scalewidth=0.2, diag=FALSE, h=NULL, v=NULL, color=NULL, add=FALSE,
                               noimage = FALSE, pval=T, lwd=0.5, cex=0.8, alpha = 1, border.width = 1)
{
    ix=x[is.finite(x)&is.finite(y)]
    iy=y[is.finite(x)&is.finite(y)]
    minx=floor(min(ix))
    maxx=floor(max(ix)+1)
    miny=floor(min(iy))
    maxy=floor(max(iy)+1)
    if(!is.null(xlim)) {
        minx=xlim[1]
        maxx=xlim[2]
    }
    if(!is.null(ylim)) {
        miny=ylim[1]
        maxy=ylim[2]
    }
    if(add) {
        minx=par('xaxp')[1]
        maxx=par('xaxp')[2]
        miny=par('yaxp')[1]
        maxy=par('yaxp')[2]
		width=par()$pin[1]
		height=par()$pin[2]
    }
    xres=width*res
    yres=height*res

    rangex=maxx-minx
    rangey=maxy-miny

	minx=minx-0.03*rangex
	maxx=maxx+0.03*rangex
	miny=miny-0.03*rangey
	maxy=maxy+0.03*rangey
	rangex=maxx-minx
	rangey=maxy-miny

    bx=(ix-minx)/rangex
    by=(iy-miny)/rangey
    cx=bx[bx>0.02&bx<0.98&by>0.02&by<0.98]
    cy=by[bx>0.02&bx<0.98&by>0.02&by<0.98]
    mat=matrix(rep(0,xres*yres),ncol=xres,nrow=yres)
    dot=matrix(rep(0,(2*rad+1)^2),ncol=2*rad+1,nrow=2*rad+1)
    for(i in 1:(2*rad+1)) for(j in 1:(2*rad+1))
        if((i-rad-0.5)^2+(j-rad-0.5)^2<rad^2) dot[i,j]=1
    for(i in 1:length(cx)) {
        px=cx[i]*xres
        py=cy[i]*yres
        mpy=round(py-rad)
        Mpy=mpy+2*rad
        mpx=round(px-rad)
        Mpx=mpx+2*rad
        
        if(mpy>0&Mpy<=yres&mpx>0&Mpx<=xres)  mat[mpy:Mpy,mpx:Mpx]=mat[mpy:Mpy,mpx:Mpx]+dot
    }

    mat[mat==0]=-1
    mat[mat>0]=log10(mat[mat>0])
    maxz=max(mat)

    if(is.null(color)) color=c('#352187','#0f5cdd','#1481d6','#06a4ca','#2eb7a4','#87bf77','#d1bb59','#fec832','#f9fb0e')
    else if(color[1]=='jet') color=c('#201040','#302070','#104090','#0060b0','#00a0a0','#80c080','#d0c010','#f0c010','#d07010','#b02030','#701020')
    else if(color[1]=='grayscale')
        color=gray(2^(1:24/8)/8*0.9)
    else if(length(color)==1) color=c(color,color)       
    
    colorLength=length(color)
    color=c(color,'#ffffff')
    color = alpha(color, alpha)
    colorfunc=colorRamp(color,alpha=T)

    if(maxz==0) maxz=1 
    mat=mat/maxz*(colorLength-1)/colorLength
    mat[mat<0]=1

    map=colorfunc(mat)/255

    library(pixmap)

    if(!add) plot(0,0,xlim=c(minx,maxx),ylim=c(miny,maxy),axes=FALSE,xlab='',ylab='',xaxs='i',yaxs='i',type='n')

    if(!add) {
        if(diag) segments(minx,miny,maxx,maxy,col='#a0a0a0')
        if(!is.null(h)) abline(h=h,col='#a0a0a0')
        if(!is.null(v)) abline(v=v,col='#a0a0a0')
    }
    if(!noimage) rasterImage(array(map,c(yres,xres,4)),minx,maxy,maxx,miny)
    
    if(!add) {
        if(!noimage) rasterImage(array(colorfunc((yres:1)/(yres-1)*(colorLength-1)/colorLength)/255,c(yres,1,4)),maxx+rangex*scalewidth/width,miny,maxx+rangex*2*scalewidth/width,maxy,xpd=TRUE)
        axis(1,lwd=0,lwd.ticks=lwd,tck=-0.03,cex.axis=cex)
        axis(2,lwd=0,lwd.ticks=lwd,tck=-0.03,cex.axis=cex)

        par(xpd=TRUE)
        rect(minx,miny,maxx,maxy,lwd=border.width)
        rect(minx+rangex*(1+scalewidth/width),miny,minx+rangex*(1+2*scalewidth/width),maxy,lwd=lwd)
        if(pval) text(minx+rangex*0.05,y=miny+rangey*0.95,labels=c(paste('r=',format(round(cor(ix,iy),4),nsmall=4))),adj=0, cex=cex)
        par(new=TRUE)

        plot(0,0,type='n',xlim=c(minx,minx+rangex),ylim=c(0,maxz),axes=FALSE,xlab='',ylab='',xaxs='i',yaxs='i')
        axis(4,pos=minx+rangex*(1+2*scalewidth/width),lwd=0,lwd.ticks=lwd,tck=-0.03,cex.axis=cex)
        par(new=TRUE)
        plot(0,0,xlim=c(minx,maxx),ylim=c(miny,maxy),axes=FALSE,xlab='',ylab='',xaxs='i',yaxs='i',type='n')
        par(xpd=FALSE)
    }
}

custom.scatter.label = function(main=NULL,xlab=NULL,ylab=NULL,scalelab=NULL,scalewidth=0.2, cex=1)
{
    mtext(main,3,1,font=2,cex=1.2*cex)
    mtext(xlab,1,1.6,cex=cex)
    mtext(ylab,2,1.6,cex=cex)
    mtext(scalelab,4,1.4+2*scalewidth/0.2,cex=0.9*cex)
}

custom.scatter.close = function()
{
    par(xpd=FALSE)
    par(new=FALSE)
    dev.off()
}
