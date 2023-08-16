custom.bigheatmap=function(Vmat,xlab=NULL,zlim=NULL,width=NULL, height=NULL, color=NULL,b=TRUE,filename='test.pdf',scalebar=1,stretch=1)
{
	library(pixmap)
	if(is.null(xlab)) xlab=colnames(Vmat)
	if(is.null(zlim)) zlim=c(min(Vmat),max(Vmat))
	if(is.null(color))
	{
		color=c('#053061','#2166ac','#4393c3','#92c5de','#d1e5f0','#f7f7f7','#fddbc7','#f4a582','#d6604d','#b2182b','#67001f')
		if(zlim[2]< -zlim[1]) zlim[2]=-zlim[1]
		if(zlim[1]> -zlim[2]) zlim[1]=-zlim[2]
	}
	else if(color[1]=='r') color=c('#f7f7f7','#fddbc7','#f4a582','#d6604d','#b2182b','#67001f')
	else if(color[1]=='b') color=c('#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')
	else if(color[1]=='p') color=c('#352187','#0f5cdd','#1481d6','#06a4ca','#2eb7a4','#87bf77','#d1bb59','#fec832','#f9fb0e')
    else if(color[1]=='j') color=c('#201040','#302070','#104090','#0060b0','#00a0a0','#80c080','#d0c010','#f0c010','#d07010','#b02030','#701020')
    else if(color[1]=='g') color=gray((1-2^(0:24/8)/8)*8/7)
    else if(length(color)==1) color=c('white',color)

	zw=zlim[2]-zlim[1]
	z0=zlim[1]
	Vmat[Vmat<zlim[1]]=zlim[1]
	Vmat[Vmat>zlim[2]]=zlim[2]
	nR=nrow(Vmat)
	nC=ncol(Vmat)
	colorfunc=colorRamp(color)
    if(is.null(width)|is.null(height)) pdf(filename,width=1.2+0.4*(nC+2*scalebar),height=1+stretch*nC)
    else pdf(filename,width=width,height=height)
	par(mar=c(1,4,4,2),mgp=c(1.5,0.3,0))
	plot(0,0,type='n',xlim=c(0,1+2/nC*scalebar),ylim=c(zlim[1],zlim[2]),axes=FALSE,xaxs='i',yaxs='i',xlab='',ylab='')
	mat=Vmat
	nnC=nC
	if(nC<100)
	{
		nCrep=ceiling(200/nC)
		mat=matrix(rep(Vmat[,1],nCrep),nrow=nR)
		for(i in 2:nC) mat=cbind(mat,matrix(rep(Vmat[,i],nCrep),nrow=nR))
		nnC=nC*nCrep
	}
	rasterImage(array(colorfunc((mat-z0)/zw)/255,c(nR,nnC,3)),0,zlim[1],1,zlim[2],interpolate=FALSE)
	for(i in 1:100)
	{
		rect(1+1/nC*scalebar,(i-1)/100*zw+z0,1+2/nC*scalebar,i/100*zw+z0,col=rgb(colorfunc((i-1)/100),max=255),border=rgb(colorfunc((i-1)/100),max=255))
	}	
	axis(4,cex.axis=1,tck=-0.2/(nC+2),lwd=0,lwd.tick=1)
	#axis(4,at=c(zlim[1],zlim[2]),labels=FALSE,lwd.tick=0)
    rect(0,zlim[1],1,zlim[2],lwd=1,xpd=T)
    rect(1+1/nC*scalebar,zlim[1],1+2/nC*scalebar,zlim[2],lwd=1,xpd=T)
	text(labels=xlab,x=seq(0.2/nC,1-0.8/nC,by=1/nC),y=rep(1.02,nR)*zw+z0,srt=45,pos=4,xpd=TRUE,cex=1)
	dev.off()
}

custom.bmpheatmap=function(vmat,basename='test',width=2.4,height=1.5,label=c('0%','100%'),at=c(0,100),
                           lwd = 0.5, scale.barwidth=5,scale.label=c('0%','100%'),scale.at=c(0,100))
{
	library(bmp)
	library(pixmap)
	image<-pixmapRGB(read.bmp(paste(basename,'.main.bmp',sep='')))
	scalebar<-pixmapRGB(read.bmp(paste(basename,'.scalebar.bmp',sep='')))
	pdf(paste(basename,'.pdf',sep=''),width=width,height=height)
	par(mar=c(1.5,1,1,1.5),mgp=c(1.5,0.2,0))
	plot(0,0,xlim=c(0,100+scalebar*2),ylim=c(0,100),axes=FALSE,xlab='',ylab='',xaxs='i',yaxs='i')
	rasterImage(array(c(image@red,image@green,image@blue),c(image@size[1],image@size[2],3)),0,0,100,100)
	rasterImage(array(c(t(scalebar@red),t(scalebar@green),t(scalebar@blue)),c(scalebar@size[2],scalebar@size[1],3)),100+scalebar,100,100+scalebar*2,0)
	axis(1,at=at,label=label,cex.axis=0.8,tck=-0.03,lwd=0.5)
	axis(4,at=scale.at,label=scale.label,cex.axis=0.7,tck=-0.03,lwd=0.5,mgp=c(1.5,0.1,0))
	rect(0,0,100,100,lwd=lwd,xpd=T)
	rect(100+scalebar,0,100+scalebar*2,100,lwd=lwd,xpd=T)
	dev.off()
}

custom.heatmap=function(Vmat,xlab=NULL,ylab=NULL,zlim=NULL,width=NULL,height=NULL,color=NULL,b=TRUE,filename='test.pdf',scalebar=1,stretch=1)
{
	dense=FALSE
	if(is.null(xlab)) xlab=colnames(Vmat)
	if(is.null(ylab)) dense=TRUE
	if(is.null(zlim)) zlim=c(min(Vmat),max(Vmat))
	if(is.null(color))
	{
		color=c('#053061','#2166ac','#4393c3','#92c5de','#d1e5f0','#f7f7f7','#fddbc7','#f4a582','#d6604d','#b2182b','#67001f')
		if(zlim[2]< -zlim[1]) zlim[2]=-zlim[1]
		if(zlim[1]> -zlim[2]) zlim[1]=-zlim[2]
	}
    else if(color[1]=='r') color=c('#f7f7f7','#fddbc7','#f4a582','#d6604d','#b2182b','#67001f')
    else if(color[1]=='b') color=c('#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')
    else if(color[1]=='p') color=c('#352187','#0f5cdd','#1481d6','#06a4ca','#2eb7a4','#87bf77','#d1bb59','#fec832','#f9fb0e')
    else if(color[1]=='j') color=c('#201040','#302070','#104090','#0060b0','#00a0a0','#80c080','#d0c010','#f0c010','#d07010','#b02030','#701020')
    else if(length(color)==1) color=c('white',color)

    else if(color=='g') color=gray((1-2^(0:24/8)/8)*8/7)
	zw=zlim[2]-zlim[1]
	z0=zlim[1]
	Vmat[Vmat>zlim[2]]=zlim[2]
	Vmat[Vmat<zlim[1]]=zlim[1]

	nR=nrow(Vmat)
	nC=ncol(Vmat)
	colorfunc=colorRamp(color)
    if(is.null(width)|is.null(height)) 
    {
    if(dense)
		pdf(filename,width=1.2+0.4*(nC+2*scalebar),height=1+stretch*nC)
	else
		pdf(filename,width=1.2+0.4*(nC+2*scalebar),height=1+0.4*nR)
	}
    else pdf(filename,width=width,height=height)

	par(mar=c(1,4,4,2),mgp=c(1.5,0.3,0))

	plot(0,0,type='n',xlim=c(0,1+2/nC*scalebar),ylim=c(zlim[1],zlim[2]),axes=FALSE,xaxs='i',yaxs='i',xlab='',ylab='')

	for(i in 1:nR)
	{
		for(j in 1:nC)
		{
            if(is.na(Vmat[i,j])) rect((j-1)/nC,(nR-i)/nR*zw+z0,j/nC,(nR-i+1)/nR*zw+z0,col='#ffffff',border=NA)
			else
            {
                nV=(Vmat[i,j]-z0)/zw
                if(nV>1) nV=1
                if(nV<0) nV=0
                if(is.na(b)) rect((j-1)/nC,(nR-i)/nR*zw+z0,j/nC,(nR-i+1)/nR*zw+z0,col=rgb(colorfunc(nV),max=255),border=NA)
                else rect((j-1)/nC,(nR-i)/nR*zw+z0,j/nC,(nR-i+1)/nR*zw+z0,col=rgb(colorfunc(nV),max=255),border=rgb(colorfunc(nV),max=255))
            }
		}
	}
	for(i in 1:100)
	{
		rect(1+1/nC*scalebar,(i-1)/100*zw+z0,1+2/nC*scalebar,i/100*zw+z0,col=rgb(colorfunc((i-1)/100),max=255),border=rgb(colorfunc((i-1)/100),max=255))
	}	
	axis(4,cex.axis=1,tck=-0.2/(nC+2),lwd=0,lwd.tick=1)
	#axis(4,at=c(zlim[1],zlim[2]),labels=FALSE,lwd.tick=0)
    rect(0,zlim[1],1,zlim[2],lwd=1,xpd=T)
    rect(1+1/nC*scalebar,zlim[1],1+2/nC*scalebar,zlim[2],lwd=1,xpd=T)
	if(!dense) text(labels=ylab,y=seq(1-0.5/nR,0.5/nR,by=-1/nR)*zw+z0,x=rep(0,nC),pos=2,xpd=TRUE,cex=1)
	text(labels=xlab,x=seq(0.2/nC,1-0.8/nC,by=1/nC),y=rep(1.02,nR)*zw+z0,srt=45,pos=4,xpd=TRUE,cex=1)
	dev.off()
}

heatmap.1col = function(mat, col='', xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), lwd=1) {
	if(col=='r') cf = colorRamp(c('#f7f7f7','#fddbc7','#f4a582','#d6604d','#b2182b','#67001f'))
	else if(col=='y') cf = colorRamp(c('#f7f7f7','#fdf7a0','#ffc16a','#f0b020','#d0a020','#907010'))
	else if(col=='g') cf = colorRamp(gray((1-2^(0:24/8)/8)*8/7))
	else if(col=='j') cf = colorRamp(c('#201040','#302070','#104090','#0060b0','#00a0a0','#80c080','#d0c010','#f0c010','#d07010','#b02030','#701020'))
	else if(col=='b') cf = colorRamp(c('#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061'))
	else cf = colorRamp(c('#053061','#2166ac','#4393c3','#92c5de','#d1e5f0','#f7f7f7','#fddbc7','#f4a582','#d6604d','#b2182b','#67001f'))
	nC = ncol(mat)
	nR = nrow(mat)
	mat[mat<zlim[1]] = zlim[1]
	mat[mat>zlim[2]] = zlim[2]
	mat = (mat-zlim[1])/(zlim[2]-zlim[1])
	rasterImage(array(cf(mat)/255,c(nR,nC,3)),xlim[1],ylim[1],xlim[2],ylim[2],interpolate=F,xpd=T)
	if(lwd) rect(xlim[1],ylim[1],xlim[2],ylim[2],lwd= lwd,xpd=T)
}

heatmap.2col = function(mat1, mat2, xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), lwd=1) {
	cf1 = colorRamp(c('#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061'))
	cf2 = colorRamp(c('#f7f7f7','#fddbc7','#f4a582','#d6604d','#b2182b','#67001f'))
	nC = ncol(mat1)
	nR = nrow(mat2)
	mat1[mat1<zlim[1]] = zlim[1]
	mat1[mat1>zlim[2]] = zlim[2]
	mat2[mat2<zlim[1]] = zlim[1]
	mat2[mat2>zlim[2]] = zlim[2]
	mat1 = (mat1-zlim[1])/(zlim[2]-zlim[1])
	mat2 = (mat2-zlim[1])/(zlim[2]-zlim[1])
#	cmap = 255 - ((255-cf1(mat1)) + (255-cf2(mat2)))
	cmap = 255 - sqrt((255-cf1(mat1))^2 + (255-cf2(mat2))^2)
	cmap[cmap>255] = 255
	cmap[cmap<0] = 0
#	cmap = 255 - sqrt(((255-cf1(mat1))^2+(255-cf2(mat2))^2)/2)
#	cmap = sqrt( (cf1(mat1)^2 + cf2(mat2)^2)/2 )
	rasterImage(array(cmap/255,c(nR,nC,3)),xlim[1],ylim[1],xlim[2],ylim[2],interpolate=F,xpd=T)
	if(lwd) rect(xlim[1],ylim[1],xlim[2],ylim[2],lwd=lwd,xpd=T)
}
