interaction.map=function(table,chr,start,end=NULL,file='test.pdf')
{
	if(is.null(end)) {end=start+500000; start=start-500000}
	table.chr=table[,1]
	table.pos=table[,2]
	r= table.chr==chr & table.pos>start & table.pos<end
	n=sum(r)
	print(n)
	if(n<2) return()
	data=table[r,-(1:2)]
	pos=table.pos[r]
	cv=cor(t(data))

	#color=colorRamp(c('#201040','#302070','#104090','#0060b0','#00a0a0','#80c080','#d0c010','#f0c010','#d07010','#b02030','#701020')
	color=colorRamp(c('#053061','#2166ac','#4393c3','#92c5de','#d1e5f0','#f7f7f7','#fddbc7','#f4a582','#d6604d','#b2182b','#67001f'))
	
	pdf(file,width=4,height=3)
	par(mar=c(3,3,1,1),mgp=c(1.5,0.3,0))
	xlim=c(start,end)
	ylim=c(0,end-start)
	xlab='Position'
	ylab=chr
	plot(0,0,xlim=xlim,ylim=ylim,axes=FALSE,xlab=xlab,ylab=ylab,xaxs='i',type='n')
	x=c()
	y=c()
	ccv=c()
	for(i in 1:(n-1))
	{
		for(j in (i+1):n)
		{
			x=c(x,(start+end)/2+(pos[i]+pos[j]-start-end)/2)
			y=c(y,pos[j]-pos[i])
			ccv=c(ccv,cv[i,j])
		}
	}
	o=order(abs(ccv))
	x=x[o]
	y=y[o]
	col=rgb(color(ccv[o]/2+0.5)/255)
	points(x,y,pch=16,cex=0.4,col=col,bg=col,lwd=0)

    axis(1,lwd=0,lwd.tick=0.5,tck=-0.015,cex.axis=0.6,at=axTicks(1),labels=formatC(axTicks(1),format='d',big.mark=','))
	box(lwd=.5)
	dev.off()
}

interaction.bitmap=function(table,chr,start,end=NULL,file='test.pdf',xres=1600,yres=1000,rad=4,cvmax=0.4,weight=T)
{
	if(is.null(end)) {end=start+500000; start=start-500000}
	table.chr=table[,1]
	table.pos=table[,2]
	r= table.chr==chr & table.pos>start & table.pos<end
	n=sum(r)
	print(n)
	if(n<2) return()
	data=table[r,-(1:2)]
	pos=table.pos[r]
	cv=cor(t(data))
	m=apply(data,1,mean)
	if(weight==F) m=rep(1,n)
	#color=colorRamp(c('#201040','#302070','#104090','#0060b0','#00a0a0','#80c080','#d0c010','#f0c010','#d07010','#b02030','#701020')
	color=colorRamp(c('#053061','#2166ac','#4393c3','#92c5de','#d1e5f0','#ffffff','#fddbc7','#f4a582','#d6604d','#b2182b','#67001f'))
	
	pdf(file,width=4,height=3)
	par(mar=c(3,3,1,1),mgp=c(1.5,0.3,0))
	xlim=c(start,end)
	ylim=c(0,end-start)
	xlab='Position'
	ylab=chr
	x=c()
	y=c()
	ccv=c()
	mm=c()
	for(i in 1:(n-1))
	{
		for(j in (i+1):n)
		{
			x=c(x,(start+end)/2+(pos[i]+pos[j]-start-end)/2)
			y=c(y,pos[j]-pos[i])
			ccv=c(ccv,cv[i,j])
			mm=c(mm,sqrt(m[i]*m[j]))
		}
	}
	mat=matrix(rep(0,xres*yres),ncol=xres,nrow=yres)
	count=matrix(rep(0,xres*yres),ncol=xres,nrow=yres)
	dot=matrix(rep(0,(2*rad+1)^2),ncol=2*rad+1,nrow=2*rad+1)
	circ=matrix(rep(0,(2*rad+1)^2),ncol=2*rad+1,nrow=2*rad+1)
	for(i in 1:(2*rad+1)) for(j in 1:(2*rad+1))
	{
		if((i-rad-0.5)^2+(j-rad-0.5)^2<rad^2) dot[i,j]=1
		if((i-rad-0.5)^2+(j-rad-0.5)^2<rad^2) circ[i,j]=1
	}
	for(i in 1:(2*rad+1)) for(j in 1:(2*rad+1))
		if((i-rad-0.5)^2+(j-rad-0.5)^2<(rad*0.8)^2) circ[i,j]=0
	

	vx=(x-start)/(end-start)*(xres-2*rad-1)+rad+1
	vy=y/(end-start)*(yres-2*rad-1)+rad+1
	
	for(i in 1:length(x))
	{
		px=vx[i]
		py=vy[i]
		mat[(py-rad):(py+rad),(px-rad):(px+rad)]=mat[(py-rad):(py+rad),(px-rad):(px+rad)]+dot*ccv[i]*mm[i]
		count[(py-rad):(py+rad),(px-rad):(px+rad)]=count[(py-rad):(py+rad),(px-rad):(px+rad)]+dot*mm[i]
	}

	dotmat=matrix(rep(0,xres*yres),ncol=xres,nrow=yres)
	for(i in 1:n)
	{
		px=(pos[i]-start)/(end-start)*(xres-2*rad-1)+rad+1
		py=rad+1
		dotmat[(py-rad):(py+rad),(px-rad):(px+rad)]=dotmat[(py-rad):(py+rad),(px-rad):(px+rad)]+circ
	}

	mat[count>0]=mat[count>0]/count[count>0]
	mat=(mat/cvmax+1)/2
	mat[mat>1]=1
	mat[mat<0]=0
	map=color(mat)/255
	map[dotmat>0,]=c(0.2,0.2,0.2)
	library(pixmap)

	plot(0,0,xlim=xlim,ylim=ylim,axes=FALSE,xlab=xlab,ylab=ylab,xaxs='i',type='n')
    rasterImage(array(map,c(yres,xres,3)),start-(rad+1)/xres*(end-start),(1+(rad+1)/yres)*(end-start),end+(rad+1)/xres*(end-start),-(rad+1)/yres*(end-start))
	    
    axis(1,lwd=0,lwd.tick=0.5,tck=-0.015,cex.axis=0.7,at=axTicks(1),labels=formatC(axTicks(1),format='d',big.mark=','))
	box(lwd=0.5)
	dev.off()
}
snp.map=function(chr,start,end=NULL,file='test.pdf')
{
	if(is.null(end)) {end=start+500000; start=start-500000}
	system(paste("tabix snp.txt.gz ",chr,":",start,"-",end," > snp.tmp",sep=""))
	snptable=read.table("snp.tmp")
	nr=nrow(snptable)
	if(nr>400) snptable=snptable[sort(sample(1:nr,400)),]
	interaction.bitmap(snptable,chr,start,end,weight=F,rad=3)
}
	
