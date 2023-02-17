source("~/Sandbox/procap/rscript/boxplot.R")
source("~/Sandbox/procap/rscript/colors.R")

# Load correlation coefficients table and split pairs
pc = read.table("cortable/Table.5c.procap.corr.5M.txt",stringsAsFactors=F,header=1)
pos2 = lapply(strsplit(pc$pos2,","),as.numeric)
cor = lapply(strsplit(pc$cor,","),as.numeric)
pos1 = sapply(1:length(pos2), function(i) rep(pc$pos1[i],length(pos2[[i]])))
chr = sapply(1:length(pos2), function(i) rep(pc$cr = sapply(1:length(pos2), function(i) rep(pc$chr[i],length(pos2[[i]])))
r[i],length(pos2[[i]])))
tc = data.frame(chr=unlist(chr),pos1=unlist(pos1),pos2=unlist(pos2),cor=unlist(cor))
tc$dist = tc$pos2 - tc$pos1

# Load ChIP-seq dataset
ChIPnames = c("ATF2","ATF3","BATF","BCL11A","BCL3","BCLAF1","BHLHE40","BRCA1",
	"CEBPB","CHD1","CHD2","CTCF","E2F4","EBF1","EGR1","ELF1","ELK1","EP300",
	"ETS1","EZH2","FOS","FOXM1","GABPA","IKZF1","IRF4","JUND","MAX","MAZ",
	"MEF2A","MEF2C","MTA3","MXI1","MYC","NFATC1","NFE2","NFIC","NFYA","NFYB",
	"NR2C2","NRF1","PAX5","PBX3","PML","POLR2A","POLR3G","POU2F2","RAD21",
	"RCOR1","RELA","REST","RFX5","RUNX3","RXRA","SIN3A","SIX5","SMC3","SP1",
	"SPI1","SRF","STAT1","STAT3","STAT5A","TAF1","TBL1XR1","TBP","TCF12","TCF3",
	"USF1","USF2","WRNIP1","YY1","ZBTB33","ZEB1","ZNF143","ZNF274","ZZZ3")
ct = lapply(paste("ChIP/Tfbs.",ChIPnames,".txt",sep=""), function(x) read.table(x,stringsAsFactors=F,header=F))
names(ct) = ChIPnames

# Function to convert chromosomal positions to numbers
conv.pos = function(p) {
	chr = substring(p[,1],4)
	chr[chr=="X"] = "23"
	chr[chr=="Y"] = "24"
	chr = as.numeric(chr)
	chr[is.na(chr)] = 0
	return(chr*1000000000 + p[,2])
}

# Function to count intersects by pos between pos1 and pos2
count.intersect = function(pos1, pos2, pos) {
	p1 = conv.pos(pos1)
	p2 = conv.pos(pos2)
	p = conv.pos(pos)
	p = sort(p)
	i1 = findInterval(p1, p)
	i2 = findInterval(p2, p)
	return(i2-i1)
}

# Function to count overlaps
count.overlap = function(pos1, pos2, pos, range=500) {
	p1 = conv.pos(pos1)
	p2 = conv.pos(pos2)
	p = conv.pos(pos)
	p = sort(p)
	i1 = findInterval(p1, p)
	i2 = findInterval(p2, p)
	o1 = rep(F,length(p1))
	o2 = rep(F,length(p2))
	n = length(p)
	o1[i1>0] = p1[i1>0] - p[i1[i1>0]] < range
	o1[i1<n] = o1[i1<n] | (p[i1[i1<n]+1] - p1[i1<n] < range)
	o2[i2>0] = p2[i2>0] - p[i2[i2>0]] < range
	o2[i2<n] = o2[i2<n] | (p[i2[i2<n]+1] - p2[i2<n] < range)
	return(o1+o2)
}
	
# Count intersect by ChIP factors
is = lapply(ct, function(x) count.intersect(tc[,c(1,2)], tc[,c(1,3)], x))
names(is) = ChIPnames

ol = lapply(ct, function(x) count.overlap(tc[,c(1,2)], tc[,c(1,3)], x))
names(ol) = ChIPnames

plot.scatter.dist = function(x1, y1, x2=NULL, y2=NULL, file="test.pdf", xlim=NULL, ylim=c(-1,1), xlab="Distance (kb)") {
	if(is.null(xlim)) xlim = c(0,1000*2^(0:12))
	ndata = split(y1,cut(x1/1000,xlim/1000))
	n = length(xlim)-1
	labels=paste(xlim[1:n]/1000,xlim[1:n+1]/1000,sep="-")
	col = NULL
	w=0.9
	xaxt=T
	at=NULL
	boxwd=0.75
	if(!is.null(x2)) {
		n = length(ndata)
		ndata2 = split(y2,cut(x2/1000,xlim/1000))
		ndata1 = ndata
		ndata = list()
		pv=c()
		for(i in 1:n) {
			ndata[[2*i-1]] = ndata1[[i]]
			ndata[[2*i]] = ndata2[[i]]
			pv[i] = wilcox.test(ndata1[[i]],ndata2[[i]])$p.value
		}
		w=0.7
		col = rep(c(color$blue(8)[3],color$red(8)[3]),n)
		labels = c(rbind(rep("",n),labels))
		xaxt=F
		at = c(rbind(2*(1:n)-0.9,2*(1:n)-0.1))
		boxwd=0.7
	}
	custom.boxplot(ndata, ylim=ylim, ylab="Correlation coefficients", xlab=xlab, filename=file, labels=labels, width=w, col=col, grid=T, lwd=1.5,close=F, xaxt=xaxt, at=at, boxwd=boxwd)
	if(!xaxt) {
		axis(1,at=0.5+(0:n)*2,tck=0.03,lwd=0,lwd.tick=0.5,labels=F)
		for(i in 1:n) {
			a = ""
			if(pv[i]<0.05) a="*"
			if(pv[i]<0.01) a="**"
			if(pv[i]<0.001) a="***"
			text(2*i-0.5,1.01,a,cex=0.6)	
		}
	}
	return(unlist(lapply(ndata,length)))
}

x = tc$dist
y = tc$cor
xlim = c(0,1000*2^(1:9)) 
col.br = c(color$blue(8)[3],color$red(8)[3])

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

# Plot CTCF intersect boxplots
plot.scatter.dist(x[is$CTCF<2],y[is$CTCF<2],x[is$CTCF>=2],y[is$CTCF>=2],file="FigS8D.pdf",ylim=c(-0.7,1), xlim=xlim, xlab="")
legend("top", inset=c(0,-0.2), legend=c(expression(CTCF<=1),expression(CTCF>=2)), pch=22, pt.cex=1.5, lty=0, pt.bg=col.br,xpd=T,cex=0.8, ncol=2, bty='n')
dev.off()

x.c1 = x[is$CTCF<2 & x<200000 & x>10000]
y.c1 = y[is$CTCF<2 & x<200000 & x>10000]
x.c2 = x[is$CTCF>=2 & x<200000 & x>10000]
y.c2 = y[is$CTCF>=2 & x<200000 & x>10000]
id = match.list(list(x.c1,x.c2))
custom.boxplot(list(y.c1,y.c2),file="Fig8A.pdf",ylim=c(-0.7,1),col=col.br,labels=c(expression(CTCF<=1),expression(CTCF>=2)),ylab="Correlation coefficients",lwd=1.5)
custom.boxplot(list(x.c1/1000,x.c2/1000),file="FigS8B.pdf",col=col.br,labels=c(expression(CTCF<=1),expression(CTCF>=2)),ylab="Distance (kb)",lwd=1.5)
custom.boxplot(list(y.c1[id[[1]]],y.c2[id[[2]]]),file="FigS8C1.pdf",ylim=c(-0.7,1),col=col.br,labels=c(expression(CTCF<=1),expression(CTCF>=2)),ylab="Correlation coefficients",lwd=1.5)
custom.boxplot(list(x.c1[id[[1]]]/1000,x.c2[id[[2]]]/1000),file="FigS8C2.pdf",col=col.br,labels=c(expression(CTCF<=1),expression(CTCF>=2)),ylab="Distance (kb)",lwd=1.5)


# Plot P300 overlap boxplots
plot.scatter.dist(x[ol$EP300==0],y[ol$EP300==0],x[ol$EP300>=2],y[ol$EP300>=2],file="FigS8F.pdf",ylim=c(-0.7,1), xlim=xlim, xlab="")
legend("top", inset=c(0,-0.2), legend=c("P300 (-)","P300 (+)"), pch=22, pt.cex=1.5,lty=0, pt.bg=col.br,xpd=T,cex=0.8, ncol=2, bty='n')
dev.off()
custom.boxplot(list(y[ol$EP==0&x<200000],y[ol$EP==2&x<200000]),file="Fig8B.pdf",ylim=c(-0.7,1),col=col.br,labels=c("P300 (-)","P300 (+)"),ylab="Correlation coefficients",lwd=1.5)
custom.boxplot(list(x[ol$EP==0&x<200000]/1000,x[ol$EP==2&x<200000]/1000),file="FigS8E.pdf",ylim=c(0,200),col=col.br,labels=c("P300 (-)","P300 (+)"),ylab="Distance (kb)",lwd=1.5)

plot.io.box = function(fac="CTCF",file="FigS8X.pdf",cutoff=2) {
	x.c1 = x[is[[fac]]<cutoff & x<200000 & x>10000]
	y.c1 = y[is[[fac]]<cutoff & x<200000 & x>10000]
	x.c2 = x[is[[fac]]>=cutoff & x<200000 & x>10000]
	y.c2 = y[is[[fac]]>=cutoff & x<200000 & x>10000]
	id = match.list(list(x.c1,x.c2))
	custom.boxplot(list(y.c1[id[[1]]],y.c2[id[[2]]]),file=file,ylim=c(-0.7,1),col=col.br,
		labels=c(expression(TF <= 1), expression(TF >= 2)),
		ylab="Correlation coefficients",lwd=1.5,close=F)
	x.c1 = x[ol[[fac]]==0 & x<200000 & x>10000]
	y.c1 = y[ol[[fac]]==0 & x<200000 & x>10000]
	x.c2 = x[ol[[fac]]>=2 & x<200000 & x>10000]
	y.c2 = y[ol[[fac]]>=2 & x<200000 & x>10000]
	id = match.list(list(x.c1,x.c2))
custom.boxplot(list(y.c1[id[[1]]],y.c2[id[[2]]]),ylim=c(-0.7,1),col=col.br,labels=c(paste(fac,"(-)"),paste(fac,"+")),ylab="Correlation coefficients",lwd=1.5)
}

plot.sc.box = function(fac="CTCF", file="FigS8X.pdf",cutoff=2) {
plot.scatter.dist(x[is[[fac]]<cutoff],y[is[[fac]]<cutoff],x[is[[fac]]>=cutoff],y[is[[fac]]>=cutoff],file=file,ylim=c(-0.7,1), xlim=xlim, xlab="")
legend("top", inset=c(0,-0.2), legend=c(as.expression(bquote(.(fac) <= .(cutoff-1))),as.expression(bquote(.(fac) >= .(cutoff)))), pch=22, pt.cex=1.5,lty=0, pt.bg=col.br,xpd=T,cex=0.8, ncol=2, bty='n')
plot.scatter.dist(x[ol[[fac]]==0],y[ol[[fac]]==0],x[ol[[fac]]>=2],y[ol[[fac]]>=2],file="FigS8X.pdf",ylim=c(-0.7,1), xlim=xlim, xlab="")
legend("top", inset=c(0,-0.2), legend=c(paste(fac,"(-)"),paste(fac,"(+)")), pch=22, pt.cex=1.5,lty=0, pt.bg=col.br,xpd=T,cex=0.8, ncol=2, bty='n')
dev.off()
}

plot.sc.box("CTCF", file="Fig8C1.pdf", cutoff=2)
plot.sc.box("EP300", file="Fig8C2.pdf", cutoff=2)
plot.sc.box("RAD21", file="Fig8C3.pdf", cutoff=2)
plot.sc.box("SMC3", file="FigS8F.pdf", cutoff=2)
