source("rscript/general/boxplot.R")
source("rscript/general/colors.R")
source("rscript/functions.R")

library(dplyr)
library(tidyr)

# Load correlation coefficients table and split pairs
pc = read.table("cortable/Table.5c.procap.corr.5M.txt",stringsAsFactors=F,header=1)
tc = pc %>%
  separate_rows(pos2, cor, sep = ",")

tc = tc %>%
  mutate(pos2 = as.numeric(pos2)) %>%
  mutate(dist = pos2 - pos1) %>%
  mutate(cor = as.numeric(cor)) %>%
  mutate(pval = cor2pvalue(cor, 78))

tc = as.data.frame(tc)

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
	chr <<- substring(p[,1],4)
	chr[chr=="X"] = "23"
	chr[chr=="Y"] = "24"
	chr = as.numeric(chr)
	chr[is.na(chr)] = 0
	return(unlist(chr*1000000000 + p[,2]))
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
	
# Count intersect by ChIP factors per each paired tTRE positions
is = lapply(ct[c("CTCF", "EP300", "RAD21", "SMC3")], function(x) count.intersect(tc[,c(1,2)], tc[,c(1,3)], x))
names(is) = c("CTCF", "P300", "RAD21", "SMC3")
#names(is) = ChIPnames

# Count overlap by ChIP factors
ol = lapply(ct[c("CTCF", "EP300", "RAD21", "SMC3")], function(x) count.overlap(tc[,c(1,2)], tc[,c(1,3)], x))
names(ol) = c("CTCF", "P300", "RAD21", "SMC3")
#names(is) = ChIPnames

plot.scatter.dist2 = function(x1, y1, x2=NULL, y2=NULL, file="pdf/test.pdf", xlim=NULL,
                              ylim=c(-1,1), xlab="Distance (kb)",
                              ylab = "Correlation coefficients",
                              cols = "#ffffff", grid = F, alpha = 1) {
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
		col = rep(cols,n)
		labels = c(rbind(rep("",n),labels))
		xaxt=F
		at = c(rbind(2*(1:n)-0.9,2*(1:n)-0.1))
		boxwd=0.7
	}
	custom.boxplot(ndata, ylim=ylim, ylab=ylab, xlab=xlab, filename=file, labels=labels, width=w, col=col, alpha = alpha,
	               grid=grid, lwd=1,close=F, xaxt=xaxt, at=at, boxwd=boxwd, outlier = T, paired = T,border.width = 1.25)
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

#######
# Finished loading and defining functions.

# extract essential data from the co-variation matrix
x = tc$dist
y = tc$cor
p = tc$pval

xlim = c(0,1000*2^(1:9)) 
col.br = c(color$blue(8)[3],color$red(8)[3])
col.gw = c("#ffffff","#e0e0e0")

######
# Functions to plot box plots of TF differences

plot.io.box = function(fac="CTCF", file="pdf/test.pdf",
                       low.tf = 1, high.tf = 2,
                       min.dist = 10000, max.dist = 500000,
                       col = col, alpha = 1) {
	x.c1 = x[is[[fac]]<=low.tf & x<max.dist & x>min.dist]
	y.c1 = y[is[[fac]]<=low.tf & x<max.dist & x>min.dist]
	x.c2 = x[is[[fac]]>=high.tf & x<max.dist & x>min.dist]
	y.c2 = y[is[[fac]]>=high.tf & x<max.dist & x>min.dist]
	id = match.list(list(x.c1,x.c2))
	plot.pv = function(y1, y2) {
	  pv = wilcox.test(y1, y2)$p.value
	  a = ""
	  if(pv<0.05) a="*"
	  if(pv<0.01) a="**"
	  if(pv<0.001) a="***"
	  text(1.5,1.01,a,cex=0.6)		  
	}
	
	custom.boxplot(list(y.c1,y.c2),file=file,ylim=c(-0.7,1),col=col, paired = T, alpha = alpha,
	               labels=list(bquote(.(fac)<=.(low.tf)), bquote(.(fac)>=.(high.tf))),
	               ylab="Correlation coefficients",lwd=1,close=F, grid = F,
	               border.width = 1.25)
	plot.pv(y.c1, y.c2)
	custom.boxplot(list(x.c1/1000,x.c2/1000),col=col, paired = T, alpha = alpha,
	               labels=list(bquote(.(fac)<=.(low.tf)), bquote(.(fac)>=.(high.tf))),
	               ylab="Distance (kb)",lwd=1,close=F, grid = F,
	               border.width = 1.25)
	custom.boxplot(list(y.c1[id[[1]]],y.c2[id[[2]]]),ylim=c(-0.7,1),col=col, paired = T, alpha = alpha,
	               labels=list(bquote(.(fac)<=.(low.tf)), bquote(.(fac)>=.(high.tf))),
	               ylab="Correlation coefficients (SRS)",lwd=1,close=F, grid = F,
	               border.width = 1.25)
	plot.pv(y.c1[id[[1]]],y.c2[id[[2]]])
	custom.boxplot(list(x.c1[id[[1]]]/1000,x.c2[id[[2]]]/1000),col=col, paired = T, alpha = alpha,
	               labels=list(bquote(.(fac)<=.(low.tf)), bquote(.(fac)>=.(high.tf))),
	               ylab="Matched distance (kb)",lwd=1,close=F, grid = F,
	               border.width = 1.25)
	x.c1 = x[ol[[fac]]==0 & x<max.dist & x>min.dist]
	y.c1 = y[ol[[fac]]==0 & x<max.dist & x>min.dist]
	x.c2 = x[ol[[fac]]>=2 & x<max.dist & x>min.dist]
	y.c2 = y[ol[[fac]]>=2 & x<max.dist & x>min.dist]
	id = match.list(list(x.c1,x.c2))
	custom.boxplot(list(y.c1,y.c2),ylim=c(-0.7,1),col=col, paired = T, alpha = alpha,
	               labels=list(bquote(.(fac)~"(-)"), bquote(.(fac)~"(+)")),
	               ylab="Correlation coefficients",lwd=1,close=F, grid = F,
	               border.width = 1.25)
	plot.pv(y.c1, y.c2)
	custom.boxplot(list(x.c1/1000,x.c2/1000),col=col, paired = T, alpha = alpha,
	               labels=list(bquote(.(fac)~"(-)"), bquote(.(fac)~"(+)")),
	               ylab="Distance (kb)",lwd=1,close=F, grid = F,
	               border.width = 1.25)
	custom.boxplot(list(y.c1[id[[1]]],y.c2[id[[2]]]),ylim=c(-0.7,1),col=col, paired = T, alpha = alpha,
	               labels=list(bquote(.(fac)~"(-)"), bquote(.(fac)~"(+)")),
	               ylab="Correlation coefficients (SRS)",lwd=1,close=F, grid = F,
	               border.width = 1.25)
	plot.pv(y.c1[id[[1]]],y.c2[id[[2]]])
	custom.boxplot(list(x.c1[id[[1]]]/1000,x.c2[id[[2]]]/1000),col=col, paired = T, alpha = alpha,
	               labels=list(bquote(.(fac)~"(-)"), bquote(.(fac)~"(+)")),
	               ylab="Matched distance (kb)",lwd=1, grid = F,
	               border.width = 1.25)
}

plot.sc.box = function(fac="CTCF", file="pdf/test.pdf",
                       low.tf = 1, high.tf = 2, col = col, alpha = 0.5) {
  n=plot.scatter.dist2(x[is[[fac]]<=low.tf],y[is[[fac]]<=low.tf],x[is[[fac]]>=high.tf],y[is[[fac]]>=high.tf],file=file,ylim=c(-0.7,1),
                       xlim=xlim, xlab="", col = col, alpha = alpha)
  legend("top", inset=c(0,-0.2), legend=c(as.expression(bquote(.(fac) <= .(low.tf))),as.expression(bquote(.(fac) >= .(high.tf)))),
         pch=22, pt.cex=1.5,lty=0, pt.bg=alpha(c("#ffffff", colorRampPalette(col)(1)), alpha),xpd=T,cex=0.8, ncol=2, bty='n')
  print(n)
  n=plot.scatter.dist2(x[ol[[fac]]==0],y[ol[[fac]]==0],x[ol[[fac]]>=2],y[ol[[fac]]>=2],file=file,ylim=c(-0.7,1),
                       xlim=xlim, xlab="", col = col, alpha = alpha)
  legend("top", inset=c(0,-0.2), legend=c(paste(fac,"(-)"),paste(fac,"(+)")),
         pch=22, pt.cex=1.5,lty=0, pt.bg=alpha(c("#ffffff", colorRampPalette(col)(1)), alpha), xpd=T,cex=0.8, ncol=2, bty='n')
  print(n)
  dev.off()
}

#####
# plot data
plot.sc.box("CTCF", "pdf/Fig3/Fig3a1.CTCF_cdistBox.pdf", col = color$greenpa(10), alpha = 0.5)
plot.io.box("CTCF", "pdf/Fig3/Fig3a2.CTCF_cSRS.pdf", col = color$greenpa(10), alpha = 0.5)
plot.sc.box("CTCF", "pdf/Fig3/Fig3a3.CTCF_distBox_0co.pdf", low.tf = 0, high.tf = 1, col = color$greenpa(10), alpha = 0.5) # only 62 pairs that has 0 intersection in 256-512 kb range
plot.io.box("CTCF", "pdf/Fig3/Fig3a4.CTCF_SRS_0co.pdf", low.tf = 0, high.tf = 1, col = color$greenpa(10), alpha = 0.5)

plot.sc.box("P300", "pdf/Fig3/Fig3b1.P300_cdistBox.pdf", col = color$orpa(10), alpha = 0.5)
plot.io.box("P300", "pdf/Fig3/Fig3b2.P300_cSRS.pdf", col = color$orpa(10), alpha = 0.5)

plot.sc.box("RAD21", "pdf/Fig3/Fig3c1.RAD21_cdistBox.pdf", col = color$pinkpa(10), alpha = 0.5)
plot.io.box("RAD21", "pdf/Fig3/Fig3c2.RAD21_cSRS.pdf", col = color$pinkpa(10), alpha = 0.5)

plot.sc.box("SMC3", "pdf/Fig3/Fig3d1.SMC3_distBox.pdf", col = color$japa(10), alpha = 0.5)
plot.io.box("SMC3", "pdf/Fig3/Fig3d2.SMC3_SRS.pdf", col = color$japa(10), alpha = 0.5)
