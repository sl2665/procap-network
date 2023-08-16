source('rscript/functions.R')
source('rscript/colors.R')

pro = read.table('readcount/procap.var.norm.txt',stringsAsFactors=F)
pro.near = read.table('window/procap.window.all-all.5M.txt',stringsAsFactors=F)

pro.window = make.window.list(pro.near)
pro.cor  = cor.pro(pro.window, pro, pro)
pro.cor.cdc = quantile.curve(pro.cor$x,pro.cor$y,n=10000,p.val=c(0.05,0.5,0.95),spar=0)
plot.scatter.dist(pro.cor$x,pro.cor$y,cdc=pro.cor.cdc$r, xlim=c(0,5000000),file="FigS7A.pdf")

# Interchromosomal correlations
ich.matrix = matrix(ceiling(runif(2000000)*nrow(pro)),ncol=2)
ich.matrix[pro[ich.matrix[,1],1]==pro[ich.matrix[,2],1],]=c(0,0)
ich.matrix = ich.matrix[ich.matrix[,1]>0,]

ich.cor = apply(ich.matrix,1,function(x) cor(unlist(pro[x[1],-(1:2)]),unlist(pro[x[2],-(1:2)])))

# Q-Q plot for normality in 10,000 samples from the 5Mb distance correlations
pdf(file="pdf/Fig1S/FigS1a.qq_5Mover.pdf",width=3,height=3.4)
par(mar=c(3,3,3,1),mgp=c(1.5,0.4,0))
y=sample(pro.cor$y[pro.cor$x>1000000],10000)
interchromosomal = y
qqnorm(y,axes=F,pch=16,cex=0.5,bty='n',
       main = "Long range (> 1 Mb)",
       xlab = "Predicted quantiles",
       ylab = "Sample quantiles")
qqline(y, col="#C0C0C0", lwd=0.5)
axis(1,lwd=0,lwd.ticks=0.5, tck=-0.03, cex.axis=0.75)
axis(2,lwd=0,lwd.ticks=0.5, tck=-0.03, cex.axis=0.75)
box(lwd=0.5)
dev.off()

# Q-Q plot for normality in 10,000 samples from interchromosomal correlations
pdf(file="pdf/Fig1S/FigS1b.qq_interchromosomal.pdf",width=3,height=3.4)
par(mar=c(3,3,3,1),mgp=c(1.5,0.4,0))
y=sample(ich.cor,10000)
qqnorm(y,axes=F,pch=16,cex=0.5,bty='n',
       main = "Interchromosomal",
       xlab = "Predicted quantiles",
       ylab = "Sample quantiles")
qqline(y, col="#C0C0C0", lwd=0.5)
axis(1,lwd=0,lwd.ticks=0.5, tck=-0.03, cex.axis=0.75)
axis(2,lwd=0,lwd.ticks=0.5, tck=-0.03, cex.axis=0.75)
box(lwd=0.5)
dev.off()

my = mean(ich.cor)
sy = sd(ich.cor)

# Histogram of the correlation coefficients
cor.h=hist(ich.cor,plot=F,breaks=seq(-1,1,by=0.02))
pdf("pdf/Fig1S/FigS1c.hist_correlation.pdf",width=3.6,height=2.4)
par(mar=c(3,3,1,1),mgp=c(1.5,0.3,0))
par(lwd=0.5)
plot(cor.h,freq=F,main="",xlab="Correlation coefficients",ylab="Density",ylim=c(0,3),col="#C0C0C0",xaxs='i',yaxs='i',axes=F,lwd=0.5)
lines(-100:100/100,dnorm(-100:100/100,my,sy),lty=2,lwd=1)
axis(1,lwd=0,tck=-0.03,lwd.ticks=0.5,cex.axis=0.8)
axis(2,lwd=0,tck=-0.03,lwd.ticks=0.5,cex.axis=0.8)
box(lwd=0.5)
legend("topright",c("Interchromosomal","Gaussian fit"),pt.bg="#C0C0C0",pch=c(22,32),bty='n',lty=c(0,2),cex=0.8,lwd=c(0.5,1))
dev.off()

# Calculate p-values from the normal distribution
my = mean(ich.cor)
sy = sd(ich.cor)
pro.cor.p = lapply(1:length(pro.cor$scd$corr),function(i) 1-pnorm(pro.cor$scd$corr[[i]],my,sy))
pro.cor.fdr = p.adjust(unlist(pro.cor.p),method="fdr")
cor.cutoff = min(pro.cor$y[pro.cor.fdr<0.1])

plot.scatter.dist(pro.cor$x,pro.cor$y,cdc=pro.cor.cdc$r, xlim=c(0,5000000),file="FigS7E.pdf", cutoff=cor.cutoff)

# Save tables
pro.cor.table = t(sapply(1:length(pro.cor$scd$corr), function(i) {
	unlist(c(pro[i,1:2],paste(pro[pro.window[[i]],2],collapse=","),paste(pro.cor$scd$corr[[i]],collapse=","))) }))
colnames(pro.cor.table) = c("chr","pos1","pos2","cor")
write.table(pro.cor.table, "../table/Table.5c.procap.corr.5M.txt",quote=F,sep="\t",col.names=T,row.names=F)

pro.cor.sig = lapply(1:length(pro.cor$scd$corr), function(i) {co = pro.cor$scd$corr[[i]]>=cor.cutoff;
	if(sum(co)>0) return(rbind(i,pro.window[[i]][co]))})
pro.cor.sig = t(matrix(unlist(pro.cor.sig),nrow=2))
pro.cor.pos = data.frame(chr=pro[pro.cor.sig[,1],1], pos1=pro[pro.cor.sig[,1],2], pos2=pro[pro.cor.sig[,2],2])
write.table(pro.cor.pos, "../table/Table.5d.procap.sigcorpos.5M.txt",quote=F, sep="\t", col.names=T, row.names=F)

# PRO-cap read count
View(pro.enh)
View(enh.pos.all)
caplist = read.csv("readcount/caplist.csv")
samplist = read.table("readcount/procap.ambr.list.txt", col.names = "id")
rawread = read.table("readcount/Table.3a.procap.readcount.ambr.txt", col.names = c("chr", "pos", samplist$id))
rawread.enhsum = rawread %>%
  inner_join(enh.pos.all %>% select("chr", "pos") %>% distinct,
             by = c("chr", "pos")) %>%
  select(-chr, -pos) %>%
  colSums()
rawread.prmsum = rawread %>%
  anti_join(enh.pos.all %>% select("chr", "pos") %>% distinct,
            by = c("chr", "pos")) %>%
  select(-chr, -pos) %>%
  colSums()
totalread = samplist %>%
  left_join(caplist %>% select("id", "mappedReads"), by = "id")

RPM_on_enhancers = data.frame(sample = totalread$id,
                              mappedReads = totalread$mappedReads,
                              enhancerReads = rawread.enhsum,
                              promoterReads = rawread.prmsum)
RPM_on_enhancers = RPM_on_enhancers %>%
  mutate(enhancerPercent = enhancerReads / mappedReads * 100,
         promoterPercent = promoterReads / mappedReads * 100)
RPM_on_enhancers[,-1] %>% colMeans
write.csv(RPM_on_enhancers, "readcount/enhancerReadSum.csv", quote = F, col.names = T, row.names = F, sep = ",")
