setwd("~/WorkS/sl2665/procap_network")
source('rscript/general/heatmap.R')
library(dplyr)
library(tidyr)


# Read table
pro = read.table('readcount/procap.ambr.norm.txt',stringsAsFactors=F)
pro.var = read.table('readcount/procap.var.norm.txt',stringsAsFactors=F)

# Locate targets
bcl2 = c(60790579, 60986613)
bcl2 = c(bcl2[1]*1.1-bcl2[2]*0.1,bcl2[2]*1.1-bcl2[1]*0.1)
slfn5 = c(33568000, 33573000)

bcl2.range = pro[,1]=="chr18" & pro[,2]>bcl2[1] & pro[,2]<bcl2[2]
slfn5.range = pro[,1]=="chr17" & pro[,2]>slfn5[1] & pro[,2]<slfn5[2]
sum(bcl2.range)

slfn5.mat = pro[slfn5.range,-(1:2)]
bcl2.mat = pro[bcl2.range,-(1:2)]

# Make a multipanel scatterplot for slfn5
pdf("pdf/Fig1/Fig1a.pdf",width=3.6,height=3.6)
par(mfrow=c(5,5))
par(cex=0.6)
par(mar=c(0.1,0.1,0.1,0.1))
par(oma=c(3,4,2,1))
par(tcl=-0.25)
par(mgp=c(1.5,.2,0))
for(i in 1:5) for(j in 1:5) {
	plot(log10(unlist(slfn5.mat[i,])+1),log10(unlist(slfn5.mat[j,])+1),pch=16,cex=0.5,xlab=NULL,ylab=NULL,axes=F)
	box(lwd=0.5)
}
dev.off()

# Make a correlation heatmap for bcl2
pdf("pdf/Fig1/Fig1b.pdf",width=3.6,height=4)
par(mar=c(3,2,2,1),mgp=c(1.5,0.2,0))
plot(0,0,type='n',xlim=bcl2,ylim=bcl2,axes=F,xaxs='i',yaxs='i',xlab="",ylab="")
# Heatmap
cor.mat = cor(t(bcl2.mat))
for(i in 1:nrow(cor.mat)) cor.mat[i,i] = 0
cor.mat = matrix(rep(cor.mat,each=20),ncol=ncol(cor.mat))
cor.mat = matrix(apply(cor.mat,2,rep,20),nrow=nrow(cor.mat))
heatmap.1col(cor.mat,xlim=bcl2,ylim=bcl2,zlim=c(-1,1))
range = 1000
heatmap.1col(matrix(-100:100/100,nrow=1),xlim=c(bcl2[2]-range*0.3,bcl2[2]),ylim=c(bcl2[1]-range*0.08,bcl2[1]-range*0.04),zlim=c(-1,1), lwd = 1.5)
# Lines from TREs to heatmap
n = nrow(bcl2.mat)
pos = pro[bcl2.range,2]
heatmappos = (1:n-0.5)*(bcl2[2]-bcl2[1])/n + bcl2[1]
range = bcl2[2] - bcl2[1]
for(i in 1:n) lines(c(pos[i],heatmappos[i]),c(bcl2[2]+0.05*range,bcl2[2]),lwd=1.5,xpd=T)
for(i in 1:n) lines(c(bcl2[1]-0.05*range,bcl2[1]),bcl2[1]+bcl2[2]-c(pos[i],heatmappos[i]),lwd=1.5,xpd=T)

dev.off()

#####
# Add Hi-C contact data

# Convert coordinates to hg38 version

write.table(data.frame(chr = c("chr17", "chr18"),
                       rbind(slfn5, floor(bcl2))),
            "HiC/tTRE/example.coord.txt",
            row.names = F, col.names = F, quote = F)

# Liftover to hg38 and read the coordinates
example.coord = read.table("HiC/tTRE/example.hg38.txt", col.names = c("chr", "start", "end"))

# write bash scripts for extracting HiC data from the coordinates
bash.command = c(paste0("cooler dump -b -t pixels -f --join -r ",
                        example.coord$chr[1], ":",
                        floor(example.coord$start[1]/1000), "k-",
                        floor(example.coord$end[1]/1000),
                        "k HiC/4DN/4DNFIXP4QG5B.mcool::/resolutions/1000 > HiC/4DN/example_contacts.txt"),
                 paste0("cooler dump -b -t pixels -f --join -r ",
                        example.coord$chr[2], ":",
                        floor(example.coord$start[2]/1000), "k-",
                        floor(example.coord$end[2]/1000),
                        "k HiC/4DN/4DNFIXP4QG5B.mcool::/resolutions/1000 >> HiC/4DN/example_contacts.txt"))
                 
writeLines(bash.command, "HiC/getExampleContacts.sh")

example.contacts = read.table("HiC/4DN/example_contacts.txt", header = F,
                          col.names = c("chr1", "pos1", "pos1e", "chr2", "pos2", "pos2e",
                                        "count", "contact"), fill = T)

slfn5.hiC.mat = example.contacts %>%
  filter(chr1 == "chr17") %>%
  select(pos1, pos2, contact) %>%
  spread(key = pos2, value = contact, fill = 0) %>%
  mutate(pos1 = as.numeric(pos1)) %>%
  arrange(pos1)

pdf("pdf/Fig1/FigS1.HiC.slfn5.pdf", width = 3.6, height = 4)
par(mar=c(3,2,2,1),mgp=c(1.5,0.2,0))
plot(0,0,type='n',xlim=slfn5,ylim=slfn5,axes=F,xaxs='i',yaxs='i',xlab="",ylab="")
# Heatmap
hiC.mat = as.matrix(slfn5.hiC.mat[,-1])
hiC.mat0 = hiC.mat
hiC.mat[hiC.mat==0] = min(hiC.mat[hiC.mat>0])/2
hiC.mat = matrix(rep(hiC.mat,each=40), ncol=ncol(hiC.mat))
hiC.mat = matrix(apply(hiC.mat, 2, rep,40), nrow=nrow(hiC.mat))
hiC.mat = log10(hiC.mat)
heatmap.1col(hiC.mat, xlim=slfn5, ylim=slfn5, zlim=c(min(hiC.mat)-2.5, max(hiC.mat)+0.5))
# Display HiC contact frequency
matpos.x = (1:ncol(hiC.mat0) - 0.5) / ncol(hiC.mat0) * (slfn5[2] - slfn5[1]) + slfn5[1]
matpos.y = rep(rev(matpos.x), ncol(hiC.mat0))
matpos.x = rep(matpos.x, each = ncol(hiC.mat0))
hiC.mat0 = sprintf("%0.3f", hiC.mat0)
text(matpos.x, matpos.y, labels = hiC.mat0)
dev.off()

bcl2.hiC.mat = example.contacts %>%
  filter(chr1 == "chr18") %>%
  select(pos1, pos2, contact) %>%
  spread(key = pos2, value = contact, fill = 0) %>%
  mutate(pos1 = as.numeric(pos1)) %>%
  arrange(pos1)


pdf("pdf/Fig1/FigS1.HiC.bcl2.pdf", width = 3.6, height = 4)
par(mar=c(3,2,2,1),mgp=c(1.5,0.2,0))
plot(0,0,type='n',xlim=bcl2,ylim=bcl2,axes=F,xaxs='i',yaxs='i',xlab="",ylab="")
# Heatmap
hiC.mat = as.matrix(bcl2.hiC.mat[,-1])
hiC.mat[hiC.mat==0] = min(hiC.mat[hiC.mat>0])/2
hiC.mat = matrix(rep(hiC.mat,each=5), ncol=ncol(hiC.mat))
hiC.mat = matrix(apply(hiC.mat, 2, rep,5), nrow=nrow(hiC.mat))
hiC.mat = log10(hiC.mat)
heatmap.1col(hiC.mat, xlim=bcl2, ylim=bcl2, zlim=c(min(hiC.mat)-2, max(hiC.mat)-1))
for(i in 1:n) lines(c(heatmappos[i], pos[i]),c(bcl2[2]+0.05*range,bcl2[2]),lwd=1.5,xpd=T)
dev.off()

