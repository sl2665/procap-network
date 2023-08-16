setwd("~/WorkS/sl2665/procap_network")
library(dplyr)
library(tidyr)


# Plot correlation scatterplot along the distance
pro = read.table('readcount/procap.var.norm.txt',stringsAsFactors=F)
pro.near = read.table('window/procap.window.all-all.1M.txt',stringsAsFactors=F)
pro.window = make.window.list(pro.near)
pro.cor  = cor.pro(pro.window, pro, pro)

pair.cor = data.frame(t1 = rep(1:length(pro.window), unlist(lapply(pro.window, length))),
                         t2 = unlist(pro.window),
                        cor = unlist(pro.cor$scd$corr), pval = unlist(pro.cor$scd$pval))

# Reattach hg19 positions in number format
pair.cor = pair.cor %>%
  left_join(pro.coord %>%
              mutate(t1 = 1:nrow(pro.coord),
                     pos1 = pos) %>%
              select(t1, chr, pos1),
            by = "t1") %>%
  left_join(pro.coord %>%
              mutate(t2 = 1:nrow(pro.coord),
                     pos2 = pos) %>%
              select(t2, pos2),
            by = "t2")
pair.cor = pair.cor %>%
  select(t1, t2, chr, pos1, pos2, pval, cor)
pair.cor$npos1 = conv.pos(pair.cor[, c(3,4)])
pair.cor$npos2 = conv.pos(pair.cor[, c(3,5)])
xlim = c(0,1000*2^(1:9)) 
pair.cor = pair.cor %>%
  mutate(dist = npos2 - npos1) %>%
  mutate(distBin = cut(dist/1000,xlim/1000))

#################
# Use mcool matrix from 4D nucleome dataset and HiContacts package

# Write ttRE coordinates in hg19
library(dplyr)
library(tidyr)
pro.bed = pro.coord %>%
  mutate(pos2 = pos + 1, id = paste0("t", sprintf("%05d", 1:nrow(pro.bed))))
write.table(pro.bed, "HiC/tTRE/hg19.bed", quote = F, sep = "\t", col.names = F, row.names = F)

# Change to hg38 coordinates
system("~/Work/tools/liftOver/liftOver HiC/tTRE/hg19.bed ~/Work/tools/liftOver/hg19ToHg38.over.chain.gz HiC/tTRE/hg38.bed HiC/tTRE/unmapped.bed")

# Round up to 5 kb resolution
pro.bed.hg38 = read.table("HiC/tTRE/hg38.bed", header = F, col.names = c("chr", "pos", "pos2", "id")) %>%
  mutate(pos2 = floor(pos/5000)*5000)

# Extract Hi-C contract matrix +/- 500 kb region in 5kb around the tTRES
# Generate a bash script using Cooler in the syntax example below
# cooler dump -b -t pixels -f --join -r chr3:20000k-20005k -r2 chr3:19000k-21000k HiC/4DNFIXP4QG5B.mcool::/resolutions/5000 > tmp

bash.command = sapply(1:nrow(pro.bed.hg38),
                      function(i) {
                        chr = pro.bed.hg38$chr[i]
                        pos = pro.bed.hg38$pos2[i]
                        return(paste0("cooler dump -b -t pixels -f --join -r ",
                                      chr, ":",
                                      pos/1000, "k-",
                                      pos/1000 + 5, "k -r2 ",
                                      chr, ":",
                                      pos/1000 - 500, "k-",
                                      pos/1000 + 500, "k HiC/4DN/4DNFIXP4QG5B.mcool::/resolutions/5000 >> HiC/4DN/tTREcontacts.txt"))})
writeLines(bash.command, "HiC/getContacts.sh")

# Run the getContacts.sh bash in the system console

# Read tTREcontacts.txt file

HiC.contacts = read.table("HiC/4DN/tTREcontacts.txt", header = F,
                          col.names = c("chr1", "pos1", "pos1e", "chr2", "pos2", "pos2e",
                                        "count", "contact"), fill = T)
HiC.contacts[is.na(HiC.contacts)] = 0

HiC.contacts = HiC.contacts %>%
  distinct(chr1, pos1, pos2, .keep_all = T)

# Retrieve coordinates of tTRE pairs in pro.window as a dataframe, and indicate pairs
pairlist.id = Reduce(rbind,
       lapply(1:length(pro.window),
       function(i) {
         if(length(pro.window[[i]]>0)) return(data.frame(t1 = i,t2 = pro.window[[i]]))
         else return(data.frame(t1=c(), t2=c()))}))
                                              
pairlist.id = data.frame(t1 = rep(1:length(pro.window), unlist(lapply(pro.window, length))),
                         t2 = unlist(pro.window))
pairlist.id = pairlist.id %>%
  mutate(i1 = paste0("t", sprintf("%05d", t1)),
         i2 = paste0("t", sprintf("%05d", t2)))

# Retrieve hg38 coordinates
pairlist.id = pairlist.id %>%
  inner_join(pro.bed.hg38 %>%
               mutate(pos1 = pos2, i1 = id) %>%
               select(chr, pos1, i1),
             by = "i1") %>%
  inner_join(pro.bed.hg38 %>%
               mutate(pos2, i2 = id) %>%
               select(pos2, i2),
             by = "i2")

# Extract HiC contact values
pair.HiC = pairlist.id %>%
  left_join(HiC.contacts %>%
               mutate(chr = chr1) %>%
               select(chr, pos1, pos2, count, contact),
             by = c("chr", "pos1", "pos2"))
pair.HiC[is.na(pair.HiC)] = 0

# Reattach hg19 positions in number format
pair.HiC = pair.HiC %>%
  left_join(pro.coord %>%
              mutate(t1 = 1:nrow(pro.coord),
                     pos1.hg19 = pos) %>%
              select(t1, pos1.hg19),
            by = "t1") %>%
  left_join(pro.coord %>%
              mutate(t2 = 1:nrow(pro.coord),
                     pos2.hg19 = pos) %>%
              select(t2, pos2.hg19),
            by = "t2")
pair.HiC = pair.HiC %>%
  mutate(pos1 = pos1.hg19, pos2 = pos2.hg19) %>%
  select(t1, t2, chr, pos1, pos2, count, contact)

# Function to convert chromosomal positions to numbers
conv.pos = function(p) {
  chr = substring(p[,1],4)
  chr[chr=="X"] = "23"
  chr[chr=="Y"] = "24"
  chr = as.numeric(chr)
  chr[is.na(chr)] = 0
  return(chr*1000000000 + p[,2])
}

# convert HiC coordinates to numeric positions
pair.HiC$npos1 = conv.pos(pair.HiC[, c(3,4)])
pair.HiC$npos2 = conv.pos(pair.HiC[, c(3,5)])

# Calculate distance between the two
xlim = c(0,1000*2^(1:9)) 
pair.HiC = pair.HiC %>%
  mutate(dist = npos2 - npos1) %>%
  mutate(distBin = cut(dist/1000,xlim/1000))
pair.HiC$contact[pair.HiC$contact==0] = min(pair.HiC$contact[pair.HiC$contact>0])/2

# Load ChIP-seq dataset
ChIPnames = c("ATF2","ATF3","BATF","BCL11A","BCL3","BCLAF1","BHLHE40","BRCA1",
              "CEBPB","CHD1","CHD2","CTCF","E2F4","EBF1","EGR1","ELF1","ELK1","EP300",
              "ETS1","EZH2","FOS","FOXM1","GABPA","IKZF1","IRF4","JUND","MAX","MAZ",
              "MEF2A","MEF2C","MTA3","MXI1","MYC","NFATC1","NFE2","NFIC","NFYA","NFYB",
              "NR2C2","NRF1","PAX5","PBX3","PML","POLR2A","POLR3G","POU2F2","RAD21",
              "RCOR1","RELA","REST","RFX5","RUNX3","RXRA","SIN3A","SIX5","SMC3","SP1",
              "SPI1","SRF","STAT1","STAT3","STAT5A","TAF1","TBL1XR1","TBP","TCF12","TCF3",
              "USF1","USF2","WRNIP1","YY1","ZBTB33","ZEB1","ZNF143","ZNF274","ZZZ3")
ct = lapply(ChIPnames, function(x) data.frame(TF=x, read.table(paste("ChIP/Tfbs.",x,".txt",sep=""),stringsAsFactors=F,header=F, col.names = c("chr", "pos"))))

ChIPnames[18] = "P300"
names(ct) = ChIPnames


# Convert ChIP-seq peaks to numeric positions
ctpos = Reduce(bind_rows, ct)
ctpos = ctpos %>%
  mutate(TF = ifelse(TF == "EP300", "P300", TF))

ctpos$npos = conv.pos(ctpos[, 2:3])
ctpos = ctpos %>%
  arrange(npos) %>%
  select(TF, npos)

# New function to determine intersecting 
count.in.tf = function(npos1, npos2, tf) {
  p = ctpos %>%
    filter(TF == tf) %>%
    select(npos) %>% unlist
  i1 = findInterval(npos1, p)
  i2 = findInterval(npos2, p)
  return(i2-i1)
}

count.ol.tf = function(p1, p2, tf, range=500) {
  p = ctpos %>%
    filter(TF == tf) %>%
    select(npos) %>% unlist
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

# Function to generate intersect boxplot
plot.HiC.4DN = function(tf = "RAD21", type = "intersect", tf.cutoff = c(1,2),
                        ylim = c(-4.5, -0.5),
                        col = "#ffffff", alpha = 1,
                        file = "test.pdf") {
  countfunc = count.in.tf
  if(type == "overlap") countfunc = count.ol.tf
  pair.HiC.tfin <<- pair.HiC %>%
    mutate(tf.in = countfunc(npos1, npos2, tf)) %>%
    mutate("{tf}" := ifelse(tf.in <= tf.cutoff[1], 0, 
                          ifelse(tf.in >= tf.cutoff[2], 1, 2))) %>%
    select(dist, distBin, {{tf}}, contact, count) %>%
    filter(get({{tf}}) <=1) %>%
    drop_na()
  tf0 = pair.HiC.tfin %>%
    filter(get({{tf}}) == 0) %>%
    select(dist, contact)
  tf1 = pair.HiC.tfin %>%
    filter(get({{tf}}) == 1) %>%
    select(dist, contact)
  n=plot.scatter.dist2(tf0$dist, log10(tf0$contact), tf1$dist, log10(tf1$contact),
                       file=file,ylim=ylim,
                       xlim=xlim, xlab="", col = col, alpha = alpha,
                       ylab = "Hi-C contact (log)")
  if(type == "overlap") {
    legend("top", inset=c(0,-0.2), legend = c(as.expression(bquote(.(tf)~"(-)" )),
                                              as.expression(bquote(.(tf)~"(+)" )) ),
                               pch=22, pt.cex=1.5,lty=0, pt.bg=alpha(c("#ffffff", colorRampPalette(col)(1)), alpha),xpd=T,cex=0.8, ncol=2, bty='n')
  } else {
    legend("top", inset=c(0,-0.2), legend=c(as.expression(bquote(.(tf) <= .(tf.cutoff[1]))),as.expression(bquote(.(tf) >= .(tf.cutoff[2])))),
         pch=22, pt.cex=1.5,lty=0, pt.bg=alpha(c("#ffffff", colorRampPalette(col)(1)), alpha),xpd=T,cex=0.8, ncol=2, bty='n')}
  print(n)
}

majorTF = c("CTCF", "P300", "RAD21")
cols = list(color$greenpa(10), color$orpa(10), color$pinkpa(10))
names(cols) = majorTF


for(i in majorTF) plot.HiC.4DN(i, file = "pdf/Fig4/Fig4a.HiC.4DN.intersect.pdf", col = cols[[i]])
dev.off()

for(i in majorTF) plot.HiC.4DN(i, "overlap", c(0, 1), file = "pdf/Fig4/Fig4b.HiC.4DN.overlap.pdf", col = cols[[i]])
dev.off()


pdf("pdf/Fig4/HiC.4DN.intersect.pdf")
for(i in ChIPnames) plot.HiC.4DN(i)
dev.off()

pdf("pdf/Fig4/HiC.4DN.overlap.pdf")
for(i in ChIPnames) plot.HiC.4DN(i, "overlap", c(0, 1))
dev.off()


# Function to calculate AUC
source("rscript/functions.R")

cdc.tf = function(tf = "CTCF", type = "intersection", tf.cutoff = c(1, 2),
                      median = F, max.dist = 200000, ylim = NULL, spar = 0.7,
                      pairdata = pair.HiC, log.transform = T) {
  countfunc = count.in.tf
  if(type == "occupancy") countfunc = count.ol.tf
  pair.tfin = pairdata %>%
    mutate(tf.in = countfunc(npos1, npos2, tf)) %>%
    mutate("{tf}" := ifelse(tf.in <= tf.cutoff[1], 0, 
                            ifelse(tf.in >= tf.cutoff[2], 1, 2))) %>%
    filter(get({{tf}}) <=1) %>%
    drop_na()
  tf.0 = pair.tfin %>%
    filter(get({{tf}}) == 0)
  tf.1 = pair.tfin %>%
    filter(get({{tf}}) == 1)
  
  yname = colnames(pairdata)[7]
  if(log.transform) {
    qc.0 = quantile.curve(x = tf.0$dist, y = log10(tf.0[,7]), n = 200,
                        xlim = c(0, max.dist), p.val = c(0.5, 0.95), spar = spar)
    qc.1 = quantile.curve(x = tf.1$dist, y = log10(tf.1[,7]), n = 200,
                        xlim = c(0, max.dist), p.val = c(0.5, 0.95), spar = spar)
  } else {
    qc.0 = quantile.curve(x = tf.0$dist, y = tf.0[,7], n = 200,
                          xlim = c(0, max.dist), p.val = c(0.5, 0.95), spar = spar)
    qc.1 = quantile.curve(x = tf.1$dist, y = tf.1[,7], n = 200,
                          xlim = c(0, max.dist), p.val = c(0.5, 0.95), spar = spar)    
  }
  qc.0a = qc.0$r
  qc.1a = qc.1$r
  qc.0s = qc.0$s
  qc.1s = qc.1$s
  ycol = 3
  if(median) ycol = 2
  
  plotdata = rbind(
    data.frame(Distance = qc.0a[,1],
               Class = 1) %>%
      mutate("{yname}" := qc.0a[,ycol]),
    data.frame(Distance = qc.1a[,1],
               Class = 0) %>%
      mutate("{yname}" := qc.1a[,ycol])) %>%
    mutate(Class = as.factor(Class)) %>%
    mutate(type = 1)
  if(F) {
  plotdata = rbind(
    data.frame(Distance = qc.0s[,1],
               Class = 1) %>%
      mutate("{yname}" := qc.0s[,ycol]),
    data.frame(Distance = qc.1s[,1],
               Class = 0) %>%
      mutate("{yname}" := qc.1s[,ycol])) %>%
    mutate(Class = as.factor(Class)) %>%
    mutate(type = 0) %>%
    rbind(plotdata) %>%
    mutate(type = as.factor(type))
  }
  
  g = ggplot(plotdata, aes(x = Distance/1000, y = get(yname), col = Class)) +
    geom_line() +
    theme_bw() +
    coord_cartesian(xlim = c(0, max.dist/1000), ylim = ylim) +
    ggtitle(paste(tf, type)) +
    scale_color_manual(values = c("#E09810", "#084FA0")) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
          legend.position = "none",
          axis.line = element_line(linewidth = 0)) +
    ylab(yname) + xlab("Distance (kb)") 
  
  print(g)
  return(cal.AUC(x = qc.1s[,1], y = qc.1s[,ycol], max.dist = max.dist) -
           cal.AUC(x = qc.0s[,1], y = qc.0s[,ycol], max.dist = max.dist))
}


# Top15, bottom 15 show

pdf("pdf/Fig4/HiC.4DN.cdc.intersect.pdf", height = 4, width = 4)
dAUC.in = data.frame(TF = c(), dAUC = c())
for(i in ChIPnames) {
  dAUC.in <<- rbind(dAUC.in, data.frame(TF = i,
                      dAUC = cdc.tf(i, median = T, pairdata = pair.HiC, ylim = c(-3, -1))))
}
dev.off()

pdf("pdf/Fig4/HiC.4DN.cdc.overlap.pdf", height = 4, width = 4)
dAUC.ol = data.frame(TF = c(), dAUC = c())
for(i in ChIPnames) {
  dAUC.ol <<- rbind(dAUC.ol, data.frame(TF = i, 
                                          dAUC = cdc.tf(i, type = "occupancy",
                                                            tf.cutoff = c(0, 1), median = T, pairdata = pair.HiC, ylim = c(-3, -1))))
}
dev.off()

#
ChIPnames1 =c("CTCF", "P300", "RAD21", "EZH2")

pdf("pdf/Fig4/FigS4.HiC.4DN.cdc.intersect.pdf", height = 3.6, width = 3.4)
dAUC.in = data.frame(TF = c(), dAUC = c())
for(i in ChIPnames1) {
  dAUC.in <<- rbind(dAUC.in, data.frame(TF = i,
                                        dAUC = cdc.tf(i, median = T, pairdata = pair.HiC, ylim = c(-3, -1))))
}
dev.off()

pdf("pdf/Fig4/FigS4.HiC.4DN.cdc.overlap.pdf", height = 3.6, width = 3.4)
dAUC.ol = data.frame(TF = c(), dAUC = c())
for(i in ChIPnames1) {
  dAUC.ol <<- rbind(dAUC.ol, data.frame(TF = i, 
                                        dAUC = cdc.tf(i, type = "occupancy",
                                                      tf.cutoff = c(0, 1), median = T, pairdata = pair.HiC, ylim = c(-3, -1))))
}
dev.off()


# Performing this kind of AUC analysis for median correlation coefficients
pdf("pdf/Fig4/cor.4DN.cdc.intersect.pdf", height = 4, width = 4)
cAUC.in = data.frame(TF = c(), dAUC = c())
for(i in ChIPnames) {
  cAUC.in <<- rbind(cAUC.in, data.frame(TF = i,
                                        dAUC = cdc.tf(i, median = T, pairdata = pair.cor, log.t = F, ylim = c(0, 0.7))))
}
dev.off()

pdf("pdf/Fig4/cor.4DN.cdc.overlap.pdf", height = 4, width = 4)
cAUC.ol = data.frame(TF = c(), dAUC = c())
for(i in ChIPnames) {
  cAUC.ol <<- rbind(cAUC.ol, data.frame(TF = i, 
                                        dAUC = cdc.tf(i, type = "overlap",
                                                      tf.cutoff = c(0, 1), median = T, pairdata = pair.cor, ylim = c(0, 0.7), log.t = F)))
}
dev.off()

pdf("pdf/Fig4/c95.4DN.cdc.intersect.pdf", height = 4, width = 4)
c95AUC.in = data.frame(TF = c(), dAUC = c())
for(i in ChIPnames) {
  c95AUC.in <<- rbind(c95AUC.in, data.frame(TF = i,
                                        dAUC = cdc.tf(i, median = F, pairdata = pair.cor, log.t = F, ylim = c(0, 0.7))))
}
dev.off()

pdf("pdf/Fig4/c95.4DN.cdc.overlap.pdf", height = 4, width = 4)
c95AUC.ol = data.frame(TF = c(), dAUC = c())
for(i in ChIPnames) {
  c95AUC.ol <<- rbind(c95AUC.ol, data.frame(TF = i, 
                                        dAUC = cdc.tf(i, type = "overlap",
                                                      tf.cutoff = c(0, 1), median = F, pairdata = pair.cor, ylim = c(0, 0.7), log.t = F)))
}
dev.off()

AUC.all = Reduce(bind_rows,
       lapply(c("dAUC.in", "dAUC.ol", "cAUC.in", "cAUC.ol", "c95AUC.in", "c95AUC.ol"),
              function(x) get(x) %>%
                mutate(type = x))) %>%
  spread(type, dAUC) %>%
  drop_na %>%
  filter(!(TF %in% c("ATF3",
                     "BRCA1",
                     "NFE2",
                     "NR2C2",
                     "POLR3G",
                     "STAT1",
                     "ZBTB33",
                     "ZNF274",
                     "ZZZ3")))  # removed TF data that has less than 2,000 pairs in one of the TF intersect or overlap classes (200 pair per interval, less than 10 interval)


pdf("pdf/Fig3/Fig3S.AUC_intersect_median_vs_95.pdf")
ggplot(data = AUC.all, aes(x = c95AUC.in, y = cAUC.in, label = TF)) +
  geom_point(col = "black") +
  geom_text(vjust = -0.5, position=position_jitter(width=0.02,height=0)) +
  xlab("AUC (95 percentile)") +
  ylab("AUC (median)") +
  ggtitle("TF intersect") +
  theme_bw() +
  annotate("text", label = paste0("r = ", cor(AUC.all$cAUC.in, AUC.all$c95AUC.in)),
                  x= 0, y = -0.1)
dev.off()

pdf("pdf/Fig3/Fig3S.AUC_occupy_median_vs_95.pdf")
ggplot(data = AUC.all, aes(x = c95AUC.ol, y = cAUC.ol, label = TF)) +
  geom_point(col = "black") +
  geom_text(vjust = -0.5) +
  xlab("AUC (95 percentile)") +
  ylab("AUC (median)") +
  ggtitle("TF occupancy") +
  theme_bw() +
  annotate("text", label = paste0("r = ", cor(AUC.all$cAUC.ol, AUC.all$c95AUC.ol)),
           x= 0.05, y = -0.05)
dev.off()

pdf("pdf/Fig3/Fig3S.intersect_overlap.pdf")
ggplot(data = AUC.all, aes(x = cAUC.in, y = cAUC.ol, label = TF)) +
  geom_point(col = "black") +
  geom_text(vjust = -0.5) +
  xlab("Intersect") +
  ylab("Occupancy") +
  ggtitle("Correlation (intersect - occupancy)") +
  theme_bw()
dev.off()


### Figures to compare contact vs correlation
TF.show = data.frame(TF = c("CTCF", "RAD21", "P300"),
               col = c(color$green(3)[2], color$pinkpa(3)[2], color$orange(3)[2]))
TF.show = TF.show %>%
  rbind(data.frame(TF = c("EZH2", "#SMC3", "#YY1"),
                   col = "#000000"))
TF.col = TF.show$col
names(TF.col) = TF.show$TF
TF.col = c(TF.col, "all" = "#000000")

AUC.all = AUC.all %>%
  mutate(TF = ifelse(TF == "EP300", "P300", TF))

AUC.new = AUC.all[, 1:7] %>%
  left_join(TF.show, by = "TF")
AUC.new = AUC.new %>%
  mutate(TF.show = ifelse(is.na(AUC.new$col), "", TF)) %>%
  mutate(col = ifelse(is.na(AUC.new$col), "all", TF))

pdf("pdf/Fig4/Fig4c.Cor_vs_HiC_intersect.pdf", width = 4, height = 4)
ggplot(data = AUC.new, aes(x = cAUC.in, y = dAUC.in, label = TF.show, col = col)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
  geom_point() +
  geom_text(vjust = -0.5, hjust = 0.2) +
  scale_color_manual(values = TF.col) +
  xlab("Coexpression (PRO-cap)") +
  ylab("Contact frequency (Hi-C)") +
  ggtitle("TF intersection") +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        legend.position = "none",
        axis.line = element_line(linewidth = 0))
dev.off()

pdf("pdf/Fig4/Fig4d.Cor_vs_HiC_overlap.pdf", width = 4, height = 4)
ggplot(data = AUC.new, aes(x = cAUC.ol, y = dAUC.ol, label = TF.show, col = col)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
  geom_point() +
  geom_text(vjust = -0.5, hjust = 0.7) +
  scale_color_manual(values = TF.col) +
  xlab("Coexpression (PRO-cap)") +
  ylab("Contact frequency (Hi-C)") +
  ggtitle("TF occupancy") +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        legend.position = "none",
        axis.line = element_line(linewidth = 0))
dev.off()



########
# Make some HiC plots for Fig 2

head(pair.HiC)
head(rna.pos)

pair.HiC.scatter = pair.HiC %>%
  filter(contact > min(pair.HiC$contact))
# Plot all Hi-C interactions
plot.scatter.dist(pair.HiC.scatter$dist,
                  log10(pair.HiC.scatter$contact), xlim=c(0,1000*2^(0:9)), box=T,
                  ylim=c(-4, -0.5),xlab="All tTRE contacts",
                  ylab = "Hi-C contact frequency (log)", file="pdf/Fig2/FigS2a.all_HiC.pdf",
                  col = color$parula(10), alpha = 0.3)

HiC.cdc = quantile.curve(pair.HiC.scatter$dist,pair.HiC.scatter$contact,n=1000,p.val=c(0.05,0.5,0.95),spar=0)
plot.scatter.dist(pair.HiC.scatter$dist, pair.HiC.scatter$contact,
                  cdc = HiC.cdc$r,
                  xlim=c(0,500000), file="pdf/Fig2/FigS2a1.all_HiC.500kb.pdf",
                  xlab = "Distance (kb)", ylab = "Hi-C contact frequency", grid = F,
                  ylim = c(0, 0.05),
                  col = color$parula(10), alpha = 0.75)

HiC.log.cdc = quantile.curve(pair.HiC.scatter$dist,log10(pair.HiC.scatter$contact),n=1000,p.val=c(0.05,0.5,0.95),spar=0)
plot.scatter.dist(pair.HiC.scatter$dist, log10(pair.HiC.scatter$contact),
                  cdc = HiC.log.cdc$r,
                  xlim=c(0,500000), file="pdf/Fig2/FigS2a2.all_HiC.500kb_log.pdf",
                  xlab = "Distance (kb)", ylab = "Hi-C contact frequency (log)",
                  ylim = c(-4.5, -0.5),
                  col = color$parula(10), alpha = 0.75)

# Plot HiC-mRNA TSS interaction with directional content
# Identify all promoters within 500 bp from mRNA TSSs and label their gene strand

pro.npos = data.frame(ppos = conv.pos(pro.pro[,1:2]))
rna.npos = data.frame(rpos = conv.pos(rna.pos[,1:2]), strand = rna.pos[,6]) %>%
  arrange(rpos)

pro.npos = pro.npos %>%
  mutate(rpos1 = rna.npos$rpos[findInterval(ppos, rna.npos$rpos)],
         strand1 = rna.npos$strand[findInterval(ppos, rna.npos$rpos)],
         rpos2 = rna.npos$rpos[findInterval(ppos, rna.npos$rpos)+1],
         strand2 = rna.npos$strand[findInterval(ppos, rna.npos$rpos)+1]) %>%
  mutate(dist1 = ppos - rpos1,
         dist2 = rpos2 - ppos) %>%
  filter(dist1 < 500 | dist2 < 500)

pro.stranded = rbind(
  pro.npos %>%
    filter(dist1 < 500) %>%
    mutate(strand = strand1) %>%
    select(ppos, strand),
  pro.npos %>%
    filter(dist2 < 500) %>%
    mutate(strand = strand2) %>%
    select(ppos, strand))
  
# Promoters in pos1
pair.HiC1 = pair.HiC %>%
  inner_join(pro.stranded %>%
               mutate(npos1 = ppos) %>%
               select(npos1, strand)) %>%
  mutate(dist = ifelse(strand == "+",
                       dist, -dist))
# Promoters in pos2
pair.HiC2 = pair.HiC %>%
  inner_join(pro.stranded %>%
               mutate(npos2 = ppos) %>%
               select(npos2, strand)) %>%
  mutate(dist = ifelse(strand == "-",
                       dist, -dist))
# All promoter-enhancer Hi-C contacts with regard to the orientation
strand.HiC = rbind(pair.HiC1 %>% select(dist, contact),
                   pair.HiC2 %>% select(dist, contact))

# Plot stranded Hi-C
HiC.stranded.cdc = quantile.curve(strand.HiC$dist,log10(strand.HiC$contact),n=1000,p.val=c(0.05,0.5,0.95),spar=0)
plot.scatter.dist(strand.HiC$dist, log10(strand.HiC$contact),
                  cdc = HiC.stranded.cdc$s,
                  xlim=c(-200000,200000), file="pdf/Fig2/FigS3a1.stranded_HiC.500kb_log.pdf",
                  xlab = "Relative position (kb)", ylab = "Hi-C contact frequency (log)",
                  ylim = c(-4.5, -0.5))

plot.scatter.dist(strand.HiC$dist,log10(strand.HiC$contact),
                  xlim=c(-1000*2^(2*(5:0)),1000*2^(2*(0:5))),
                  box=T,ylim=c(-6.5,-0.5),
                  file="pdf/Fig2/FigS3a2.stranded_HiC.boxplot.pdf",
                  xlab="Relative position (kb)")


# Read RNA position and HiC around the TSSs

rna.pos=read.table('window/rnaseq.pos.bed', stringsAsFactors=F) %>%
  mutate(V4 = paste0("r", sprintf("%05d", 1:nrow(rna.pos)))) %>%
  mutate(V2 = ifelse(V2<0, 1, V2)) %>%
  mutate(V3 = V2 + 1)

write.table(rna.pos, "HiC/tTRE/rna.hg19.bed", row.names = F, col.names = F, quote = F)
rna.hg38 = read.table("HiC/tTRE/rna.hg38.txt", header = F, col.names = c("chr", "start", "end", "id", "score", "strand") ) %>%
  mutate(end = floor(end/5000)*5000)
# Random position control
rnd.hg38 = read.table("HiC/tTRE/hg38.random.txt", header = F, col.names = c("chr", "start", "end", "id", "score", "strand") ) %>%
  mutate(end = floor(end/5000)*5000)

# Extract Hi-C contract matrix +/- 500 kb region in 5kb around the tTRES
# Generate a bash script using Cooler in the syntax example below
# cooler dump -b -t pixels -f --join -r chr3:20000k-20005k -r2 chr3:19000k-21000k HiC/4DNFIXP4QG5B.mcool::/resolutions/5000 > tmp

bash.command = sapply(1:nrow(rna.hg38),
                      function(i) {
                        chr = rna.hg38$chr[i]
                        pos = rna.hg38$end[i]
                        return(paste0("cooler dump -b -t pixels -f --join -r ",
                                      chr, ":",
                                      pos/1000, "k-",
                                      pos/1000 + 5, "k -r2 ",
                                      chr, ":",
                                      pos/1000 - 1000, "k-",
                                      pos/1000 + 1000, "k HiC/4DN/4DNFIXP4QG5B.mcool::/resolutions/5000 >> HiC/4DN/RNAcontacts.txt"))})
writeLines(bash.command, "HiC/getRNAContacts.sh")

bash.command = sapply(1:nrow(rna.hg38),
                      function(i) {
                        chr = rnd.hg38$chr[i]
                        pos = rnd.hg38$end[i]
                        return(paste0("cooler dump -b -t pixels -f --join -r ",
                                      chr, ":",
                                      pos/1000, "k-",
                                      pos/1000 + 5, "k -r2 ",
                                      chr, ":",
                                      pos/1000 - 1000, "k-",
                                      pos/1000 + 1000, "k HiC/4DN/4DNFIXP4QG5B.mcool::/resolutions/5000 >> HiC/4DN/RNDcontacts.txt"))})
writeLines(bash.command, "HiC/getRandomContacts.sh")


# Read HiC contacts around the gene, with respect to the direction

RNA.contacts = read.table("HiC/4DN/RNAcontacts.txt", header = F,
                          col.names = c("chr1", "pos1", "pos1e", "chr2", "pos2", "pos2e",
                                        "count", "contact"), fill = T)
RNA.contacts[is.na(RNA.contacts)] = 0
RNA.contacts = RNA.contacts %>%
  distinct(chr1, pos1, pos2, .keep_all = T)

RNA.contacts.npos = RNA.contacts %>%
  mutate(npos1 = conv.pos(RNA.contacts[,c(1,2)]),
         npos2 = conv.pos(RNA.contacts[,c(1,5)])) %>%
  select(npos1, npos2, contact)

# Read random HiC contacts, used as the baseline
random.contacts = read.table("HiC/4DN/RNDcontacts.txt", header = F,
                          col.names = c("chr1", "pos1", "pos1e", "chr2", "pos2", "pos2e",
                                        "count", "contact"), fill = T)
random.contacts[is.na(random.contacts)] = 0
random.contacts = random.contacts %>%
  distinct(chr1, pos1, pos2, .keep_all = T)

random.contacts.npos = random.contacts %>%
  mutate(npos1 = conv.pos(random.contacts[,c(1,2)]),
         npos2 = conv.pos(random.contacts[,c(1,5)])) %>%
  select(npos1, npos2, contact) %>%
  mutate(dist = npos2 - npos1)

random.contacts.npos[is.na(random.contacts.npos)] = 0
View(random.contacts.npos)
random.contacts.distinct = random.contacts.npos %>%
  distinct(npos1, npos2)
View(random.contacts.distinct)



# Extract all TSS-tTRE contacts with respect to the gene strand
RNA.contacts.TSS = RNA.contacts.npos %>%
  inner_join(rna.hg38 %>%
               mutate(npos1 = conv.pos(rna.hg38[, c(1,3)]),
                      pos1 = start) %>%
               select(npos1, strand, pos1),
             by = "npos1")
View(RNA.contacts.TSS)

RNA.contacts.tre = RNA.contacts.TSS %>%
  inner_join(pro.bed.hg38 %>%
               mutate(npos2 = conv.pos(pro.bed.hg38[, c(1, 3)]),
                      pos2 = pos) %>%
               select(npos2, pos2),
             by = "npos2")
View(RNA.contacts.tre)

RNA.contacts.tre = RNA.contacts.tre %>%
  mutate(dist = ifelse(strand == "+",
                       pos2 - pos1,
                       pos1 - pos2))

RNA.strand = RNA.contacts.TSS %>%
  mutate(dist = ifelse(strand == "+",
                       npos2 - npos1,
                       npos1 - npos2))

# Plot Hi-C relative to mRNA genes
cdc.HiC.TSS_tre = quantile.curve(
  RNA.contacts.tre$dist,
  RNA.contacts.tre$contact,
  n = 1000,
  p.val = c(0.05, 0.5, 0.95),
  spar = 0.1)
cdc.cor.TSS_tre = quantile.curve(
  rna.cor$x,
  rna.cor$y,
  n = 1000,
  p.val = c(0.05, 0.5, 0.95),
  spar = 0.1)
cdc.HiC.TSS_all = quantile.curve(
  RNA.strand$dist,
  RNA.strand$contact,
  n = 1000,
  p.val = c(0.05, 0.5, 0.95),
  spar = 0.1)
cdc.HiC.random = quantile.curve(
  random.contacts.npos$dist,
  random.contacts.npos$contact,
  n = 1000,
  p.val = c(0.05, 0.5, 0.95),
  spar = 0.1)

# Quick check plot
plotdata = data.frame(cdc.HiC.random$r) %>%
  drop_na %>%
  mutate(Type = "Background", X1 = X1 - 2500) %>%
  bind_rows(data.frame(cdc.HiC.TSS_all$r) %>%
      drop_na %>%
      mutate(Type = "TSS-local", X1 = X1 - 2500)) %>%
  bind_rows(data.frame(cdc.HiC.TSS_tre$r) %>%
      drop_na %>%
      mutate(Type = "TSS-tTRE")) %>%
  mutate(Type = factor(Type))
colnames(plotdata) = c("dist", "Q5", "Median", "Q95", "Type")

plotdata2 = data.frame(cdc.cor.TSS_tre$r) %>%
  drop_na %>%
  mutate(Type = "Coexpression")
colnames(plotdata2) = c("dist", "Q5", "Median", "Q95", "Type")

cdc.interpolate = function(cdc, x, shift2 = -2500) {
  r = data.frame(cdc$r) %>% drop_na()
  return(r[findInterval(x, r[,1] + shift2), ])
}

cdc.diff = function(cdc1, cdc2, shift2 = 0) {
  r = data.frame(cdc1$r) %>% drop_na()
  r2 = r[seq(1, nrow(r), by = 2), ]
  b = cdc.interpolate(cdc2, r2[,1], shift2)
  y = rep(1:nrow(r2), each = 2)
  res = cbind(r[,1], r[, -1] - b[y,-1])
  colnames(res) = colnames(r)
  return(res)
}

plotdata3 = cdc.diff(cdc.HiC.TSS_tre, cdc.HiC.random, 0)
idx1 = findInterval(-7500, hiC.TSS_tre.norm[,1])
idx2 = findInterval(7500, hiC.TSS_tre.norm[,1])
plotdata3[idx1:(idx2+1),-1] = NA
colnames(plotdata3) = c("dist", "Q5", "Median", "Q95")
plotdata3$Type = "Hi-C (background\nsubtracted)"


pdf("pdf/Fig2/FigS2d.HiC_all.pdf", width = 4, height = 3)
ggplot(plotdata, aes(x = dist/1000, y = Median, col = Type)) +
  geom_line() +
#  geom_line(data = plotdata2, aes(x = dist/1000, y = Median, col = Type)) +
  coord_cartesian(xlim = c(-200, 200)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = c(0.8, 0.7)) +
  xlab("Position from TSS (kb)") +
  ylab("Median Hi-C contact") +
  scale_color_manual(values = color$adjust(c("#c0c0c0", color$greenpa(2)[2], color$japa(3)[2]), 0.95 , 1))
dev.off()

pdf("pdf/Fig2/FigS2d2.HiC_all_log.pdf", width = 4, height = 3)
ggplot(plotdata, aes(x = dist/1000, y = log10(Median), col = Type)) +
  geom_line() +
  #  geom_line(data = plotdata2, aes(x = dist/1000, y = Median, col = Type)) +
  coord_cartesian(xlim = c(-200, 200), ylim = c(-3.5, -1)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = c(0.8, 0.7)) +
  xlab("Position from TSS (kb)") +
  ylab("log Median Hi-C contact") +
  scale_color_manual(values = color$adjust(c("#c0c0c0", color$greenpa(2)[2], color$japa(3)[2]), 0.95 , 1))
dev.off()

pdf("pdf/Fig2/FigS2e.HiC_re.pdf", width = 4, height = 3)
ggplot(plotdata3, aes(x = dist/1000, y = Median*220, col = Type, fill = Type)) +
  geom_ribbon(aes(ymax = Median*220, ymin = -1)) +
  geom_ribbon(data = plotdata2, aes(x = dist/1000, ymax = Median*10, ymin = -1), col = "white", fill = "white") +
  geom_ribbon(data = plotdata2, aes(x = dist/1000, ymax = Median*10, col = Type, fill = Type, ymin = -1)) +
  xlim(-200, 200) +
  coord_cartesian(xlim = c(-185, 185), ylim = c(-0.1, 1.6)) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1.125),
        legend.position = c(0.75, 0.7),
        axis.line=element_line(linewidth = 0)) +
  xlab("Position from TSS (kb)") +
  ylab("Coexpression or contact (AU)") +
  scale_color_manual(values = color$adjust(color$bluered(10)[c(2,9)], 0.95 , 1)) +
  scale_fill_manual(values = color$adjust(color$bluered(10)[c(2,9)], 1 , 0.2)) 
dev.off()

# Plot Hi-C relative to mRNA genes
plot.scatter.dist(RNA.contacts.tre$dist,
                  log10(RNA.contacts.tre$contact),
                  cdc = RNA.contacts.cdc,
                  xlim=c(-1000*2^(2*(5:0)),1000*2^(2*(0:5))),box=T, grid = F,
                  ylim=c(-4.5,-0.5),file="pdf/Fig2/FigS2c.HiC.TRE-mRNA.pdf",xlab="tTRE-mRNA Hi-C contact",
                  ylab = "Hi-C contact (log)",
                  col = color$pinkmel(10), alpha = 0.5)

plot.scatter.dist(RNA.contacts.tre$dist,
                  log10(RNA.contacts.tre$contact),
                  cdc = RNA.contacts.cdc$r,xlim=c(-1000000,1000000), col = color$bluered(10),
                  ylim = c(-4.5, -0.5),
                  file="pdf/Fig2/FigS2c1.TRE-mRNA.1Mb.pdf",xlab="Position from TSS (kb)")

plot.scatter.dist(RNA.strand$dist,
                  log10(RNA.strand$contact),
                  xlim=c(-1000*2^(2*(5:0)),1000*2^(2*(0:5))),box=T, col = c("green", "orange", "purple"), alpha = 0.5,
                  ylim=c(-4.5,-0.5),file="pdf/Fig2/FigS2d.HiC.mRNA.pdf",xlab="tTRE-mRNA covariation")

plot.cdc.fold(RNA.contacts.tre$dist,
                  RNA.contacts.tre$contact,
                  xlim=c(-200000,200000), col = color$bluered(2), 
                  ylim = c(0, 0.02),
                  file="pdf/Fig2/FigS2d1.HiC_mRNA.1Mb.pdf",xlab="Position from TSS (kb)")

plot.cdc.fold(rna.cor$x,
              rna.cor$y,
              xlim=c(-200000,200000), col = color$bluered(2),
              ylim = c(-0.3, 0.6),
              file="pdf/Fig2/FigS2d2.cor_mRNA.1Mb.pdf",xlab="Position from TSS (kb)")

plot.cdc.fold(RNA.strand$dist,
              RNA.strand$contact,
              xlim=c(-200000,200000), col = color$bluered(2), 
              ylim = c(0, 0.02),
              file="pdf/Fig2/FigS2d3.HiC_mRNA_all.1Mb.pdf",xlab="Position from TSS (kb)")
