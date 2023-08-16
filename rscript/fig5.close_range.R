# Read strans specific normalzied reads
pl = read.table("readcount/Table.4d.procap.ambr.norm.pl.txt",header = 1)
mn = read.table("readcount/Table.4e.procap.ambr.norm.mn.txt",header = 1)

# Select enhancers only
enh.pos.all = read.table("window/procap.eTSSc.bed6",
                         col.names = c("chr", "pos", "pos2", "id", "score", "strand"))
pl.enh = pl %>%
  inner_join(enh.pos.all %>%
               select(chr, pos)) %>%
  distinct(chr, pos, .keep_all = T)
mn.enh = mn %>%
  inner_join(enh.pos.all %>%
               select(chr, pos)) %>%
  distinct(chr, pos, .keep_all = T)

# Distance measures between adjascent tTREs 
l = nrow(pl) - 1
dist = pl$pos[-1] - pl$pos[-(l+1)]
# Enhancers
l.enh = nrow(pl.enh) - 1
dist.enh = pl.enh$pos[-1] - pl.enh$pos[-(l.enh+1)]

# Correlation coefficients between adjascent tTREs
cor.div = sapply(1:l, function(i) cor(unlist(mn[i,-(1:2)]), unlist(pl[i+1,-(1:2)])))
cor.conv = sapply(1:l, function(i) cor(unlist(pl[i,-(1:2)]), unlist(mn[i+1,-(1:2)])))
cor.pl = sapply(1:l, function(i) cor(unlist(pl[i,-(1:2)]), unlist(pl[i+1,-(1:2)])))
cor.mn = sapply(1:l, function(i) cor(unlist(mn[i,-(1:2)]), unlist(mn[i+1,-(1:2)])))
# Enhancers
cor.e.div = sapply(1:l.enh, function(i) cor(unlist(mn.enh[i,-(1:2)]), unlist(pl.enh[i+1,-(1:2)])))
cor.e.conv = sapply(1:l.enh, function(i) cor(unlist(pl.enh[i,-(1:2)]), unlist(mn.enh[i+1,-(1:2)])))
cor.e.pl = sapply(1:l.enh, function(i) cor(unlist(pl.enh[i,-(1:2)]), unlist(pl.enh[i+1,-(1:2)])))
cor.e.mn = sapply(1:l.enh, function(i) cor(unlist(mn.enh[i,-(1:2)]), unlist(mn.enh[i+1,-(1:2)])))


d.co.e = dist.enh > 500 & dist.enh < 10000
cor.e.ss = c(cor.e.pl[d.co.e], cor.e.mn[d.co.e])
# Distance between skipping tTREs
l.s = nrow(pl) - 2
dist.s = pl$pos[-(1:2)] - pl$pos[-((l.s+1):(l.s+2))]

# Correlation coefficients between skipping tTREs
cor.div.s = sapply(1:l.s, function(i) cor(unlist(mn[i,-(1:2)]), unlist(pl[i+2,-(1:2)])))
cor.conv.s = sapply(1:l.s, function(i) cor(unlist(pl[i,-(1:2)]), unlist(mn[i+2,-(1:2)])))
cor.pl.s = sapply(1:l.s, function(i) cor(unlist(pl[i,-(1:2)]), unlist(pl[i+2,-(1:2)])))
cor.mn.s = sapply(1:l.s, function(i) cor(unlist(mn[i,-(1:2)]), unlist(mn[i+2,-(1:2)])))

# Distance categories
d.o = order(dist)
d.o.s = order(dist.s)
d.o = d.o[dist[d.o]>250]
d.o.s = d.o.s[dist.s[d.o.s]>250]

d.s = split(d.o, ceiling(seq_along(d.o)/1000))
d.s.s = split(d.o.s, ceiling(seq_along(d.o.s)/1000))
#Enhancers
d.o.e = order(dist.enh)
d.o.e = d.o.e[dist[d.o.e]>250]
d.e.s = split(d.o.e, ceiling(seq_along(d.o.e)/1000))

# quantiles
q = c(0.05, 0.25, 0.5, 0.75, 0.95)
cor.div.q = t(matrix(unlist(lapply(d.s, function(x) quantile(cor.div[x], q, na.rm=T))), nrow=5))
cor.conv.q = t(matrix(unlist(lapply(d.s, function(x) quantile(cor.conv[x], q, na.rm=T))), nrow=5))
cor.ss.q = t(matrix(unlist(lapply(d.s, function(x) quantile(c(cor.pl[x],cor.mn[x]), q, na.rm=T))), nrow=5))
dist.q = t(sapply(d.s, function(x) c(dist[x[1]], dist[x[length(x)]])))

cor.div.q.s = t(matrix(unlist(lapply(d.s.s, function(x) quantile(cor.div.s[x], q, na.rm=T))), nrow=5))
cor.conv.q.s = t(matrix(unlist(lapply(d.s.s, function(x) quantile(cor.conv.s[x], q, na.rm=T))), nrow=5))
cor.ss.q.s = t(matrix(unlist(lapply(d.s.s, function(x) quantile(c(cor.pl.s[x],cor.mn.s[x]), q, na.rm=T))), nrow=5))
dist.q.s = t(sapply(d.s.s, function(x) c(dist.s[x[1]], dist.s[x[length(x)]])))

plot.q = function(dist.q, cor, xdev=0) {
	lines( unlist(t(dist.q))+xdev, rep(cor, each=2) )
}


#pdf("FigS7a.pdf")
#plot(0,0, xlim=c(500,10000), ylim=c(0.1,0.4), xaxs='i', type='n')
#plot.q(dist.q, cor.conv.q[,3], -250)
#plot.q(dist.q, cor.div.q[,3], 250)
#plot.q(dist.q, cor.ss.q[,3])
#dev.off()

source('rscript/general/scatterplot.R')
d.1 = dist >= 250 & dist <= 1000
d.2 = dist > 5000 & dist < 10000
boxp.data = rbind(
  data.frame(Type = "0.25 - 1 kb", Orientation = "Sense", cor = c(cor.pl[d.1], cor.mn[d.1])),
  data.frame(Type = "> 5 kb", Orientation = "Sense", cor = c(cor.pl[d.2], cor.mn[d.2])),
  data.frame(Type = "0.25 - 1 kb", Orientation = "Convergent", cor = cor.conv[d.1]),
  data.frame(Type = "> 5 kb", Orientation = "Convergent", cor = cor.conv[d.2]),
  data.frame(Type = "0.25 - 1 kb", Orientation = "Divergent", cor = cor.div[d.1]),
  data.frame(Type = "> 5 kb", Orientation = "Divergent", cor = cor.div[d.2])) %>%
  mutate(Type = factor(Type, levels = c("0.25 - 1 kb", "> 5 kb")),
         Orientation = factor(Orientation, levels = c("Sense", "Convergent", "Divergent")))

pdf("pdf/Fig5/FigS5x.boxplot.pdf", height = 3, width = 2.5)
ggplot(boxp.data, aes(x = Orientation, y = cor, fill = Orientation)) +
  geom_boxplot(outlier.size = 0.2, color = "black") +
  facet_wrap(.~Type) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
        legend.position = "none",
        axis.line = element_line(linewidth = 0),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, color = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.text.y=element_text(color="black")) +
  scale_fill_manual(values = color$adjust(c('#7094CD','#80C080','#CC72AE'), 1 , 0.8)) +
  ylab("Correlation coefficients") +
  ggtitle("Distance") +
  coord_cartesian(ylim = c(-0.75, 1.15))

dev.off()

custom.scatter.open(file="pdf/Fig5/FigS5a.sense_strand.pdf",xlab="Distance (bp)",ylab="Correlation coefficients", main = "Sense strand", 
                    height=3,width=4)
d.co = dist > 500 & dist < 10000
cor.ss = c(cor.pl[d.co],cor.mn[d.co])
custom.scatter.plot(rep(dist[d.co],2), cor.ss, xlim=c(500,10000), color=color$bluered(10), pval=F, height=3, width=4, lwd=1.5, border.width = 1.5)
loess.ss = loess.smooth(rep(dist[d.co],2), cor.ss, span=0.2, eval=200)
lines(loess.ss$x, loess.ss$y, lwd=2*100/37, col='white')
lines(loess.ss$x, loess.ss$y, lwd=100/37, col='black')
custom.scatter.label(xlab="Distance (bp)", ylab="Correlation coefficients", main = "Sense strand")
custom.scatter.close()
sum(d.co)
custom.scatter.open(file="pdf/Fig5/FigS5b.convergent.pdf",xlab="Distance (bp)",ylab="Correlation coefficients", main = "Convergent", height=3,width=4)
custom.scatter.plot(dist[d.co], cor.conv[d.co], xlim=c(500,10000), color=color$bluered(10), pval=F, height=3, width=4, lwd=1.5, border.width = 1.5)
loess.conv = loess.smooth(dist[d.co], cor.conv[d.co], span=0.2, eval=200)
lines(loess.conv$x, loess.conv$y, lwd=2*100/37/3*4, col='white')
lines(loess.conv$x, loess.conv$y, lwd=100/37/3*4, col='black')
custom.scatter.label(xlab="Distance (bp)", ylab="Correlation coefficients", main = "Convergent")
custom.scatter.close()
sum(d.co.e)

custom.scatter.open(file="pdf/Fig5/FigS5c.divergent.pdf",xlab="Distance (bp)",ylab="Correlation coefficients", main = "Divergent", height=3,width=4)
custom.scatter.plot(dist[d.co], cor.div[d.co], xlim=c(500,10000), color=color$bluered(10), pval=F, height=3, width=4, lwd=1.5, border.width = 1.5)
loess.div = loess.smooth(dist[d.co], cor.div[d.co], span=0.2, eval=200)
lines(loess.div$x, loess.div$y, lwd=2*100/37/3*4, col='white')
lines(loess.div$x, loess.div$y, lwd=100/37/3*4, col='black')
custom.scatter.label(xlab="Distance (bp)", ylab="Correlation coefficients", main = "Divergent")
custom.scatter.close()

pdf("pdf/Fig5/FigS5b.orientation_compared.pdf",height=1.4, width=2)
par(mar=c(1,1,0.5,0.5), mgp=c(0.5,0.2,0))
plot(0,xlim=c(500,5000), ylim=c(0.1,0.4), axes=F, xlab='', ylab='', xaxs='i')
lines(loess.ss$x, loess.ss$y, lwd=1.5, col='#7094CD')
lines(loess.conv$x, loess.conv$y, lwd=1.5, col='#80C080')
lines(loess.div$x, loess.div$y, lwd=1.5, col='#CC72AE')
legend("topright",legend=c("Sense","Convergent","Divergent"),col=c('#7094CD','#80C080','#CC72AE'),lwd=1.5, bty='n',cex=7/12)
axis(1,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5, mgp=c(0.5,-0.1,0))
axis(2,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5)
box()
dev.off()

pdf("pdf/Fig5/Fig5b.orientation_distance_shifted.pdf",height=1.25, width=2)
par(mar=c(1,1,0.5,0.5), mgp=c(0.5,0.2,0))
plot(0,xlim=c(500,5000), ylim=c(0.1,0.35), axes=F, xlab='', ylab='', xaxs='i')
lines(loess.ss$x-250, loess.ss$y, lwd=1.5, col='#7094CD')
lines(loess.conv$x-500, loess.conv$y, lwd=1.5, col='#80C080')
lines(loess.div$x, loess.div$y, lwd=1.5, col='#CC72AE')
legend("topright",legend=c("Sense","Convergent","Divergent"),col=c('#7094CD','#80C080','#CC72AE'),lwd=1.5, bty='n',cex=7/12)
axis(1,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5, mgp=c(0.5,-0.1,0))
axis(2,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5)
box(lwd = 1.25)
dev.off()

# Enhancer comparison only
loess.e.ss = loess.smooth(rep(dist.enh[d.co.e], 2), cor.e.ss, span=0.2, eval=200)
loess.e.conv = loess.smooth(dist.enh[d.co.e], cor.e.conv[d.co.e], span=0.2, eval=200)
loess.e.div = loess.smooth(dist.enh[d.co.e], cor.e.div[d.co.e], span=0.2, eval=200)

pdf("pdf/Fig5/Fig5b2.enhancer_distance_shifted.pdf",height=1.4, width=2)
par(mar=c(1,1,0.5,0.5), mgp=c(0.5,0.2,0))
plot(0,xlim=c(500,5000), ylim=c(0.1,0.3), axes=F, xlab='', ylab='', xaxs='i')
lines(loess.e.ss$x-250, loess.e.ss$y, lwd=1.5, col='#7094CD')
lines(loess.e.conv$x-500, loess.e.conv$y, lwd=1.5, col='#80C080')
lines(loess.e.div$x, loess.e.div$y, lwd=1.5, col='#CC72AE')
legend("topright",legend=c("Sense","Convergent","Divergent"),col=c('#7094CD','#80C080','#CC72AE'),lwd=1.5, bty='n',cex=7/12)
axis(1,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5, mgp=c(0.5,-0.1,0))
axis(2,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5)
box()
dev.off()

pdf("pdf/Fig5/Fig5b3.enhancer_orientation.pdf",height=1.4, width=2)
par(mar=c(1,1,0.5,0.5), mgp=c(0.5,0.2,0))
plot(0,xlim=c(250,5000), ylim=c(0.1,0.4), axes=F, xlab='', ylab='', xaxs='i')
lines(loess.e.ss$x, loess.e.ss$y, lwd=1.5, col='#7094CD')
lines(loess.e.conv$x, loess.e.conv$y, lwd=1.5, col='#80C080')
lines(loess.e.div$x, loess.e.div$y, lwd=1.5, col='#CC72AE')
legend("topright",legend=c("Sense","Convergent","Divergent"),col=c('#7094CD','#80C080','#CC72AE'),lwd=1.5, bty='n',cex=7/12)
axis(1,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5, mgp=c(0.5,-0.1,0))
axis(2,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5)
box()
dev.off()

# Loess for skipping
d.co.s = dist.s > 500 & dist.s < 10000
cor.ss.s = c(cor.pl.s[d.co.s],cor.mn.s[d.co.s])
loess.ss.s = loess.smooth(rep(dist.s[d.co.s],2), cor.ss.s, span=0.2, eval=200)
loess.conv.s = loess.smooth(dist.s[d.co.s], cor.conv.s[d.co.s], span=0.2, eval=200)
loess.div.s = loess.smooth(dist.s[d.co.s], cor.div.s[d.co.s], span=0.2, eval=200)
sum(d.co.s)

pdf("pdf/Fig5/Fig5d1.sense.skipped.pdf",height=1.4, width=2)
par(mar=c(1,1,0.5,0.5), mgp=c(0.5,0.2,0))
plot(0,xlim=c(500,5000), ylim=c(0.1,0.4), axes=F, xlab='', ylab='', xaxs='i')
lines(loess.ss$x-250, loess.ss$y, lwd=1.5, col='#7094CD')
lines(loess.ss.s$x-250, loess.ss.s$y, lwd=1.5, col='#7094CD', lty=2)
legend("topright",legend=c("Neighboring","Interlacing"),col=c('#7094CD','#7094CD'),lwd=1.5, bty='n',cex=7/12, lty=c(1,2))
axis(1,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5, mgp=c(0.5,-0.1,0))
axis(2,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5)
box()
dev.off()

pdf("pdf/Fig5/Fig5d2.convergent_skipped.pdf",height=1.4, width=1.8)
par(mar=c(1,1,0.5,0.5), mgp=c(0.5,0.2,0))
plot(0,xlim=c(500,5000), ylim=c(0.1,0.4), axes=F, xlab='', ylab='', xaxs='i')
lines(loess.conv$x-500, loess.conv$y, lwd=1.5, col='#80C080')
lines(loess.conv.s$x-500, loess.conv.s$y, lwd=1.5, col='#80C080', lty=2)
legend("topright",legend=c("Adjacent","Interleaved"),col=c('#80C080','#80C080'),lwd=1.5, bty='n',cex=7/12, lty=c(1,2))
axis(1,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5, mgp=c(0.5,-0.1,0))
axis(2,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5)
box(lwd = 1.25)
dev.off()

pdf("pdf/Fig5/Fig5d2.divergent_skipped.pdf",height=1.4, width=2)
par(mar=c(1,1,0.5,0.5), mgp=c(0.5,0.2,0))
plot(0,xlim=c(500,5000), ylim=c(0.1,0.4), axes=F, xlab='', ylab='', xaxs='i')
lines(loess.div$x, loess.div$y, lwd=1.5, col='#CC72AE')
lines(loess.div.s$x, loess.div.s$y, lwd=1.5, col='#CC72AE', lty=2)
legend("topright",legend=c("Neighboring","Interlacing"),col=c('#CC72AE','#CC72AE'),lwd=1.5, bty='n',cex=7/12, lty=c(1,2))
axis(1,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5, mgp=c(0.5,-0.1,0))
axis(2,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5)
box()
dev.off()

######
# HiC analysis of interlaced contacts
#

# Read tTREs in hg38 coordinates again
# Round up to 1 kb resolution
pro.1kb.hg38 = read.table("HiC/tTRE/hg38.bed", header = F, col.names = c("chr", "pos", "pos2", "id")) %>%
  mutate(pos2 = floor(pos/1000)*1000)

# Extract Hi-C contract matrix +/- 10 kb region in 1kb resolution around the tTRES
# Generate a bash script using Cooler in the syntax example below
# cooler dump -b -t pixels -f --join -r chr3:20000k-20005k -r2 chr3:19000k-21000k HiC/4DNFIXP4QG5B.mcool::/resolutions/5000 > tmp

bash.command.1kb = sapply(1:nrow(pro.1kb.hg38),
                      function(i) {
                        chr = pro.1kb.hg38$chr[i]
                        pos = pro.1kb.hg38$pos2[i]
                        return(paste0("cooler dump -b -t pixels -f --join -r ",
                                      chr, ":",
                                      pos/1000, "k-",
                                      pos/1000 + 1, "k -r2 ",
                                      chr, ":",
                                      pos/1000 - 10, "k-",
                                      pos/1000 + 10, "k HiC/4DN/4DNFIXP4QG5B.mcool::/resolutions/1000 >> HiC/4DN/tTRE.1kb.txt"))})
writeLines(bash.command.1kb, "HiC/get.1kbContacts.sh")

# Read tTREcontacts.txt file

HiC.1kb = read.table("HiC/4DN/tTRE.1kb.txt", header = F,
                          col.names = c("chr1", "pos1", "pos1e", "chr2", "pos2", "pos2e",
                                        "count", "contact"), fill = T)
HiC.1kb[is.na(HiC.1kb)] = 0

HiC.1kb = HiC.1kb %>%
  distinct(chr1, pos1, pos2, .keep_all = T) %>%
  mutate(chr = chr1) %>%
  select(chr, pos1, pos2, count, contact)

# list of adjacent and skipped tTRE pairs in the hg38 list
closepair.next = pro.1kb.hg38 %>%
  mutate(pos1 = pos, id1 = id) %>%
  select(chr, pos1, id1) %>%
  filter(id1 != tail(id1, 1)) %>%
  bind_cols(pro.1kb.hg38 %>%
              mutate(pos2 = pos, id2 = id) %>%
              select(pos2, id2) %>%
              filter(id2 != head(id2, 1))) %>%
  mutate(dist = abs(pos2-pos1)) %>%
  filter(dist < 5000) %>%
  mutate(pos1 = floor(pos1/1000)*1000,
         pos2 = floor(pos2/1000)*1000)
closepair.skip = pro.1kb.hg38 %>%
  mutate(pos1 = pos, id1 = id) %>%
  select(chr, pos1, id1) %>%
  filter(!(id1 %in% tail(id1, 2))) %>%
  bind_cols(pro.1kb.hg38 %>%
              mutate(pos2 = pos, id2 = id) %>%
              select(pos2, id2) %>%
              filter(!(id2 %in% head(id2, 2)))) %>%
  mutate(dist = abs(pos2-pos1)) %>%
  filter(dist < 5000) %>%
  mutate(pos1 = floor(pos1/1000)*1000,
        pos2 = floor(pos2/1000)*1000)

HiC.close.next = closepair.next %>%
  left_join(HiC.1kb, by = c("chr", "pos1", "pos2")) %>%
  select(chr, id1, id2, dist, contact) %>%
  drop_na()

HiC.close.skip = closepair.skip %>%
  left_join(HiC.1kb, by = c("chr", "pos1", "pos2")) %>%
  select(chr, id1, id2, dist, contact) %>%
  drop_na()

# Compared the loess fits of contact by distance between next-pairs and skip-pairs
# Use the same plotting format

loess.hiC = loess.smooth(HiC.close.next$dist, HiC.close.next$contact, span=0.2, eval=200)
loess.hiC.s = loess.smooth(HiC.close.skip$dist, HiC.close.skip$contact, span=0.2, eval=200)

pdf("pdf/Fig5/FigS5.HiC_skipped.pdf",height=1.4, width=1.8)
par(mar=c(1,1,0.5,0.5), mgp=c(0.5,0.2,0))
plot(c(-Inf, -Inf),xlim=c(500,5000), axes=F, ylim = c(0.00, 0.03), xlab='', ylab='', xaxs='i')
lines(loess.hiC$x, loess.hiC$y, lwd=1.5, col='#808080')
lines(loess.hiC.s$x, loess.hiC.s$y, lwd=1.5, col='#808080', lty=2)
legend("topright",legend=c("Adjacent","Interleaved"),col=c('#808080','#808080'),lwd=1.5, bty='n',cex=7/12, lty=c(1,2))
axis(1,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5, mgp=c(0.5,-0.1,0))
axis(2,lwd=0, lwd.tick=1, tck=-0.03, cex.axis=0.5)
box(lwd = 1.25)
dev.off()


# Statistical testing of convergent vs divergent
conv.pv.test = cor.conv[dist >= 1000 & dist <= 2000]
div.pv.test = cor.div[dist >= 500 & dist <= 1500]
sense.pv.test = cor.ss[dist >= 750 & dist <= 1750]
t.test(conv.pv.test, div.pv.test)
t.test(conv.pv.test,sense.pv.test)
t.test(div.pv.test, sense.pv.test)


# Statisticall testing of convergent vs HiC at 500-1500 whether interlaced are lower
conv.pval.table = list(cor.conv[dist >= 500 & dist <= 1500],
                       cor.conv.s[dist.s >= 500 & dist.s <= 1500])
HiC.pval.table = list(HiC.close.next %>%
                        filter(dist >= 500 & dist <= 1500) %>%
                        select(contact) %>%
                        unlist,
                      HiC.close.skip %>%
                        filter(dist >= 500 & dist <= 1500) %>%
                        select(contact) %>%
                        unlist)

t.test(conv.pval.table[[1]], conv.pval.table[[2]]) # Correlation: p < 2.2 x 10^-16
t.test(HiC.pval.table[[1]], HiC.pval.table[[2]]) # Hi-C: p = 0.078



