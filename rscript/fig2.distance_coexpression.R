source('rscript/functions.R')
source('rscript/general/colors.R')

# Read PRO-cap data
pro = read.table('readcount/procap.var.norm.txt',stringsAsFactors=F)
pro.near = read.table('window/procap.window.all-all.1M.txt',stringsAsFactors=F)
geno = read.table('genotype/genotype.var.txt', stringsAsFactors=F)
pro.pro = read.table('readcount/procap.var.promoter.norm.txt', stringsAsFactors=F)
pro.enh = read.table('readcount/procap.var.enhancer.norm.txt', stringsAsFactors=F)
pro.enh.near = read.table('window/procap.window.pro-enh.1M.txt', stringsAsFactors=F)

rna=read.table('readcount/rnaseq.expr.norm.txt',stringsAsFactors=F)
rna.near=read.table('window/rnaseq.window.gene-all.1M.txt', stringsAsFactors=F)
rna.pro.near=read.table('window/rnaseq.window.gene-pro.1M.txt', stringsAsFactors=F)
rna.enh.near=read.table('window/rnaseq.window.gene-enh.1M.txt', stringsAsFactors=F)
rna.list=read.table('readcount/rnaseq.list.txt',stringsAsFactors=F)[,1]
pro.list=read.table('readcount/procap.ambr.list.txt',stringsAsFactors=F)[,1]
rna.pos=read.table('window/rnaseq.pos.bed', stringsAsFactors=F)
rna[,2] = rna.pos[,2]

# Make interactions (between TREs) 
pro.window = make.window.list(pro.near) # pair list
pro.cor  = cor.pro(pro.window, pro, pro) # cor.pro : function to calculate correlation between tTRE pair lists

# Plot all interactions
plot.scatter.dist(pro.cor$x,pro.cor$y,xlim=c(0,1000*2^(0:10)),box=T,ylim=c(-0.7,1),xlab="All tTRE covariation", file="pdf/Fig2/Fig2a.pdf",
                  col = color$jet(10), alpha = 0.1)
plot.scatter.dist(pro.cor$x,pro.cor$y,cdc = pro.cor$cdc, xlim=c(0,1000000),file="pdf/Fig2/Fig2a1.all_tTRE.1Mb.pdf",
                  xlab = "Distance (kb)", col = color$bluered(10))
plot.scatter.dist(pro.cor$x,pro.cor$y,cdc = pro.cor$cdc, xlim=c(0,500000),file="pdf/Fig2/Fig2a1.all_tTRE.500kb.pdf",
                  xlab = "Distance (kb)", col = color$bluered(10))
plot.scatter.dist(pro.cor$x,pro.cor$y,cdc = pro.cor$cdc, xlim=c(0,200000),file="pdf/Fig2/Fig2a1.all_tTRE.200kb.pdf",
                  xlab = "Distance (kb)", col = color$bluered(10))
plot.fdr.dist(pro.cor, file = "pdf/Fig2/Fig2a2.all_tTRE.sig_percent.pdf", main = "All tTRE",
              col = color$jet(10), alpha = 0.5)


# Interactions with genotype correlation less than 0.05 which are not genetically associated tTREs
geno.cor = cor.pro(pro.window, geno, geno)
plot.scatter.dist(pro.cor$x[abs(geno.cor$y)<0.05],pro.cor$y[abs(geno.cor$y)<0.05],xlim=c(0,1000*2^(0:10)),box=T,ylim=c(-0.7,1),
                  file="pdf/Fig2f.indep_geno.pdf", col = color$grpa(10), alpha = 0.5)

# Percent of correlations that are highly associated genetically
geno.pval = unlist(geno.cor$scd$pval)
geno.fdr = p.adjust(geno.pval, method = "fdr")
sum(geno.fdr < 0.1) / length(geno.fdr) # About 5% of the pairs are significanctly associated genetically
sum(geno.cor$y < 0.05) / length(geno.cor$y) # About 70% of the not highly associated, correlation of genotypes less than 0.05

# Interactions between promoter and enhancer # of enhancers 4006 variable promoters and 25688 variable enhancers
pro.enh.window = make.window.asym.list(pro.enh.near)
pro.enh.cor  = cor.pro(pro.enh.window, pro.pro, pro.enh, asym=T)
plot.scatter.dist(pro.enh.cor$x,pro.enh.cor$y,xlim=c(0,1000*2^(0:10)),box=T,ylim=c(-0.7,1),xlab="Enhancer-promoter covariation", file="pdf/Fig2/Fig2b.enh-prm.pdf",
                  col = color$pinkmel(10), alpha = 0.1)
plot.scatter.dist(pro.enh.cor$x,pro.enh.cor$y,cdc = pro.enh.cor$cdc, xlim=c(0,200000),file="pdf/Fig2/Fig2b1.enh-prm.200kb.pdf", xlab = "Distance (kb)", col = color$bluered(10))
plot.scatter.dist(pro.enh.cor$x,pro.enh.cor$y,cdc = pro.enh.cor$cdc, xlim=c(0,500000),file="pdf/Fig2/Fig2b2.enh-prm.500kb.pdf", xlab = "Distance (kb)", col = color$bluered(10))
plot.scatter.dist(pro.enh.cor$x,pro.enh.cor$y,cdc = pro.enh.cor$cdc, xlim=c(0,1000000),file="pdf/Fig2/Fig2b3.enh-prm.1Mb.pdf", xlab = "Distance (kb)", col = color$bluered(10))
plot.fdr.dist(pro.enh.cor, file = "pdf/Fig2/Fig2b2.enh-prm.pdf", main = "Enhancer-promoter pairs",  col = color$pinkmel(10), alpha = 0.1)

# Interactions between mRNA and TRE
rna.window = make.rna.window.list(rna.near)
rna.cor = cor.rna_pro(rna.window, pro)
plot.scatter.dist(rna.cor$x,rna.cor$y,xlim=c(-1000*2^(2*(5:0)),1000*2^(2*(0:5))),box=T,ylim=c(-0.5,1),file="pdf/Fig2/Fig2c.TRE-mRNA.pdf",xlab="tTRE-mRNA covariation",
                  col = color$japa(10), alpha = 0.5)
plot.scatter.dist(rna.cor$x,rna.cor$y,rna.cor$cdc,xlim=c(-200000,200000),file="pdf/Fig2/Fig2c1.TRE-mRNA.200kb.pdf",xlab="Position from TSS (kb)", col = color$bluered(10))
plot.scatter.dist(rna.cor$x,rna.cor$y,rna.cor$cdc,xlim=c(-500000,500000),file="pdf/Fig2/Fig2c2.TRE-mRNA.500kb.pdf",xlab="Position from TSS (kb)", col = color$bluered(10))
plot.scatter.dist(rna.cor$x,rna.cor$y,rna.cor$cdc,xlim=c(-1000000,1000000),file="pdf/Fig2/Fig2c3.TRE-mRNA.1Mb.pdf",xlab="Position from TSS (kb)", col = color$bluered(10))

# Interactions between mRNA and promoter
rna.pro.window = make.rna.window.list(rna.pro.near)
rna.pro.cor = cor.rna_pro(rna.pro.window, pro.pro)
plot.scatter.dist(rna.pro.cor$x,rna.pro.cor$y,xlim=c(-1000*2^(2*(5:0)),1000*2^(2*(0:5))),box=T,ylim=c(-0.5,0.75),file="pdf/Fig2/Fig2d.prm-mRNA.pdf",xlab="Promoter-mRNA covariation",
                  col = color$pinkpa(10), alpha = 0.5)
plot.scatter.dist(rna.pro.cor$x,rna.pro.cor$y,rna.pro.cor$cdc,xlim=c(-1000000,1000000),file="pdf/Fig2/Fig2d1.prm-mRNA.1Mb.pdf",xlab="Position from TSS (kb)", col = color$bluered(10))

# Interaction between mRNA and enhancer
rna.enh.window = make.rna.window.list(rna.enh.near)
rna.enh.cor = cor.rna_pro(rna.enh.window, pro.enh)
plot.scatter.dist(rna.enh.cor$x,rna.enh.cor$y,xlim=c(-1000*2^(2*(5:0)),1000*2^(2*(0:5))),box=T,ylim=c(-0.5,0.75),file="pdf/Fig2/Fig2e.enh-mRNA.pdf",
                  xlab="Enhancer-mRNA covariation", col = color$pinkmel(10), alpha = 0.5)
plot.scatter.dist(rna.enh.cor$x,rna.enh.cor$y,rna.enh.cor$cdc,xlim=c(-1000000,1000000),file="pdf/Fig2/Fig2e1.enh-mRNA.1Mb.pdf",xlab="Position from TSS (kb)", col = color$bluered(10))

# Comparison of RNA-enhancer correlation between upstream and downstream
rna.updn.pval = function(start = 4, end = 8) {
  rna.enh.cor.up = rna.enh.cor$y[rna.enh.cor$x > -end*1000 & rna.enh.cor$x < -start*1000]
  rna.enh.cor.dn = rna.enh.cor$y[rna.enh.cor$x > start*1000 & rna.enh.cor$x < end*1000]
  return(t.test(rna.enh.cor.up, rna.enh.cor.dn))  
}
# pvalue of 1-4 kb up and down
rna.updn.pval(1,4) # downstream is higher, pval = 0.022
# pvalue of 1-2 kb up and down
rna.updn.pval(1,2) # downstream is higher, pval = 0.0001
# pvalue of 2-4 kb up and down
rna.updn.pval(2,4) # downstream/upstream almost same, pval = 0.83
# pvalue of 4-8 kb up and down
rna.updn.pval(4,8) # downstream is lower, pval = 0.003
# pvalue of 4-16 kb up and down
rna.updn.pval(4,16) # downstream is lower, pval = 0.0000003
# pvalue of 16-32 kb up and down
rna.updn.pval(16,32) # downstream is lower, pval = 0.0000000009
# pvalue of 16-64 kb up and down
rna.updn.pval(16,64) # downstream is lower, pval = 0

# p value and fdr of all correlations
pro.cor.fdr = p.adjust(unlist(pro.cor$scd$pval), method = "fdr")
pro.cor.fdr %>% head


pdf("pdf/test.pdf")
dddd = abs(as.numeric(pro.bed$pos[-1]) - as.numeric(pro.bed$pos[-nrow(pro.bed)]))
hist(dddd, breaks = c(-Inf, seq(0, 200000, by = 5000), Inf),
     xlim = c(0, 200000), freq = T)
dev.off()

length(dddd)
sum(dddd > 300)/length(dddd)

# P-value compared to interchromosomal distribution

xlim=c(0,1000*2^(0:10))
ic.pval = function(xlim) t.test(pro.enh.cor$y[pro.enh.cor$x > xlim[1] & pro.enh.cor$x <= xlim[2]], interchromosomal)$p.val
ic.pval(c(xlim[1:2]))
xlim[5:6]/1000
