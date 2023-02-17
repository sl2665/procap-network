source('rscript/functions.R')

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


# Make interactions
pro.window = make.window.list(pro.near)
pro.cor  = cor.pro(pro.window, pro, pro)

# Plot all interactions
plot.scatter.dist(pro.cor$x,pro.cor$y,xlim=c(0,1000*2^(0:10)),box=T,ylim=c(-0.7,1),file="Fig6C.pdf")
plot.scatter.dist(pro.cor$x,pro.cor$y,cdc = pro.cor$cdc, xlim=c(0,1000000),file="FigS6B.pdf")

# Interactions with genotype correlation less than 0.05 or greater than 0.5
geno.cor = cor.pro(pro.window, geno, geno)
plot.scatter.dist(pro.cor$x[abs(geno.cor$y)<0.05],pro.cor$y[abs(geno.cor$y)<0.05],xlim=c(0,1000*2^(0:10)),box=T,ylim=c(-0.7,1),file="FigS6C.pdf")

# Interactions between promoter and enhancer
pro.enh.window = make.window.asym.list(pro.enh.near)
pro.enh.cor  = cor.pro(pro.enh.window, pro.pro, pro.enh, asym=T)
plot.scatter.dist(pro.enh.cor$x,pro.enh.cor$y,xlim=c(0,1000*2^(0:10)),box=T,ylim=c(-0.7,1),file="FigS6D.pdf")

# Interactions between mRNA and TRE
rna.window = make.rna.window.list(rna.near)
rna.cor = cor.rna_pro(rna.window, pro)
plot.scatter.dist(rna.cor$x,rna.cor$y,xlim=c(-1000*2^(2*(5:0)),1000*2^(2*(0:5))),box=T,ylim=c(-0.5,0.75),file="Fig6D.pdf",xlab="Position from TSS (kb)")
plot.scatter.dist(rna.cor$x,rna.cor$y,rna.cor$cdc,xlim=c(-1000000,1000000),file="FigS6E.pdf",xlab="Position from TSS (kb)")

# Interactions between mRNA and promoter or enhancer
rna.pro.window = make.rna.window.list(rna.pro.near)
rna.pro.cor = cor.rna_pro(rna.pro.window, pro.pro)
plot.scatter.dist(rna.pro.cor$x,rna.pro.cor$y,xlim=c(-1000*2^(2*(5:0)),1000*2^(2*(0:5))),box=T,ylim=c(-0.5,0.75),file="FigS6F.pdf",xlab="Position from TSS (kb)")

rna.enh.window = make.rna.window.list(rna.enh.near)
rna.enh.cor = cor.rna_pro(rna.enh.window, pro.enh)
plot.scatter.dist(rna.enh.cor$x,rna.enh.cor$y,xlim=c(-1000*2^(2*(5:0)),1000*2^(2*(0:5))),box=T,ylim=c(-0.5,0.75),file="FigS6G.pdf",xlab="Position from TSS (kb)")

