do=T

# Read data files
if(do) {
	pro=read.table('readcount/procap.ambr.norm.txt')
	rna=read.table('readcount/rnaseq.norm.txt')
	rna.list=read.table('readcount/rnaseq.list.txt',stringsAsFactors=F)[,1]
	pro.list=read.table('readcount/procap.ambr.list.txt',stringsAsFactors=F)[,1]
}
# Function to get the distribution of the difference betweeen replicates
get.rep.dev=function(mat,list,suffix="r"){ 
	rep.list=intersect(list,paste(list,suffix,sep="")) 
	rep1.id=match(substr(rep.list,1,nchar(rep.list)-nchar(suffix)),list) 
	rep2.id=match(rep.list,list) 
	n=length(rep.list) 
	get.diff=function(i) {
		return(abs(mat[,rep1.id[i]]-mat[,rep2.id[i]]))
	}
	s=sapply(1:n,get.diff)
	return(s)
}

# Function to get pairwise distribution of the difference betweeen samples
get.dev=function(mat,list,prefixlen=7){ 
	l=substring(list,1,prefixlen)
	p=c()
	for(i in 1:(length(l)-1))
		for(j in (i+1):length(l))
			if(l[i]!=l[j]) p=rbind(p,c(i,j))
	n=nrow(p)
	get.diff=function(i){
		return(abs(mat[,p[i,1]]-mat[,p[i,2]]))
	}
	s=sapply(1:n,get.diff)
	return(s)
}

if(do) {
	pro.rep.dev=get.rep.dev(pro[,-(1:2)],pro.list)
	pro.dev=get.dev(pro[,-(1:2)],pro.list)
	rna.rep.dev=get.rep.dev(rna[,-(1:2)],rna.list,suffix="r1")
	rna.dev=get.dev(rna[,-(1:2)],rna.list)

	n=nrow(pro)
	pro.t=sapply(1:n,function(i) {return(t.test(pro.rep.dev[i,],pro.dev[i,])$p.value)})
	pro.fdr=p.adjust(pro.t,method="fdr")

	n=nrow(rna)
	rna.t=sapply(1:n,function(i) {return(t.test(rna.rep.dev[i,],rna.dev[i,])$p.value)})
	rna.exp=apply(rna[,-(1:2)],1,mean)>0
	rna.fdr=p.adjust(rna.t,method="fdr")
}

pro.h=hist(pro.t,breaks=seq(0,1,by=0.05),plot=F)
pdf("pdf/Fig1S/expression variation_p_value_histogram.pdf",width=4,height=3)
par(mar=c(4,4,1,1),mgp=c(2,0.5,0))
plot(pro.h,freq=F,main="",xlab="p-value",ylim=c(0,6),col='#c0c0c0')
pro.null=mean(pro.h$density[8:20])
lines(c(0,1),rep(pro.null,2),lty=2,lwd=1.5,col='firebrick')
lines(c(0,1),c(1,1),lty=2,lwd=1.5)
dev.off()

# write.table(pro[pro.fdr<0.25,1:2],'~/Sandbox/procap/fig/table/Table.5b.procap.nTSS.var.txt',sep="\t",quote=F,col.names=F,row.names=F)


pro.cor.h=hist(unlist(pro.cor$scd$pval),breaks=seq(0,1,by=0.05),plot=F)
pdf("pdf/Fig1S/Coexpression_p_value_histogram.pdf",width=4,height=3)
par(mar=c(4,4,1,1),mgp=c(2,0.5,0))
plot(pro.cor.h,freq=F,main="",xlab="p-value",ylim=c(0,6),col='#c0c0c0')
pro.null=mean(pro.cor.h$density[8:20])
lines(c(0,1),rep(pro.null,2),lty=2,lwd=1.5,col='firebrick')
lines(c(0,1),c(1,1),lty=2,lwd=1.5)
dev.off()

pro.cor.fdr = p.adjust(unlist(pro.cor$scd$pval), method = "fdr") * sign(unlist(pro.cor$scd$corr))
sum(pro.cor.fdr < 0.0006 & pro.cor.fdr > 0)/length(pro.cor.fdr)
