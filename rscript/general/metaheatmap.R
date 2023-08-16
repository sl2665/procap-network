get_meta_reads=function(bed,bgPl,bgMn=NULL,window=2000)
{
	script.path="~/Sandbox/procap/rscript/"
	system(paste("awk '$6==",'"+"{a=$2-',window,';if(a<1) a=1; print $1"\\t"a"\\t"$2+',window,';next}{b=$3-',window,'-1;if(b<1) b=1; print $1"\\t"b"\\t"$3+',window-1,"}' ",bed," > bed.tmp",sep=""))
	system(paste(script.path,"bedreads.sh ",bgPl," bed.tmp > reads.tmp",sep=""))
}
