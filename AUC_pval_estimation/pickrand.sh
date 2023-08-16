awk -v ns=$1 '{
	chr[NR] = $1;
	start[NR] = $2;
	end[NR] = $3;
	relpos[NR] = cumpos;
	cumpos += $3 - $2;
	++n;
	for(i=int(relpos[NR]/100);i<=int(cumpos/100);++i) bpos[i]=NR;
}
END {
	srand();
	for(i=1;i<=100;++i) {
		fn = sprintf("shuf.%02dk.%03d.txt",ns,i)
		for(j=0;j<ns*1000;++j) {
			a = int(cumpos*rand());
			l = bpos[int(a/100)];
			print chr[l]"\t"start[l]+a-relpos[l] > fn;
		}
	}
}' ../Tfbs.merge.bed

