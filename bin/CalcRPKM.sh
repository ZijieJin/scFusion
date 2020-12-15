bed=$1
dir=$2
outdir=$3
mystart=$4
myend=$5
for((i=${mystart};i<=${myend};i++))
do
	file=`ls ${dir}/${i}/*Aligned.sortedByCoord.out.bam`
	samtools index ${file}
	export total_reads=`samtools idxstats ${file}|awk -F '\t' '{s+=$3}END{print s}'`
	bedtools multicov -bams ${file} -bed $bed |perl -alne '{$len=$F[2]-$F[1];if($len <1 ){print "$F[3]\t$F[4]\t0" }else{$rpkm=(1000000000*$F[4]/($len* $ENV{total_reads}));print "$F[3]\t$F[4]\t$rpkm"}}' >  ${outdir}/${i}.rpkm.txt
done