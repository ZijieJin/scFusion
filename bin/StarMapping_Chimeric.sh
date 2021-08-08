#!/bin/bash
filedir=$1
mystart=$2
myend=$3
outdir=$4
genomedir=$5
ncore=$6
for ((i=${mystart};i<=${myend};i++))
do
	if [ -f ${filedir}/${i}_1.fastq ];then
		mkdir -p ${outdir}/${i}
		STAR --runThreadN ${ncore} --genomeDir ${genomedir} --readFilesIn ${filedir}/${i}_1.fastq ${filedir}/${i}_2.fastq --outSAMtype BAM SortedByCoordinate --chimOutType SeparateSAMold --outSAMunmapped Within KeepPairs --quantMode GeneCounts --outFileNamePrefix ${outdir}/${i}/human --chimSegmentMin 12 --chimJunctionOverhangMin 8 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5 --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimNonchimScoreDropMin 10 --peOverlapMMp 0.1
	fi
done
