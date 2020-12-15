#!/bin/bash

FilePath=$1
mystart=$2
myend=$3
prefix=$4
weightfile=$5
hg19file=$6
gtf=$7
codedir=$8

if [ "${prefix}" = "." ]
then
	prefix=""
fi

mkdir -p ${FilePath}/ChiDist/

python ${codedir}/MyPredict.py  ${FilePath}/ChiDist/${prefix}Prob.txt ${weightfile} ${prefix}
paste ${FilePath}/ChiDist/${prefix}ChiDist_middle.txt ${FilePath}/ChiDist/${prefix}Prob.txt > ${FilePath}/ChiDist/${prefix}ChiDist.txt
python ${codedir}/FilterChiDist.py ${FilePath}/ChiDist/${prefix}ChiDist.txt > ${FilePath}/ChiDist/${prefix}ChiDist_filtered.txt