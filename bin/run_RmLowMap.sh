#!/bin/bash
start=$1
end=$2
filedir=$3
outdir=$4
thres=$5
mappabilityfile=$6
for ((i=${start};i<=${end};i++))
do
	file=`ls ${filedir}/${i}/*Chimeric.out.sam`
	if [[ -n ${file} ]]; then
		python codes/RmLowMappibility_ChimericRead.py ${file} ${outdir}/${i}.sam ${mappabilityfile} ${thres}
	fi
done