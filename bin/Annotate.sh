#!/bin/bash
path=$1
mystart=$2
myend=$3
gtf=$4
codedir=$5
for((i=${mystart};i<=${myend};i++))
do
	file=${path}/${i}.sam
	python ${codedir}/Annotate.py ${file} ${gtf}
done
