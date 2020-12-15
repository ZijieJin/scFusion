#!/bin/bash
path=$1
mystart=$2
myend=$3
codedir=$4
for((i=${mystart};i<=${myend};i++))
do
	file=${path}/${i}_geneanno.sam
	python ${codedir}/FindFusionSupport.py ${file} ${file%_*}_FusionSupport.txt
done
