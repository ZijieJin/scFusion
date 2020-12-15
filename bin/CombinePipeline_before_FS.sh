#!/bin/bash

FilePath=$1
mystart=$2
myend=$3
gtffile=$4
mappabilityfile=$5
exonfile=$6
codedir=$7
mkdir -p ${FilePath}/StarMapping
mkdir -p ${FilePath}/ChimericOut
mkdir -p ${FilePath}/Expr/


for ((i=${mystart};i<=${myend};i++))
do
	file=`ls ${FilePath}/StarMapping/${i}/*Chimeric.out.sam`
	if [[ -n ${file} ]]; then
		python ${codedir}/RmLowMappibility_ChimericRead.py ${file} ${FilePath}/ChimericOut/${i}.sam ${mappabilityfile} 1
	fi
done

sh ${codedir}/Annotate.sh ${FilePath}/ChimericOut/ ${mystart} ${myend} ${gtffile} ${codedir}
sh ${codedir}/FindFusionSupport.sh ${FilePath}/ChimericOut/ ${mystart} ${myend} ${codedir}
sh ${codedir}/CalcRPKM.sh ${exonfile} ${FilePath}/StarMapping ${FilePath}/Expr/ ${mystart} ${myend}