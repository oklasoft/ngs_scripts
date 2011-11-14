#!/bin/bash
source /etc/profile.d/*.sh
module load vcftools

PID=$$
base=`basename ${1} .vcf`

vcftools --vcf ${1} --get-INFO AB --get-INFO AC --get-INFO AF --get-INFO AN --get-INFO BaseQRankSum --get-INFO DP --get-INFO Dels --get-INFO FS --get-INFO HRun --get-INFO HaplotypeScore --get-INFO InbreedingCoeff --get-INFO MQ --get-INFO MQ0 --get-INFO MQRankSum --get-INFO QD --get-INFO ReadPosRankSum --get-INFO SB --get-INFO VQSLOD --site-quality --out ${PID}
awk '{print $3}' ${PID}.lqual > temp1.${PID}
paste ${PID}.INFO temp1.${PID} > temp2.${PID}
echo "TAG	dbSNP" > temp3.${PID}
awk '{OFS="\t"; if (!/^#/){print $7, $3}}' ${1} >> temp3.${PID}
paste temp2.${PID} temp3.${PID} > temp4.${PID}
grep -v \? temp4.${PID} > ${base}_all_info.txt

rm temp1.${PID}
rm temp2.${PID}
rm temp3.${PID}
rm temp4.${PID}
rm ${PID}.INFO
rm ${PID}.lqual
rm ${PID}.log
