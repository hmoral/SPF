#!/bin/bash

VCFlist=$1

LOGseed=`echo $RANDOM`
LOGdate=$(date '+%d_%m_%Y_%H_%M_%S')

LOG=1_logs_snpEFF/LOG_${LOGdate}_LOGseed${LOGseed}.log
ERR=1_logs_snpEFF/LOG_${LOGdate}_LOGseed${LOGseed}.err.log

mkdir -p 1_logs_snpEFF

echo `date` "start file prep for snpEff"   >> $LOG 2>> $ERR

while read VCF;do

OUT=`echo $VCF | sed s'/.vcf.gz//'`
echo "java -Xmx12G -jar snpEff.jar eff -lof -treatAllAsProteinCoding -v -s snpeff_stats_${OUT}.html -csvStats snpeff_stats_${VCF}.csv SPF ${VCF} | gzip -c > ${OUT}_snpeff.vcf.gz" >> launch_jobs_snpEff_${LOGseed}.txt
done < ${VCFlist}

echo `date` "start parallel  for snpEff"   >> $LOG 2>> $ERR

parallel --jobs 10 < launch_jobs_snpEff_${LOGseed}.txt

echo `date` "start multiQC" >> $LOG 2>> $ERR

multiqc -n snpEFF_multiqc_report_${LOGdate} -o ./ .  >> $LOG 2>> $ERR


echo `date` "END" >> $LOG 2>> $ERR

