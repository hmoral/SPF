#!/bin/bash

LOGseed=`echo $RANDOM`
LOGdate=$(date '+%d_%m_%Y_%H_%M_%S')
LOG=LOG_snpEFF_extractFields_merged_${LOGdate}_LOGseed${LOGseed}.log
ERR=LOG_snpEFF_extractFields_merged_${LOGdate}_LOGseed${LOGseed}.err.log


VCFlist=$1

OUT=$2

echo `date` " - START" $VCFlist  >> $LOG 2>> $ERR

while read VCF;do

echo `date` " - indexing" $VCF  >> $LOG 2>> $ERR
zcat $VCF | bgzip -c  > tmp_seed${LOGseed}_${VCF}
bcftools index tmp_seed${LOGseed}_${VCF}

done < ${VCFlist}

echo `date` "- merging" $VCFlist  >> $LOG 2>> $ERR

ls tmp_seed${LOGseed}_*gz > tmp_vcfLIST_seed${LOGseed}.txt

bcftools merge --threads 10 -l tmp_vcfLIST_seed${LOGseed}.txt -Oz -o all_vcfMergeBCFtools_${OUT}.vcf.gz

echo `date` "- extractFields" $VCFlist  >> $LOG 2>> $ERR

java -Xmx40G -jar SnpSift.jar extractFields -s "," -e "."  \
all_vcfMergeBCFtools_${OUT}.vcf.gz \
CHROM POS REF ALT DP "GEN[*].GT" ANN[*].ERRORS ANN[*].EFFECT EFF[*].EFFECT EFF[*].IMPACT LOF[*].GENE LOF[*].PERC ANN[*].GENE | \
awk '{ if ( $7 != "." ) { print $0; } }' | gzip -c >  1_out_snpeff_vcfMergeBCFtools_${OUT}.tsv.gz

rm tmp_seed${LOGseed}*_${OUT}*.vcf.gz

echo `date` " - done" $VCFlist  >> $LOG 2>> $ERR


