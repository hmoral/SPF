#!/bin/bash

#SBATCH --job-name Subsample_modern   # Job name
#SBATCH -c 5                   # Use cpu
#SBATCH --mem-per-cpu 5G
#SBATCH -t 3:00:00             # Time limit hrs:min:sec 
#SBATCH --mail-user georgette.femerling@sund.ku.dk
#SBATCH -o subsample_modern_%A-%a.out   # Standard output and error log
#SBATCH --array=13-18             # Array range
#SBATCH -p hologenomics

module load htslib samtools

BAMlist=$1
#coveragefile=$2
final_cov=$2
outfile="modern"

readarray bams < $BAMlist
bam=${bams[$SLURM_ARRAY_TASK_ID]}
name="$(basename -- $bam | cut -f1 -d".")"

# outfile=$1
# get coverage
samtools depth -q 20 -Q 30 -a $bam | awk -v name=$name '{sum+=$3} END { print name"\t"sum/NR}' >> $outfile.coverage
coveragefile=$outfile.coverage

coverage=$(grep $name $coveragefile | cut -f2)
frac=$(echo $coverage | awk -v cov=$final_cov '{val = cov / $1; print val}')

echo "name: $name, coverage: $coverage, fraction: $frac"

# sub sample to aprox 3x
samtools view -s $frac -b -@ 5 $bam > $name.merged.dedup.sorted.realigned.down.${final_cov}x.bam
samtools index $name.merged.dedup.sorted.realigned.down.${final_cov}x.bam