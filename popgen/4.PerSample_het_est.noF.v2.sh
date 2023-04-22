#!/bin/bash
#SBATCH --job-name Het_Est   # Job name
#SBATCH -c 20                   # Use cpu
#SBATCH --mem-per-cpu 50G
#SBATCH -t 3-00:00:00             # Time limit hrs:min:sec
#SBATCH --mail-type TIME_LIMIT_80 
#SBATCH --array=1-20            # Array range
#SBATCH --mail-user georgette.femerling@sund.ku.dk
#SBATCH -o HetEst_Angsd_%A-%a.out   # Standard output and error log

module load htslib angsd

bamlist=$1
ref=$2
regions_file=$3
trans=$4 #only 0:with trans or 1:no trans

var=`date`;
echo -e "Start --> ${var} \n\n";

readarray bams < $bamlist
bam=${bams[$SLURM_ARRAY_TASK_ID]}
#bam=$bamlist

# readarray xjobs
# outname=$(echo $bam | rev | cut -f1 -d"/" | rev |cut -f1 -d".")
outname="$(basename -- $bam | cut -f1 -d".")"

echo -e "$outname, bam: $bam"
echo -e "$ref"

# Get per sample SAF
angsd -i $bam -GL 2 -dosaf 1 -anc $ref -ref $ref \
    -minQ 20 -minmapq 30 -rf $regions_file \
    -noTrans $trans \
    -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
    -nThreads 20 -out ${outname}.angsdput

# Run RealSFS

realSFS ${outname}.angsdput.saf.idx -P 20 -fold 1 -bootstrap 300 > ${outname}.folded.b300.est.ml 

var=`date`;
echo -e "End --> ${var} \n\n";