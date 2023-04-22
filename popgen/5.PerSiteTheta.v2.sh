#!/bin/bash
#SBATCH --job-name SFS-Theta  # Job name
#SBATCH -c 25                   # Use cpu
#SBATCH --mem-per-cpu 20G
#SBATCH -t 3-00:00:00             # Time limit hrs:min:sec
#SBATCH -o SFS-theta.folded.pi.out

module load samtools angsd

outname=$1

bamlist=$1
ref=$2
regions=$3
outname=$4
depth_file=$5
trans=$6

n1=$(wc -l $bamlist | cut -f1 -d" ")
n=$(awk -v n=$n1 'BEGIN{ printf "%.0f", n*0.75 }')

min=$(cut -f1 $depth_file)
max=$(cut -f2 $depth_file)

echo "Min Ind set to $n"
echo "-noTrans set to $trans"

angsd -bam $bamlist -GL 2 -dosaf 1 -anc $ref -ref $ref \
    -minQ 20 -minmapq 30 -rf $regions \
    -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
    -setMinDepth $min -setMaxDepth $max -doCounts 1 \
    -minInd $n -noTrans $trans \
    -nthreads 25 -out ${outname}

realSFS ${outname}.saf.idx -fold 1 -P 25 > ${outname}.sfs 
#-sites $sites

realSFS saf2theta ${outname}.saf.idx -fold 1 -sfs ${outname}.sfs -P 25 -outname ${outname}

thetaStat do_stat ${outname}.thetas.idx -win 1000 -step 1000 -outnames ${outname}.1kbwin-nonoverlapping.theta.thetasWindow.gz

thetaStat do_stat ${outname}.thetas.idx -win 50000 -step 1000 -outnames ${outname}.50kbwin-1kbstep.theta.thetasWindow.gz
