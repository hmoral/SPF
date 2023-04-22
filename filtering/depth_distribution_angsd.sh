#!/bin/bash
#SBATCH --job-name Depth-distribution   # Job name
#SBATCH -c 10                   # Use cpu
#SBATCH --mem-per-cpu 20G
#SBATCH -t 3-00:00:00             # Time limit hrs:min:sec
#SBATCH --mail-type end,fail 
#SBATCH --mail-user georgette.femerling@sund.ku.dk
#SBATCH -o Depth-distribution.out

module load htslib angsd

echo "angsd v0.940-dirty (htslib: 1.16)"

bamlist=$1
outname=$2

angsd -bam $bamlist -nThreads 10 -doDepth 1 -out ${outname} -doCounts 1 -minMapQ 30 -minQ 20 -maxdepth 3000 -howoften 1000000

Rscript ~/storage/scripts/depth_distribution.r ${outname}.depthGlobal $outname