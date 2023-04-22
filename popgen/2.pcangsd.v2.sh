#!/bin/bash
#SBATCH --job-name PCAngsd_minmaf  # Job name
#SBATCH -c 20                   # Use pu
#SBATCH --mem-per-cpu 10G
#SBATCH -t 3-10:00:00             # Time limit hrs:min:sec
#SBATCH -o PCAngsd_modern_historical.out    # Standard output and error log
#SBATCH --mail-type TIME_LIMIT_80,DONE 
#SBATCH --mail-user georgette.femerling@sund.ku.dk

module load htslib angsd gcc/11.2.0 python/v3.6.9 pcangsd R/4.0.3

BAMlist=$1
outname=$2
popinfo=$3
ref=$4
regions_file=$5
depth_file=$6
trans=$7 # can only be 0 or 1

n1=$(wc -l $BAMlist | cut -f1 -d" ")
n=$(awk -v n=$n1 'BEGIN{ printf "%.0f", n*0.75 }')

min=$(cut -f1 $depth_file)
max=$(cut -f2 $depth_file)

# getting GLs and Snp call in BEAGLE format
echo -e "Calculating Genotypelikelihoods and calling snps"
angsd -bam $BAMlist -GL 2 -doMajorMinor 1 -doMaf 2 -SNP_pval 1e-6 -minmaf 0.05 \
    -doGlf 2 \
    -ref $ref \
    -noTrans $trans -minInd $n \
    -minMapQ 30 -minQ 20 -rf $regions_file \
    -setMinDepth $min -setMaxDepth $max -doCounts 1 \
    -nThreads 20 -out $outname

echo -e "Running PCA with PCAngsd"
echo -e "START: $(date)"

# Run PCangsd
pcangsd -beagle ${outname}.beagle.gz -o $outname.pcangsd -threads 20
# Plot PCA
Rscript /home/qsv231/storage/scripts/2.pcangsd_visualize.r $outname.pcangsd.cov $popinfo $outname.PCA

# Run Inbreeding using pcangsd 
# echo -e "Estimating Individual Inbreeding coefficients"
# # Using EM algorithm based on ngsF from -1 to 1
# pcangsd -beagle ${outname}.beagle.gz -o $outname.pcangsd -inbreedSamples -threads 20

# # Run Admixture estimation using pcangsd
# for K in $(seq 2 6);
# do
# pcangsd -beagle ${outname}.beagle.gz -o $outname.pcangsd -admix -admix_K $K -threads 20
# done
echo -e "END: $(date)"
