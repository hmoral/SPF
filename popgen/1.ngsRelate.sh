#!/bin/bash
#SBATCH --job-name ngsRelate   # Job name
#SBATCH -c 20                   # Use n cpu
#SBATCH -t 3-00:00:00             # Time limit days-hrs:min:sec
#SBATCH -o ngsRelate.out    # Standard output and error log
#SBATCH --mail-type ALL,TIME_LIMIT_50 
#SBATCH --mail-user georgette.femerling@sund.ku.dk
#SBATCH -p hologenomics

module load htslib angsd python/v3.6.9 ngsRelate

# Read arguments
bamlist=$1
ids=$2
outname=$3
ref=$4
regions_file=$5
trans=$6 # can be 0 or 1 

# ids="/groups/hologenomics/gfemer/data/SPF_modern/sample_ids.txt"
n=$(wc -l $ids | cut -f1 -d" ")
echo -e "Number of individuals: $n"

# Filters: 
# No sexual chrs
# No sites with depth under or above 1% and 99% 

# ref="/groups/hologenomics/gfemer/data/Ref_genomes/Terpsiphone_cinnamomea/D1907004422.gapcloser.fasta.bg.gz"
# sites_file="/groups/hologenomics/gfemer/data/SPF_modern/Depth_filter_Sites/DepthFilteredsites.sex_depth.all.sites"
# regions_file="/groups/hologenomics/gfemer/data/Ref_genomes/Terpsiphone_cinnamomea/contigs_greatereq_10k.txt"

## First we generate a file with allele frequencies (angsdput.mafs.gz) and a file with genotype likelihoods (angsdput.glf.gz).
angsd -b $bamlist -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 2 -minmaf 0.05 -doGlf 3 \
    -minMapQ 30 -minQ 20 -minInd 17 -noTrans $trans \
    -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
    -setMinDepth 20 -setMaxDepth 242 -doCounts 1 \
    -ref $ref -rf $regions_file \
    -nThreads 20 -out $outname
    # -sites $sites_file

### Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
zcat $outname.mafs.gz | cut -f6 | sed 1d > $outname.freq

### run NgsRelate
ngsRelate -g $outname.glf.gz -n $n -f $outname.freq -O $outname.ngsRelate -z $ids -p 20

### Estimate Inbreeding
ngsRelate -g $outname.glf.gz -F 1 -f $outname.freq -n $n -O $outname.ngsRelate.inbreeding -p 20 