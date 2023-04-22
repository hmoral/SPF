#!/bin/bash
#SBATCH --job-name Admix   # Job name
#SBATCH -c 5                   # Use ten cpu
#SBATCH --mem-per-cpu 10G
#SBATCH -t 2-10:00:00             # Time limit hrs:min:sec
#SBATCH --mail-user georgette.femerling@sund.ku.dk
#SBATCH -o NGSadmix_%A-%a.out   # Standard output and error log
#SBATCH --array=2-6

module load htslib angsd evalAdmix

beagle=$1 # GL beagle file
b=$2 # bootstrap replicates
outname=$3 
minInd=$4

mkdir bestReps

K=$SLURM_ARRAY_TASK_ID # change this in the arrayx
mkdir K$K 
cd K$K

echo -e "START $K: $(date)"

# Run NGSAdmix
for j in $(seq 1 $b)
do
    echo -e "K $K, Run $j" 
    NGSadmix -likes ../${beagle} -K $K -P 5 -o $outname.admix.K${K}.run${j} -minMaf 0.05 -minInd $minInd
    # evalAdmix -beagle $beagle -fname $outname.admix.K${K}.run${j}.fopt.gz -qname $outname.admix.K${K}.run${j}.qopt -o $outname.admix.K${K}.run${j} -P 20
done

# Get bestrepetition (From Emily Cavil)
bestrep=$(for rep in $(seq 1 $b) ;do echo $rep $(tail -n 1 $outname.admix.K${K}.run${rep}.log | cut -f2 -d= | cut -f1 -d" " ); done | sort -k2,2nr | head -1 | cut -f1 -d" ");
cp $outname.admix.K${K}.run${bestrep}.qopt ../bestReps/$outname.admix.K${K}.run${bestrep}.bestrep.qopt

#(for log in `ls *.log`; do grep -Po 'like=\K[^ ]+' $log; done) > logfile

echo -e "END $K: $(date)"