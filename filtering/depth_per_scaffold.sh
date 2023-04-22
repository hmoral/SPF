#!/bin/bash
#SBATCH --job-name Depth_perScaffold   # Job name
#SBATCH -c 10                   # Use cpu
#SBATCH -t 3-00:00:00             # Time limit hrs:min:sec
#SBATCH --mail-type TIME_LIMIT_50,TIME_LIMIT_80 
#SBATCH --mail-user georgette.femerling@sund.ku.dk
#SBATCH -o Depth-%j.out

module load htslib angsd samtools

BAMLIST=$1
OUTDIR=$2

BAMS=(`cat $BAMLIST`)
mkdir -p ${OUTDIR}
echo -e "START: $(date)"

SCAFS=`samtools view -H ${BAMS[1]} | grep @SQ | cut -f 2,3 | sed -e 's/SN://g' -e 's/LN://g' | cut -f 1 ` # awk '$2>500' por si se quiere restringir a cierto length
#SCAFS=(`cat Scaffolds_que_faltan.txt`)

for i in ${SCAFS[@]}
do
    echo -e "..... Scaffold $i ....."
    if [ -s "${OUTDIR}/${i}.depthGlobal" ];
    then
        echo -e "Scaffold $i files already exist"
        continue
    fi
    echo -e "at: $(date)"
    angsd -nThreads 10 -howoften 1000000 -minMapQ 30 -minQ 20 -doDepth 1 -doCounts 1 -dumpCounts 2 -maxdepth 200 -b $BAMLIST -out ${OUTDIR}/${i} -r $i
done 

echo -e "END: $(date)"

echo "Done calculating depth per scaffold"

echo -e "START: $(date)"
script /home/qsv231/storage/scripts/Depth_filter.r
echo -e "END: $(date)"

echo "Done with depth analysis"