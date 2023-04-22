#!/bin/bash
#SBATCH --job-name PSMC-bootstrap   # Job name
#SBATCH -c 10                   # Use cpu
#SBATCH --mem-per-cpu 100G
#SBATCH -t 3-00:00:00             # Time limit hrs:min:sec 
#SBATCH --mail-user georgette.femerling@sund.ku.dk
#SBATCH -o PSMC_75x_2.out  # Standard output and error log

threads=${SLURM_CPUS_PER_TASK}

echo "using $threads cpus."

# module load htslib
# module load bcftools 
# module load psmc

module load samtools
# conda activate /projects/mjolnir1/apps/conda/psmc-0.6.5
# conda activate /projects/mjolnir1/apps/conda/gnuplot-5.4.3

bam=$1
Ref="/home/qsv231/storage/Ref_genomes/Terpsiphone_corvina_B10K/D2102046629.gapcloser.fasta"
#/groups/hologenomics/hernan/data/GNRD_data/ref_genomes/target_species/Terpsiphone_corvina_B10K/D2102046629.gapcloser.fasta"
#scaffs="/groups/hologenomics/gfemer/data/Data/HighCoverage_Modern/Scaffolds_newref.txt"
scaffs="/groups/hologenomics/gfemer/data/Real_Reference_results/regions.txt"

base=$(basename $bam .merged.dedup.realigned.bam)
#coverage=$(cut -f2 $base.coverage)
coverage=75

echo -e "$base --> $coverage"
avcov=${coverage%.*}

d=$(expr $avcov / 3)
D=$(expr $avcov \* 2)
echo -e "d = $d , D = $D"

PSMC="/projects/mjolnir1/apps/conda/psmc-0.6.5/bin/"

# echo -e "Getting consensus diploid sequence"
# samtools mpileup -C50 -uf $Ref $bam | bcftools call -c | vcfutils.pl vcf2fq -d $d -D $D | gzip > $base.consensus.diploid.fq.gz
# echo -e "finished getting consensus sequence"

# echo -e "Running PSMC"
# echo "fq2psmcfa"
# fq2psmcfa -q20 $base.consensus.diploid.fq.gz > $base.diploid.psmcfa

# echo "split fasta"
# splitfa $base.diploid.psmcfa > $base.split.psmcfa

# echo "psmc - Avian pops values" 
# psmc -N30 -t5 -r5 -p "4+30*2+4+6+10" -o $base.diploid.psmc $base.diploid.psmcfa

echo "psmc - bootstrap 100 runs" 
# N30 –t5 –r5 –p 4+30∗2+4+6+10 based on KrystynaNadachowska-Brzyska paper
# seq 100 | xargs -i echo ${PSMC}psmc -N30 -t5 -r5 -b -p "4+30*2+4+6+10" -o round-{}.psmc $base.split.psmcfa | sh

seq 100 | xargs -P 10 -Iz ${PSMC}psmc -N30 -t5 -r5 -b -p "4+30*2+4+6+10" -o round-{z}.psmc $base.split.psmcfa

echo "Combine all runs"
cat $base.diploid.psmc round-*.psmc > $base.combined.psmc

echo "History"
${PSMC}psmc2history.pl $base.diploid.par | ${PSMC}history2ms.pl > ms-cmd.sh

echo "plot"
${PSMC}psmc_plot.pl -u 4.6e-9 -g 2 -p $base.combined $base.combined.psmc

# echo "step 4. plot"
# psmc_plot.pl -u 4.6e-9 -g 2 -T $base -G $base2.diploid $base2.diploid.psmc

echo -e "Done PSMC for $base"
