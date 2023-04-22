#!/bin/bash
#SBATCH --job-name Rmdups_realing   # Job name
#SBATCH -c 1                  # Use two cpu
#SBATCH --mem-per-cpu 40G
#SBATCH -t 1-00:00:00             # Time limit hrs:min:sec
#SBATCH -o Rmdup_%A-%a.out    # Standard output and error log
#SBATCH --array=0-12                # Array range
#SBATCH -p hologenomics

# Load Required Modules
module load htslib samtools java picard GATK/v3.8.1

# Get arguments
ref="/groups/hologenomics/hernan/data/GNRD_data/ref_genomes/target_species/Terpsiphone_corvina_B10K/D2102046629.gapcloser.fasta"

# Prepare reference files needed
samtools faidx $ref
samtools dict $ref 

bams=($(ls -d $PWD/Bamfiles/Merged/*.bam))
bam=${bams[$SLURM_ARRAY_TASK_ID]}
sample="$(basename -- $bam | cut -f1 -d".")"

if [[ -f "$PWD/Bamfiles/Realigned/${sample}.merged.dedup.realigned.bai" ]]; then
    echo "$sample exists."
    exit 0;
fi

mkdir Bamfiles/Deduped Bamfiles/Realigned
echo -e $sample

echo -e "Marking duplicates"
picard MarkDuplicates INPUT=$bam \
    OUTPUT=Bamfiles/Deduped/${sample}.merged.dedup.bam METRICS_FILE=Bamfiles/Deduped/${sample}.metrics \
    ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true \

samtools index Bamfiles/Deduped/${sample}.merged.dedup.bam
echo -e "Done removing duplicates from $sample bam."

echo -e "Realingning around indels"
gatk -T RealignerTargetCreator \
    -R $ref -I Bamfiles/Deduped/${sample}.merged.dedup.bam \
    -o Bamfiles/Realigned/${sample}_realignertargetcreator.intervals

gatk -T IndelRealigner \
    -R $ref -I Bamfiles/Deduped/${sample}.merged.dedup.bam \
    -targetIntervals Bamfiles/Realigned/${sample}_realignertargetcreator.intervals \
    -o Bamfiles/Realigned/${sample}.merged.dedup.realigned.bam

# samtools index Bamfiles/Realigned/${sample}.merged.dedup.realigned.bam
