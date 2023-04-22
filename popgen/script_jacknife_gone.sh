#!/bin/bash
#SBATCH --job-name GONE  # Job name
#SBATCH -c 10                   # Use cpu
#SBATCH --mem-per-cpu 20G
#SBATCH -t 10:00:00             # Time limit hrs:min:sec
#SBATCH --array=0-17
#SBATCH -o GONE_Recent_demography-3cM-%A-%a.out

## Jackknife refers to removing one observation in the data set
## goal is to run the GONE analysis removing one sample at a time

module load plink/1.9.0 gcc/11.2.0 R/4.0.3

## GL of scaffolds of > 1MB with ANGSD - per sample
plinkfile=$1 # Has to have chr names

#cut -f1,2 -d" " $plinkfile.ped > indv_list.txt

readarray ilist < /home/qsv231/storage/Final_results/Paper_Revision_analyses/GONE_jackknife/SPF_modern_individuals.txt
indv=${ilist[$SLURM_ARRAY_TASK_ID]}

echo $indv

mkdir jackknife_${SLURM_ARRAY_TASK_ID}

cp -r PROGRAMMES jackknife_${SLURM_ARRAY_TASK_ID}
cp INPUT_PARAMETERS_FILE jackknife_${SLURM_ARRAY_TASK_ID} 

cd jackknife_${SLURM_ARRAY_TASK_ID}

printf "%s" "$indv" > remove.txt

plink --file ../$plinkfile --remove remove.txt --aec --recode --out plink.$SLURM_ARRAY_TASK_ID.tmp

awk 'BEGIN{chr="";n=0} {if($1!=chr){chr=$1;n+=1;} print n"\t"$2"\t"$3"\t"$4}' plink.$SLURM_ARRAY_TASK_ID.tmp.map > plink.$SLURM_ARRAY_TASK_ID.map
mv plink.$SLURM_ARRAY_TASK_ID.tmp.ped plink.$SLURM_ARRAY_TASK_ID.ped

/home/qsv231/storage/tools/GONE/script_GONE.sh plink.$SLURM_ARRAY_TASK_ID
