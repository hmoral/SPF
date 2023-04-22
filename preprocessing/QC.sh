#!/bin/bash
#SBATCH --job-name QC  # Job name
#SBATCH -c 10                   # Use cpu
#SBATCH -t 10:00:00             # Time limit hrs:min:sec
#SBATCH -o QC

module load openjdk perl fastqc

output_dir=$1
fq_file=$2

list=$(tr -s '\n ' ' ,' < $fq_file)

fastqc -o $output_dir -t 10 $list

multiqc $output_dir

