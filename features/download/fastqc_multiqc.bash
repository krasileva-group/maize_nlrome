#!/bin/bash
#SBATCH --job-name=fastqc_raw
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load python 
module load fastqc 
source activate e16 

#rna 
OUTPUT_DIR=$base/$sample/fastqc
mkdir -p $OUTPUT_DIR
cd $base/$sample/
FILES=$(find . -type f -wholename './rna_*/*fq')
fastqc -o $OUTPUT_DIR -t 24 $FILES 

#em_raw 
#OUTPUT_DIR=$base/$sample/em/fastqc
#mkdir -p $OUTPUT_DIR
#cd $base/$sample/em/
#FILES=$(find . -type f -maxdepth 1 -name '*fastq' -print)

#fastqc -o $OUTPUT_DIR -t 24 $FILES 

