#!/bin/bash
#SBATCH --job-name=sort_index
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

#export an input file to be sorted and indexed 
#pigz compress sorted and indexed file for easy scp 
module load samtools 
module load pigz

basename=$(basename $input .bam) 
dirname=$(dirname $input) 
cd $dirname 

samtools sort -o $basename.sorted.bam $input 
samtools index $basename.sorted.bam

pigz $basename.sorted.bam 
