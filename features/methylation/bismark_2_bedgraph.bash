#!/bin/bash
#SBATCH --job-name=bismark_extractor
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=04:30:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load python 
module load Bismark
module load samtools
module load parallel 

#initiate directories 
input=/global/scratch/users/chandlersutherland/e16/${sample}/em
output=$input/bedGraph_highcov
mkdir -p $output

#function takes two positional arguments: $1 being the .txt file to extract from, and $2 being the output directory
#cutoff set to 5 means that only cytosines with >5 reads are outputted, reducing file size and simplifying future processing 
BISMARK_BEDGRAPH () {
	basename=$(basename $1 _1_val_1_bismark_bt2_pe.deduplicated.txt)
	bismark2bedGraph --output $basename.bed \
	--dir $output \
	--cutoff 5 \
	--CX \
	$1
    echo "finished ${1}"
}
export output=$output
export -f BISMARK_BEDGRAPH

CHX_in=$input/extract/CH*_context_*_1_val_1_bismark_bt2_pe.deduplicated.txt

parallel BISMARK_BEDGRAPH ::: $CHX_in

#need a separate function for CpG since it cannot take the --CX tag 
BISMARK_BEDGRAPH_CpG () {
	basename=$(basename $1 _1_val_1_bismark_bt2_pe.deduplicated.txt)
	bismark2bedGraph --output $basename.bed \
	--dir $output \
	--cutoff 5 \
	$1
    echo "finished ${1}"
}

export -f BISMARK_BEDGRAPH_CpG

CpG_in=$input/extract/CpG_context_*_1_val_1_bismark_bt2_pe.deduplicated.txt

parallel BISMARK_BEDGRAPH ::: $CpG_in 

cd $output
gunzip *.cov.gz 

echo "finished generating context specific coverage files for ${sample}"
