#!/bin/bash
#SBATCH --job-name=coverage_py_runner
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --time=00:10:00
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

#this script passes a bam file to the python script samtools_coverage, which outputs a file $BASENAME_clean_coverage.tsv which gives the mean depth over the NLRs
#can easily be made into a for loop to pass multiple files through 
#pass variable coverage_input, which should have a directory of sam files ready for coverage processing 

module load python
module load samtools/1.14
module load parallel
#very important to have this version 


coverage_input=/global/scratch/users/chandlersutherland/e16/${sample}/rna_tip/STAR
cd $coverage_input
bams=$(find . -type f -name '*sortedByCoord.out.bam')

#input is a bam file 
processing () {
	#first, convert to bam, sort and index for samtools coverage to run appropriately 
	basename=$(basename ${1}) 
	samtools index $basename
	
	#run the coverage file 
	#pass two arguments, the NLR bed file and the sample
	NLR_bed=/global/scratch/users/chandlersutherland/e16/${sample}/genome/*_NLR.bed
	python $HOME/e16/methylation/samtools_coverage.py  $NLR_bed $basename

	#clean up working directory 
	mkdir -p $coverage_input/coverage
	mv *_clean_coverage.tsv $coverage_input/coverage
	rm *coverage.tsv 
}

export coverage_input=$coverage_input
export sample=$sample
export -f processing

parallel processing ::: $bams