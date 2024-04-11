#!/bin/bash
#SBATCH --job-name=bismark_deduplicate_n_extract
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=06:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load bowtie2
module load samtools
module load python 
source activate e14

#make appropriate directories for output of deduplicate 
deduplicate_input=/global/scratch/users/chandlersutherland/e16/${sample}/em
mkdir -p $deduplicate_input/deduplicate

#establish the extraction output directory 
extraction_output=$deduplicate_input/extract
mkdir -p $extraction_output

#define the input files 
cd $deduplicate_input/bismark 
a=$(find . -type f -name '*_1_val_1_bismark_bt2_pe.bam')
prefix=$(basename -s .bam $a)

#define cores to partition to each replicate
num_reps=$(find . -type f -name '*_1_val_1_bismark_bt2_pe.bam' | wc -l)
cores=$(expr $SLURM_NTASKS / $num_reps)

#define deduplicator/extractor function 
DEDUP_EXTRA () {

	#deduplicate bismark 
	deduplicate_bismark -p \
		--output_dir  $deduplicate_input/deduplicate \
		./${1}.bam
	echo "finished deduplication of ${1}" 
	
	bismark_methylation_extractor -p \
		--output $extraction_output \
		--ignore_r2 2 \
		--comprehensive \
		--bedGraph \
		--CX \
		--parallel $cores \
		$deduplicate_input/deduplicate/${1}.deduplicated.bam
		
	echo "finished extraction of ${1}" 
}

export deduplicate_input=$deduplicate_input
export extraction_output=$extraction_output
export cores=$cores
export -f DEDUP_EXTRA 

time parallel DEDUP_EXTRA ::: $prefix
