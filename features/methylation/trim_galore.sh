#!/bin/bash
#SBATCH --job-name=trim
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8 
#SBATCH --time=04:30:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

#currently written on savio, but could easily be done on savio4_htc with small core support (only uses 8 per sample)
module load cutadapt
module load fastqc 
module load python 
module load parallel
module load pigz

#set the trim directory, output directory, and input directory, containing the untrimmed fastq files
TRIM_DIR=/global/home/users/chandlersutherland/programs/TrimGalore-0.6.6
trim_input=$base/$sample/em
trim_output=$trim_input/trimmed

mkdir -p $trim_output
cd $trim_input

#unpigz the directory 
#unpigz * 

#run trim galore in default, paired end mode 
FILES=$(find . -maxdepth 1 -type f -name '*_1.fastq' -print)
REPS=$(basename -a -s _1.fastq $FILES)

#write function to pass to parallel, input is accession number  
#trim galore doesn't perform well after 4 cores 
trim(){
	$TRIM_DIR/trim_galore -o $trim_output \
		--fastqc \
		--illumina \
		--cores 4 \
		--paired $trim_input/$1_1.fastq $trim_input/$1_2.fastq
	echo "finished trimming $1"
}

#export necessary path variables and the trim function 
export TRIM_DIR=$TRIM_DIR
export trim_input=$trim_input
export trim_output=$trim_output
export -f trim

time parallel trim ::: $REPS

echo 'finished bisulfite trimming of $sample'
