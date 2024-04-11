#!/bin/bash
#SBATCH --job-name=bismark
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load python 
module load bowtie2
module load Bismark
module load samtools
module load parallel 

trim_output=$base/$sample/em/trimmed
genome=$base/$sample/genome/bismark
bismark_output=$base/$sample/em/bismark
mkdir -p $bismark_output

#clean up anything from previous runs 
cd $trim_output 
rm -r *prefix

#need to export three variables, trim_output directory, which is our input directory, bismark_output directory, and the genome file 

#Check for completed files 
#First define the finished files
cd $bismark_output
finished=$(find . -type f -size +1b -name '*.txt') #if the report is written, bismark is complete 
prefix=$(basename -s "_1_val_1_bismark_bt2_PE_report.txt" $finished)

cd $trim_output
for i in $prefix; do echo $i >> finished_prefix; done 

#get all of the files 
a=$(find . -type f -name '*_1_val_1.fq')
accession=$(basename -s _1_val_1.fq $a)
for i in $accession; do echo $i >> all_prefix; done 

#define the complement, aka all unknown files 
unfinished=$(comm -23 <(sort all_prefix) <(sort finished_prefix))


BISMARK() {
    echo "beginning bismark on sample ${sample} ${1}"
    bismark --genome $genome \
	--temp_dir $SCRATCH \
	--output_dir $bismark_output \
	-p 4 \
	-1 "${1}"_1_val_1.fq -2 "${1}"_2_val_2.fq
}

export genome=$genome
export bismark_output=$bismark_output
export -f BISMARK

parallel BISMARK ::: $unfinished
