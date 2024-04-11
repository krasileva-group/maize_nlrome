#!/bin/bash
#SBATCH --job-name=bismark_genome_preparation
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load bowtie2
module load samtools
module load Bismark 
module load python 

genome_dir=$base/$sample/genome
cd $genome_dir

#add lambda and puc19 controls to the genome before bismark conversion 
lambda=/global/scratch/users/chandlersutherland/e16/em_control/lambda.fa
puc19=/global/scratch/users/chandlersutherland/e16/em_control/pUC19.fa

genome=$(find . -maxdepth 1 -type f -name "*.fa")
base=$(basename $genome .fa)

#bismark requires a genome directory as an input, so make a sub directory to put lambda/pUC19 genome in 
mkdir -p $genome_dir/bismark

cat $genome $lambda $puc19 > $genome_dir/bismark/${base}_meth.fa 
echo 'added control sequences' 

#run genome preparation
time bismark_genome_preparation \
     --verbose \
     --large-index \
     --bowtie2 \
     --parallel 12 \
     --path_to_aligner /global/home/groups/consultsw/sl-7.x86_64/modules/bowtie2/2.3.4.1 \
     $genome_dir/bismark

echo "finished genome preparation for sample $sample"
