#!/bin/bash
#SBATCH --job-name=per_gene_meth_runner
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

#this script first checks to see how many bio reps and contexts are unprocessed into per gene (by exon) counts
#then, passes unfinished per position methylation files (.cov Bismark output) and exon coordintes (in .bed format) to the python script per_exon_methylation.py


module load python

#point to the all exon bed file 
exon_positions=/global/scratch/users/chandlersutherland/e16/${sample}/genome/*_all_exon.bed

#where are the coverage files that have not been processed yet
cov_dir=/global/scratch/users/chandlersutherland/e16/${sample}/em/bedGraph_highcov

cd $cov_dir
rm -r *prefix

#First define the finished files
finished=$(find . -type f -name '*.tsv')
prefix=$(basename -s "_per_exon_met_${sample}.tsv" $finished)
for i in $prefix; do echo $i >> finished_prefix; done 



#get all of the files 
all=$(find . -type f -name '*.cov')
all_prefix=$(basename -s '.bed.gz.bismark.cov' $all)
for i in $all_prefix; do echo $i >> all_prefix; done 

#define the complement, aka all unknown files 
unfinished=$(comm -23 <(sort all_prefix) <(sort finished_prefix))

echo "$unfinished have not been completed. Starting meth extraction now" 

for f in $unfinished
do 
	time python $HOME/e16/methylation/per_exon_methylation.py $cov_dir/$f.bed.gz.bismark.cov $exon_positions 
done
