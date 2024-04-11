#!/bin/bash
#SBATCH --job-name=methylation_processing
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:01:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

cd /global/home/users/chandlersutherland/e16/

#Define the input directory, with fastq files ready for trimming
base="/global/scratch/users/chandlersutherland/e16"

#if have not run before, prepare each genome for bismark 
while read sample; do sbatch --job-name=$sample.genome --export=sample=$sample,base=$base \
-A fc_kvkallow methylation/bismark_genome_preparation.bash; done < sample.txt 

#trim the reads using trim galore 
while read sample; do sbatch --job-name=$sample.trim_galore --export=base=$base,sample=$sample \
 -A fc_kvkallow methylation/trim_galore.sh; done < sample.txt

#run bismark alignment 
while read sample; do sbatch --job-name=$sample.bismark --export=base=$base,sample=$sample -A co_minium \
methylation/bismark.bash; done < sample.txt 

#check read coverage of NLRs, generating a clean_coverage.tsv file for the bismark output  
while read sample; do sbatch --job-name=$sample.nlr_coverage --export=sample=$sample -A co_minium \
methylation/samtools_coverage_runner.bash; done < sample.txt

#deduplicate reads 
while read sample; do sbatch --job-name=$sample.dedup --export=sample=$sample -A co_minium \
methylation/bismark_deduplicate.bash; done < sample.txt

#finish extraction to coverage files, and filter to >5 reads per cytosine 
while read sample; do sbatch --job-name=$sample.extract --export=sample=$sample\
 -A co_minium methylation/bismark_2_bedgraph.bash; done < samples.txt

#calculate the per gene methylation percentage for each context 
while read sample; do sbatch --job-name=$sample.per_exon_meth --sample=$sample -A co_minium \
methylation/per_exon_methylation_runner.bash; done < samples.txt 

#generate bigwig files for viewing in IGV of the high coverage extracted cytosines 
# while read sample; do sbatch --job-name=$sample.bw --export=sample=$sample -A co_minium \
# methylation/bed_2_bw.bash; done < samples.txt 