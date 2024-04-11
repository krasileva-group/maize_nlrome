#!/bin/bash
#SBATCH --job-name=rna_download
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load python 
module load parallel 
#module load pigz 


#python /global/home/users/chandlersutherland/e16/rna_download.py $sample 
sample=$(cat sample.txt)

#define download function, which takes in a accession name and begins rna_download.py, which downloads everything of all tissue types 
DWNLD (){
    python -u /global/home/users/chandlersutherland/e16/rna_download.py ${1} >> /global/scratch/users/chandlersutherland/e16/${1}/rna_download.log
    echo "finished $1"
	#unpigz /global/scratch/users/chandlersutherland/e16/${1}/rna_*/*.gz
}

#Download 26 accessions 
export -f DWNLD

parallel DWNLD ::: $sample
#while read sample; do DWNLD $sample; done < sample.txt 
