#!/bin/bash
#SBATCH --job-name=expression_processing
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:30:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

cd /global/home/users/chandlersutherland/e16/
base="/global/scratch/users/chandlersutherland/e16"
module load python 

echo "beginning re-processing of ${sample}"

#first, download problematic rna files (don't have a STAR output) 
#python rna_download.py $sample >> /global/scratch/users/chandlersutherland/e16/$sample/rna_download.log
#actually do this on a download node instead

echo "finished download script"

#unpigz files 
sbatch --job-name=$sample.unpigz --export=base=$base,sample=$sample -A co_minium -p savio4_htc \
--qos minium_htc4_normal unpigz.sh

echo "unpigzed"

#wait until complete to launch next jobs 
until [ ! -f $base/$sample/rna_root/*.fastq ]
do

    sleep 15
done
echo "proceeding to STAR, unpigz finished" 

#run STAR in quant mode 
while read tissue; do sbatch --job-name=$sample.$tissue.STAR --export=base=$base,sample=$sample,tissue=$tissue \
-A co_minium -p savio4_htc --qos minium_htc4_normal expression/STAR.bash; done < tissues.txt 

#wait until complete 
#until [ ! -f $base/$sample/rna_root/*.tab ]
#do
#    sleep 5
#done
#echo "STAR finished, make conglom output file"

#create output file 
#python expression/tpm_calc.py $sample 

#QC coverage of each NLR 
#while read sample; do sbatch --job-name=$sample.coverage --export=base=$base,sample=$sample \
#-A co_minium --partition savio4_htc expression/samtools_coverage_runner.bash; done < sample.txt 
