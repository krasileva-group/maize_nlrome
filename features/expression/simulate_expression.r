#!/usr/bin/env Rscript
library(polyester)
library(Biostrings)

#read in the list of samples 
x <- scan("/global/home/users/chandlersutherland/e16/sample.txt", what="", sep="\n")

#create a function that takes the sample name as input, and simulates RNAsequences of the NLRs, outputing them to a simulated_nlr folder within the rna subfolder
simulator <- function(sample) {
  #get the fasta file and read in as a DNAString Set 
  primary_fasta_file=Sys.glob(paste('/global/scratch/users/chandlersutherland/e16/',sample,'/genome/*.canonical.cds.fa.nlrs.fa', sep=""))
  primary_fasta=readDNAStringSet(primary_fasta_file)
  
  #calculate the number of transcripts to set the number of rows 
  sequences=length(primary_fasta)
  print(paste("There are", sequences, 'in sample', sample, 'to be simulated.'))
  
  #The polyester library simulates RNAseq reads based on an input fold change matrix and a fasta file containing transcripts to simulate from
  #Generate a fold change matrix, with only "1" values since I am just interested in coverage, not differential expression 
  #initialize the fold change matrix and the depth, aiming for 20x coverage 
  primary_fold_changes=matrix(c(1), nrow=sequences, ncol=1)
  primary_readspertx=round(20*width(primary_fasta)/75)
  
  #simulate the experiment and output into a directory called simulated 
  simulate_experiment(primary_fasta_file, reads_per_transcript=primary_readspertx, 
                      num_reps=c(2), fold_changes=primary_fold_changes, 
                      readlen=75, paired=TRUE, 
                      outdir=paste('/global/scratch/users/chandlersutherland/e16/', sample, '/rna_tip/simulated_nlr', sep=""))

}


for (i in 1:length(x)){
  simulator(x[i])
}

	