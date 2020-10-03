## ---------------------------
##
## Script name: Fasta_Subset.R
##
## Purpose of script: Get Subset of sequences based on hmmsearch output
##
## Author: Daniil Prigozhin
##
## Date Created: 2020-08-21
##
## Copyright (c) Daniil Prigozhin, 2020
## Email: daniilprigozhin@lbl.gov
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------
## load packages

require(tidyverse)
require(Biostrings)

## set working directory

setwd("~/Dropbox/MaizeNLRome/Proteomes/")  
ID <-"../NLRome/Zm_NBARC_HMM.list"

read_delim(ID,delim = ' ', col_names = c("Gene"))->My_Genes
My_Genes %>% filter(grepl("P001",Gene)) -> My_Genes


files <- list.files(pattern = ".*.fasta$")
#files <- Vertebrates$File 
#files <- Invertebrates$File

subset<-vector("list",length = length(files))
ii<-1
for (ii in seq_along(files)){
  print(paste0("Looking at file ",ii, " of ", length(files), " (",files[[ii]],")..."))
  a<-readAAStringSet(files[[ii]])
  a@ranges@NAMES %>% str_remove(pattern = "^..\\|") %>% str_remove("\\|.*$") ->a@ranges@NAMES
  a[a@ranges@NAMES %in% My_Genes$Gene]->subset[[ii]]
}

All_NLR<-subset[[1]]
for (jj in 2:length(subset)){All_NLR<-c(All_NLR,subset[[jj]])}
All_NLR
rm(subset,a,b,c)
getwd()
My_Genes %>% distinct()

writeXStringSet(All_NLR,"../NLRome/Zm_NBARC_Hits_P001.fasta")
