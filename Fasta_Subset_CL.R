## ---------------------------
##
## Script name: Fasta_Subset_CL.R
##
## Purpose of script: Get Subset of sequences based on hmmsearch output with command line inputs
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
## load packages ------------------------
package_list<-c("dplyr","Biostrings")

load_pack <- function(x){
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE )
      require( i , character.only = TRUE , quietly = T)
    }
  }
}
load_pack(package_list)
#library(odseq)
cat("=======================================\nLoaded packages\n=======================================\n")

## get options ---------------------------
option_list = list(
  make_option( c("-f", "--file") , type = "character" , default=NULL, 
               help="string that matches dataset file names", metavar="character"),
  make_option( c("-n", "--name") , type = "character" , default=NULL, 
               help="file containing protein names", metavar="character"),
  make_option(c("-d", "--directory"), type="character", default=".", 
              help="output directory name (optional)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="TestFastaSubset.fa", 
              help="output file name (optional)", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Provide protein sequences with -f option", call.=FALSE)
}
if (is.null(opt$name)){
  print_help(opt_parser)
  stop("Provide protein names with -n option", call.=FALSE)
}

# opt$f <- "/Users/prigozhin/Dropbox/NLRomes/Maize_NLRome/Proteomes/*.fasta"
# opt$n <- "/Users/prigozhin/Dropbox/NLRomes/Maize_NLRome/NLRome/Zm_NBARC_HMM.list"
## set working directory
files <- list.files(path = dirname(opt$f),pattern = basename(opt$f))
ID <- opt$n

read_delim(ID,delim = ' ', col_names = c("Gene")) -> My_Genes
#My_Genes %>% filter(grepl("P001",Gene)) -> My_Genes

## write a function for finding index of a list of names in an alignment
get_seqs <- function(list, ali){
  return(which(grepl(ali@range@NAMES,list)))
}
get_seqs(My_Genes,a)
grepl(pattern = a@ranges@NAMES,x = My_Genes)
sapply(a@ranges@NAMES, grepl, My_Genes, ignore.case = T)
which(a@ranges@NAMES %in% My_Genes$Gene)

#files <- Vertebrates$File 
#files <- Invertebrates$File
read_file(files[1])
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
