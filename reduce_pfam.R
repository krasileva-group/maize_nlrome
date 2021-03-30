## ---------------------------
##
## Script name: reduce_pfam.R 
##
## Purpose of script:
##
## Author: Daniil Prigozhin
##
## Date Created: 2021-03-04
##
## Copyright (c) Daniil Prigozhin, 2021
## Email: daniilprigozhin@lbl.gov
##
## ---------------------------
##
## Notes: This script works inside snakemake routine to remove redundant hits in hmmsearch --domtblout output.
##   
##
## ---------------------------
## load packages

#snakemake@input[["pfam"]]
require(tidyverse)
#require(parallel)

e_val_max <- 10
hmm_frac_min <- 0.3
#threads <- 4
max_overlap <-10

max_overlap <- snakemake@params[["max_overlap"]]
e_val_max <- snakemake@params[["e_val_max"]]
hmm_frac_min <- snakemake@params[["hmm_frac_min"]]
#threads <- snakemake@threads

#setwd("~/Dropbox/NLRomes/Maize_NLRome/HMM_search_pbNB-ARC/")
tblin <- read_delim(snakemake@input[["pfam"]],
                    delim = " ",
                    col_names = c("target_name","t_accession","tlen","query_name","q_accession","qlen",
                                  "fullseq_Evalue","fullseq_score","fullseq_bias",
                                  "dom_N","dom_of","dom_cEvalue","dom_iEvalue","dom_score","dom_bias",
                                  "hmm_from","hmm_to","ali_from","ali_to","env_from","env_to","acc","description_of_target"),
                    comment = "#",
                    trim_ws = T)

# tblin <- read_delim("all_samples.pbNB-ARC.Pfam_scan.out",
#                     delim = " ",
#                     col_names = c("target_name","t_accession","tlen","query_name","q_accession","qlen",
#                                   "fullseq_Evalue","fullseq_score","fullseq_bias",
#                                   "dom_N","dom_of","dom_cEvalue","dom_iEvalue","dom_score","dom_bias",
#                                   "hmm_from","hmm_to","ali_from","ali_to","env_from","env_to","acc","description_of_target"),
#                     comment = "#",
#                     trim_ws = T)
cat("==============================================================\nRead table:\n")
 tblin
## Apply e-value cutoff and hmm fraction minimum:
dom_tbl <- tblin %>% 
          select(c("target_name","tlen","query_name","qlen","fullseq_Evalue",
                   "dom_N","dom_of","dom_cEvalue",
                   "hmm_from","hmm_to","ali_from","ali_to","env_from","env_to")) %>%
          mutate(hmm_frac = (hmm_to-hmm_from)/qlen) %>% 
          arrange(target_name,env_from) %>%
          filter(fullseq_Evalue < e_val_max, hmm_frac > hmm_frac_min) %>%
          mutate(target_name = as.factor(target_name))
cat("==============================================================\nFilterred table by e-value and hmm coverage:\n")
dom_tbl

## Check length of resulting table >1???

reduced<- vector("list",length = nrow(dom_tbl))

line1 <- dom_tbl[15,]
line2 <- dom_tbl[16,]
jj <- 16
line1 <- dom_tbl[1,]
for (jj in 2:nrow(dom_tbl)){
  line2 <- dom_tbl[jj,]
  if (line1[[1,1]]==line2[[1,1]]){
    Ar = line1[[1,14]]
    Bl = line2[[1,13]]
    if (Ar - Bl > max_overlap){
      ##pick winner and store in line1
      Ae = line1[[1,5]]
      Be = line2[[1,5]]
      if (Be < Ae) { line1 <- line2 }
    }else{
      reduced[[jj-1]]<-line1
      line1 <- line2
    }
  }else{
    reduced[[jj-1]]<-line1
    line1 <- line2
  }
}
reduced_hits <- reduced[-which(sapply(reduced, is.null))] %>% map_df(.f=identity)
cat("==============================================================\nRemoved lower scoring overlapping hmm hits:\n")
reduced_hits
reduced_hits %>% write_delim(snakemake@output[["pfam"]],delim = "\t")

# remove_nesting <- function(gene_name,table){
#   gene_tbl <- table %>% filter(target_name==gene_name)
#   list <- vector()
#   for (ii in 1:gene_tbl[[1,2]]){
#     feature <- gene_tbl %>% filter(gene_tbl$env_from<=ii & gene_tbl$env_to>=ii)
#     if (nrow(feature)>0) {list<- rbind(list,feature %>% filter(fullseq_Evalue == min(fullseq_Evalue)))}
#     }
#   return(list%>%distinct())
# }
# #remove_nesting("ZM00001EB030390_P001",dom_tbl)
# 
# reduced_list <- mclapply(levels(dom_tbl$target_name), remove_nesting, table=dom_tbl,mc.cores = threads)
# reduced <- reduced_list[[1]]
# for (i in 2:length(reduced_list)){reduced <- rbind(reduced,reduced_list[[i]])}
# reduced
# reduced %>% write_delim(snakemake@output[["pfam"]],delim = "\t")
# #save.image("DEBUG.RData")
