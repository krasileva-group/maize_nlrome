## ---------------------------
##
## Script name: Zm_NLRome_PrintInitClades.R
##
## Purpose of script: snakemake routine for printing clades after a good partition is found
##
## Author: Daniil Prigozhin
##
## Date Created: 2021-03-23
##
## Copyright (c) Daniil Prigozhin, 2021
## Email: daniilprigozhin@lbl.gov
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------


## load packages
require("tidyverse")
require("treeio")
require("phangorn")
## load snakemake pointers ---------
raxml <- read.raxml(snakemake@input[["tree"]])
#raxml <- read.raxml("../RAxML_tree_pbNB-ARC/RAxML_bipartitionsBranchLabels.pbNB-ARC.Raxml.out")

partition <- read_delim(snakemake@input[["tsv"]],delim = "\t", col_names = TRUE)
#partition <- read_delim("../RAxML_tree_pbNB-ARC/Initial_Clades.pbNB-ARC.tsv",delim = "\t", col_names = TRUE)
list <- snakemake@output[["list"]]
dir <- snakemake@output[["dir"]]
if(!dir.exists(dir)){dir.create(dir)}

########Annotate Internal Nodes and save Bootstrap values---------------
cat ("========================================\nReading tree\n")
z <- mutate(as_tibble(raxml), label=if_else(is.na(label),paste0("Int",node),label))
bs <- z %>% select(label,bootstrap)
#write.tree(as.phylo(z),"ZMonly.100.ztree.txt")
########Root midpoint -----------------------
cat ("========================================\nRooting at midpoint\n")

tree <- midpoint(as.phylo(z))
x <- left_join(as_tibble(tree),bs)

###Export text files for every label in partition and populate with properly formatted gene id's for automatic retreaval---------
cat ("========================================\nWriting Clade Lists\n")
joint_clade_list<-vector()
for (n in 1:(nrow(partition))) {
  (clade <- partition[n,]$label)
  (node <- partition[n,]$node)
  (tips <- offspring(x,node, tiponly = T, self_include = T))
  (tipnames <- unlist(strsplit(tips$label,"/",fixed = T))[2*1:nrow(tips)-1] %>% unique() )
  joint_clade_list <- rbind(joint_clade_list,tibble(Clade = paste0(clade,"_",length(tipnames)), Tips = tipnames))
  write_delim(x = as.data.frame(tipnames), path = paste0(dir,"/",clade, "_",length(tipnames),".txt"), delim = "\t",quote_escape = "double",append = F,col_names = F)
}
cat ("========================================\nWriting Joint Clade List\n")
write_delim(joint_clade_list,list,delim = "\t",col_names = FALSE)
cat ("========================================\nDone\n")
