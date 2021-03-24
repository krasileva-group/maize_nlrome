## ---------------------------
##
## Script name: Zm_NLRome_InitialAssignment_v2.R
##
## Purpose of script:
##
## Author: Daniil Prigozhin
##
## Date Created: 2020-10-03
##
## Copyright (c) Daniil Prigozhin, 2020
## Email: daniilprigozhin@lbl.gov
##
## ---------------------------
##
## Notes: Take one tree (unrooted) from RAxML, min_tips and max_tips parameters and produce an annotation for iTOL with 
##        the proposed split. Extra points for mouse over with clade info, color coded bootstrap values.
##
## ---------------------------
## load packages

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("ggtree")
#BiocManager::install("treeio")
#BiocManager::install("tidytree")
# install.packages("tidyverse")
# install.packages("phangorn")
#source("~/BioInf/GitHub/Gm_NLRome/Phylo_functions.R") #macbook
#source("~/GitHub/Gm_NLRome/Phylo_functions.R")        #imac
require("parallel")
require("tidyverse")
require("ggtree")
require("treeio")
require("tidytree")
require("phangorn")

###### Useful functions -------------
'%ni%' <- Negate('%in%')

get_mrca_name <- function(tree,node){
  require(treeio)
  if (!inherits(tree,"tbl_tree")) stop("Agrument must be a 'tbl_tree' object, use as_tibble() on treedata or phylo objects before passing to this function")
  if (isTip(tree,node)){return(NA)}else{
    tip_1 <- offspring(tree,child(tree,node)[1,]$node,tiponly = TRUE, self_include = T)[1,] 
    tip_2 <- offspring(tree,child(tree,node)[2,]$node,tiponly = TRUE, self_include = T)[1,] 
    return(paste0(tip_1$label,"|", tip_2$label))
  }
}

count_tips <- function(tree,node){
  require(treeio)
  if (!inherits(tree,"tbl_tree")) stop("Agrument must be a 'tbl_tree' object, use as_tibble() on treedata or phylo objects before passing to this function")
  return(length(offspring(tree,node, tiponly = T,self_include = T)$node))
}

######## Import snakemake pointers ---------------
raxml <- read.raxml(snakemake@input[["tree"]])
min_seq <- snakemake@params[["min_seq"]]
max_seq <- snakemake@params[["max_seq"]]
cladestar <- snakemake@output[["cladestar"]]
tsv <- snakemake@output[["tsv"]]
cladestrip <- snakemake@output[["cladestrip"]]
ncores <- snakemake@threads

########Annotate Internal Nodes and save Bootstrap values---------------
cat ("========================================\nReading tree\n")
z <- mutate(as_tibble(raxml), label=if_else(is.na(label),paste0("Int",node),label))
bs <- z %>% select(label,bootstrap)
#write.tree(as.phylo(z),"ZMonly.100.ztree.txt")
########Root midpoint -----------------------
cat ("========================================\nRooting at midpoint\n")

tree <- midpoint(as.phylo(z))
x <- left_join(as_tibble(tree),bs)

###Get number of leaves per node------
cat ("========================================\nCounting leaves per node\n")

x <- x %>% mutate(N_tips = mclapply(node , tree = x, count_tips, mc.cores = ncores)%>%unlist)

## Get a MRCA definition for every node in the tree ----------------
cat ("========================================\nGetting node names by MRCA method\n")

x <- x %>% mutate(mrca_id = mclapply(x$node, get_mrca_name,tree  = x,mc.cores = ncores)%>%unlist)


## Find clades that are a good size and have best available support------
cat ("========================================\nFinding best scoring clades given min/max parameters\n")

good_size_clades<-vector()
for (m in tree$tip.label){
  (ancestry <- ancestor(x,m) %>% filter(N_tips > min_seq, N_tips < max_seq))  
  n <- ancestry %>% filter(bootstrap == max(ancestry$bootstrap))
  good_size_clades<-rbind(good_size_clades,n)
}
good_size_clades <- good_size_clades %>% distinct()

##check that the assignment is unique
a <- good_size_clades$node
offs_pool<-vector()
for (nd in a){
  offs <- offspring(x,nd)$node
  offs_pool <- c(offs_pool,offs)
}
offs_pool <- unique(offs_pool)
partition <- x[a[which(a %ni% offs_pool)],]

cat ("========================================\nFound clades for these many tips:\n")

partition %>% select(N_tips) %>% sum()

cat ("========================================\nTotal number tips:\n")

x %>% filter(is.na(bootstrap),!is.na(label)) %>% select(N_tips) %>% sum()

###Where are the missing tips and why are they missing?
tips<-tree$tip.label
missing<-vector()
for (a in 1:length(tips)){ if (x[a,]$node %ni% offs_pool){missing<-rbind(missing,x[a,])}}
if (!is_empty(missing)){
          compl<-vector()
          for (n in missing$node){
            (ancestry <- ancestor(x,n) %>% filter(N_tips <= min_seq))
            c <- ancestry %>% filter(N_tips == max(ancestry$N_tips))
            compl<-rbind(compl,c)
          }
          compl <- unique(compl)
          partition <-rbind(partition,compl)
          
          ##check that the assignment is unique
          a <- partition$node
          offs_pool<-vector()
          for (nd in a){
            offs <- offspring(x,nd)$node
            offs_pool <- c(offs_pool,offs)
          }
          offs_pool <- unique(offs_pool)
          partition <- x[a[which(a %ni% offs_pool)],]
          cat ("========================================\nAdding smaller clades for missing tips. Now covering these many tips:\n")
          partition %>% select(N_tips) %>% sum()
          cat ("========================================\nTotal number of tips:\n")
          length(tips)
}

cat ("========================================\nCurrent partition:\n")
partition %>% arrange(bootstrap) %>% print(n=300)
cat ("========================================\nClades below 70 bs value, inspect these in iTOL:\n")
partition %>% filter(bootstrap < 70) %>% arrange(N_tips) %>% print(n=50) ## 4 poor BS clades

# ggplot(partition, aes(x=N_tips))+geom_histogram(binwidth = 1)
# ggplot(partition, aes(x=N_tips))+geom_histogram(binwidth = 25)
# ggplot(partition, aes(x=bootstrap))+geom_histogram(binwidth = 1)

### Export text files for every label in partition40 and populate with properly formatted gene id's for automatic retrieval---------
# for (n in 1:(nrow(partition))) {
#   clade <- partition[n,]$label
#   node <- partition[n,]$node
#   tips <- offspring(x,node, tiponly = T, self_include = T)
#   tipnames <- unlist(strsplit(tips$label,"/",fixed = T))[2*1:nrow(tips)-1]
#   write_delim(x = as.data.frame(tipnames), path = paste0("./Autoclades_70/",clade, "_",length(tips$label),".txt"), delim = "\t",quote_escape = "double",append = F,col_names = F)
# }

########################################################
###Export partition, gene lists for initial clades -----
########################################################
write_delim(partition,tsv, delim = "\t")

# ###Export text files for every label in partition and populate with properly formatted gene id's for automatic retreaval---------
# for (n in 1:(nrow(partition))) {
#   clade <- partition[n,]$label
#   node <- partition[n,]$node
#   tips <- offspring(x,node, tiponly = T, self_include = T)
#   tipnames <- unlist(strsplit(tips$label,"/",fixed = T))[2*1:nrow(tips)-1]
#   write_delim(x = as.data.frame(tipnames), path = paste0("./Autoclades_70/",clade, "_",length(tips$label),".txt"), delim = "\t",quote_escape = "double",append = F,col_names = F)
# }

########################################################
###Export annotations for iTOL------
########################################################
### Initial Clades as Color Strips
# At this point nodes have names, which should simplify things

colourCount = length(partition$label)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
partition %>% select(mrca_id,label) %>% mutate(color = sample(getPalette(colourCount)), Annotation = label,label=mrca_id) %>%select(-mrca_id) ->export

sink(cladestrip,append = F)
cat(
  "DATASET_COLORSTRIP
SEPARATOR SPACE
DATASET_LABEL InitialClades")
cat(paste0("_",min_seq,"_",max_seq,"\n"))
cat("COLOR #ff0000
COLOR_BRANCHES 1
BORDER_WIDTH 0
BORDER_COLOR #0000ff
SHOW_INTERNAL 1
DATA
")
for (ii in 1:length(export$label)){
  cat(export[ii,] %>% unlist())
  cat("\n")
}
sink()
#getwd()

## Reference genes as text labels
# getwd()
# Names <- read_delim("Zea_Known_genes.txt", delim = "\t", col_names = c("ID","Name","Description"))
# Names <- Names %>% mutate(ID = toupper(ID))
# x %>% mutate(ID = str_remove(label,"_.*$")) %>% left_join(Names) %>% filter(!is.na(Name))%>%
#   select(label,Name)->export_names
# 
# sink("iTOL.NamedGenes.txt",append = F)
# cat(
#   "DATASET_TEXT
# SEPARATOR SPACE
# DATASET_LABEL NamedGenes
# COLOR #ff0000
# DATA
# ")
# for (ii in 1:length(export_names$label)){
#   cat(export_names[ii,] %>% unlist())
#   cat(" 1 #000000 bold-italic 1\n")
# }
# sink()

## Cut nodes as node symbols

sink(cladestar, append = F)
cat("DATASET_SYMBOL
SEPARATOR SPACE

#label is used in the legend table (can be changed later)
DATASET_LABEL CutNodes")
cat(paste0("_",min_seq,"_",max_seq,"\n"))
cat("

#dataset color (can be changed later)
COLOR #ffff00
MAXIMUM_SIZE 10
BORDER_WIDTH 2
BORDER_COLOR #000000
DATA
")
for (ii in 1:length(export$label)){
  cat(export[ii,1] %>% unlist())
  if (partition[ii,5]>90){cat(" 3 10 #00FF00 1 0.25\n")}else
    if (partition[ii,5]>70){cat(" 3 10 #FFFF00 1 0.25\n")}else
    {cat(" 3 10 #ff0000 1 0.25\n")}
}
sink()

