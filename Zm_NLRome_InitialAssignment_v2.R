## ---------------------------
##
## Script name: Zm_NLRome_InitialAssignment.R
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
source("~/GitHub/Gm_NLRome/Phylo_functions.R")
library("parallel")
library("tidyverse")
library("ggtree")
library("treeio")
library("tidytree")
library("phangorn")

'%ni%' <- Negate('%in%')
setwd("~/Dropbox/NLRomes/Maize_NLRome/NLRome/")

########Import RAXML tree---------------
raxml <- read.raxml("RAxML_bipartitionsBranchLabels.ZMonly.100.Raxml.out")

########Annotate Internal Nodes and save Bootstrap values---------------
z <- mutate(as_tibble(raxml), label=if_else(is.na(label),paste0("Int",node),label))
bs <- z %>% select(label,bootstrap)

########Root midpoint -----------------------
tree <- midpoint(as.phylo(z))
x <- left_join(as_tibble(tree),bs)
x %>% print(n=4000)
as.phylo(x)

x %>% filter(label == "Int4000")
z %>% filter(label == "Int4000")

is.rooted(as.phylo(x))

###Get number of leaves per node------
N_tips<-vector(length = nrow(x))
for (i in 1:nrow(x)){
  nd <- x$node[[i]]
  N_tips[[i]] <- length(offspring(x,nd, tiponly = T,self_include = T)$node)
}
x <- mutate(x, N_tips = N_tips)
x %>% print(n=4000)

## Get a MRCA definition for every node in the tree ----------------
ncores <- 8
x <- x %>% mutate(mrca_id = mclapply(x$node, get_mrca_name,tree  = x,mc.cores = ncores)%>%unlist)
x %>% filter(!is.na(mrca_id))

## Find clades that are a good size and have best available support------
good_size_clades<-vector()
for (m in tree$tip.label){
  #m <- messy[5,]$node
  #m <- "7416_T436-R1/25-406"
  print(m)
  (ancestry <- ancestor(x,m) %>% filter(N_tips >14, N_tips <250))  ### modified given 27 ecotypes
  n <- ancestry %>% filter(bootstrap == max(ancestry$bootstrap))
  good_size_clades<-rbind(good_size_clades,n)
}
good_size_clades <- good_size_clades %>% distinct()
good_size_clades %>% arrange(bootstrap) %>% print(n=200)
good_size_clades %>% filter(bootstrap <70) %>% arrange(N_tips) %>% print(n=50) 

##check that the assignment is unique
a <- good_size_clades$node
anc_pool<-vector()
offs_pool<-vector()
for (nd in a){
  anc <- ancestor(x,nd)$node
  anc_pool <- c(anc_pool,anc)
  offs <- offspring(x,nd)$node
  offs_pool <- c(offs_pool,offs)
}
anc_pool <- unique(anc_pool)
anc_pool
offs_pool <- unique(offs_pool)
offs_pool

partition <- good_size_clades
partition
(partition <- x[a[which(a %ni% offs_pool)],])
partition %>% select(N_tips) %>% sum()
x %>% filter(is.na(bootstrap))
partition %>% arrange(bootstrap) %>% print(n=300)
partition %>% filter(bootstrap < 70) %>% arrange(N_tips) %>% print(n=50) ## 4 poor BS clades
partition %>% arrange(N_tips) %>% print(n=300)


###Where are the missing tips and why are they missing?
tips<-tree$tip.label
missing<-vector()
for (a in 1:length(tips)){ if (x[a,]$node %ni% offs_pool){missing<-rbind(missing,x[a,])}}
missing %>% print(n=nrow(missing)) ## nothing missing

# compl<-vector()
# for (n in missing$node){
#   (ancestry <- ancestor(x,n) %>% filter(N_tips <51))
#   c <- ancestry %>% filter(N_tips == max(ancestry$N_tips))
#   compl<-rbind(compl,c)
# }
# compl <- unique(compl)
# compl
# #compl <- missing
# partition <-rbind(partition,compl)
# 
# ##check that the assignment is unique
# a <- partition$node
# anc_pool<-vector()
# offs_pool<-vector()
# for (nd in a){
#   anc <- ancestor(x,nd)$node
#   anc_pool <- c(anc_pool,anc)
#   offs <- offspring(x,nd)$node
#   offs_pool <- c(offs_pool,offs)
# }
# anc_pool <- unique(anc_pool)
# offs_pool <- unique(offs_pool)
# 
# partition
# x[a[which(a %ni% offs_pool)],]
# partition <- x[a[which(a %ni% offs_pool)],]
# partition %>% select(N_tips) %>% sum()
# length(tips)
# partition %>% arrange(N_tips) %>% print(n=300)

ggplot(partition, aes(x=N_tips))+geom_histogram(binwidth = 1)
ggplot(partition, aes(x=N_tips))+geom_histogram(binwidth = 25)
ggplot(partition, aes(x=bootstrap))+geom_histogram(binwidth = 1)

###Export text files for every label in partition40 and populate with properly formated gene id's for automatic retreaval---------
for (n in 1:(nrow(partition))) {
  clade <- partition[n,]$label
  node <- partition[n,]$node
  tips <- offspring(x,node, tiponly = T, self_include = T)
  tipnames <- unlist(strsplit(tips$label,"/",fixed = T))[2*1:nrow(tips)-1]
  write_delim(x = as.data.frame(tipnames), path = paste0("./Autoclades_70/",clade, "_",length(tips$label),".txt"), delim = "\t",quote_escape = "double",append = F,col_names = F)
}
