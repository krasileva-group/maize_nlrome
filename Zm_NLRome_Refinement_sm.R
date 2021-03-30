## ---------------------------
##
## Script name: Zm_NLRome_Refinement_sm.R
##
## Purpose of script: Refine clades by splitting into smaller sub-clades
##
## Author: Daniil Prigozhin
##
## Date Created: 2021-03-29
##
## Copyright (c) Daniil Prigozhin, 2020
## Email: daniilprigozhin@lbl.gov
##
## ---------------------------
##
## Notes: This is a script where variable parameters are set on top including input/output directories. 
##        This allows the same script to be used in successive refinements.
##        Snakemake version.
##
## ---------------------------
##Installing Packages for alignment manipulation-----------
# install.packages("entropy")
# install.packages("tidyverse")
# if (!requireNamespace("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install("msa")
# BiocManager::install("odseq")
# BiocManager::install("ggtree")

#Loading libraries-----------------------------------------
require("tidyverse")
require("ggtree")
require("treeio")
require("msa")
require("entropy")
#library(odseq)

##########################################
## Snakemake parameters ------------------
MinGapFraction <- snakemake@params[["MinGapFraction"]]
MinGapBlockWidth <- snakemake@params[["MinGapBlockWidth"]]
hvSiteEntCutoff <-  snakemake@params[["hvSiteEntCutoff"]]
OutputDirectory <- dirname(snakemake@input[["afa"]][[1]])
cladestar <- snakemake@output[["cladestar"]]
ncores <- snakemake@threads
##########################################
## Parameters for first refinement-------
# 
# setwd("~/Dropbox/MaizeNLRome/NLRome/Autoclades_70/")
# MinGapFraction <- 0.9
# MinGapBlockWidth <- 1
# hvSiteEntCutoff <-  1.5
# OutputDirectory <- "../Autoclades_70_Refinement_1/"
# ##########################################

##########################################
# ### Parameters for second refinement-------
# 
# setwd("~/Dropbox/MaizeNLRome/NLRome/Autoclades_70_Refinement_1/")
# MinGapFraction <- 1
# MinGapBlockWidth <- 1
# hvSiteEntCutoff <-  1.5
# OutputDirectory <- "../Autoclades_70_Refinement_2/"
# #########################################

##########################################
# ### Parameters for third refinement-------
# 
# setwd("~/BioInf/Brachy/Protein/RAxML_237aa/Autoclades_70_Refinement_2/")
# MinGapFraction <- 1
# MinGapBlockWidth <- 1
# hvSiteEntCutoff <-  1.5
# OutputDirectory <- "../Autoclades_70_Refinement_3"
# ##########################################

##########################################
# # ### Parameters for fourth refinement-------
# 
# setwd("~/BioInf/Brachy/Protein/RAxML_237aa/Autoclades_70_Refinement_3/")
# MinGapFraction <- 1
# MinGapBlockWidth <- 1
# hvSiteEntCutoff <-  1.5
# OutputDirectory <- "../Autoclades_70_Refinement_4/"
##########################################

##########################################
# ### Parameters for fifth refinement-------
# 
# setwd("~/BioInf/Brachy/Protein/RAxML_237aa/Autoclades_70_Refinement_4/")
# MinGapFraction <- 1
# MinGapBlockWidth <- 1
# hvSiteEntCutoff <-  1.5
# OutputDirectory <- "../Autoclades_70_Refinement_5/"
# ##########################################

#if (!dir.exists(OutputDirectory)){dir.create(OutputDirectory)}

##Collect subdirectories and aligment files ending in *.best.fa
files <- snakemake@input[["afa"]] ### This now includes all subsequent subclades

##Calculate entropy with and without gap character, gather alignment stats---------------------
Alph_21 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","-")
Alph_20 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

EntropyNG <- vector("list",length(files))
Entropy <- vector("list",length(files))
stats<-vector()
cat("======================================\nReading .afa files\n")
## Entropy calculation
for (i in seq_along(files)) {
  ## Read alignment
  maa <- readAAMultipleAlignment(files[[i]])
  
  ## Extracting folder name
  folder <- str_split(files[[i]], "/")[[1]][3] %>% str_remove(".afa")
  cat(paste0(folder,"\n"))
  ## Filter out gappy columns
  
  if ("-" %in% rownames(consensusMatrix(maa))){
    autoMasked <- maskGaps(maa, min.fraction = MinGapFraction, min.block.width = MinGapBlockWidth) ##KEY FILTERING PARAMETERS
    MinAli <- as(autoMasked, "AAStringSet")
    }else{MinAli<-as(maa, "AAStringSet")}
  MinAli
  
  ## Calculating Consensus Matrix
  (Tidy_CM<-as_tibble(t(consensusMatrix(MinAli, baseOnly = T))))
  ## Compensating for consensus matrix not keeping full alphabet in output
  for (a in setdiff(Alph_21,colnames(Tidy_CM))){
    vec <- as_tibble(0*(1:nrow(Tidy_CM)))
    colnames(vec) <- paste(a)
    Tidy_CM <- as_tibble(cbind(Tidy_CM,vec))
  } 
  ##Selecting relevant columns
  (Tidy_CM_Gaps <- select(Tidy_CM,all_of(Alph_21)))
  (Tidy_CM_NoGaps <- select(Tidy_CM,all_of(Alph_20)))
  
  ##Entropy Calculation
  ent <- apply(Tidy_CM_Gaps, 1, entropy,unit="log2") %>% as_tibble()
  colnames(ent)<-paste0("Entropy_",folder)
  ent
  
  ##Entropy Calculation Ignoring Gaps
  entNG <- apply(Tidy_CM_NoGaps, 1, entropy,unit="log2") %>% as_tibble()
  colnames(entNG)<-paste0("EntropyNoGaps_",folder)
  entNG
  
  ##Save fraction of invariant positions with/without gaps and number of highly variable positions (without gaps)
  frbl <- length(which(ent == 0))/nrow(ent)
  frblng <- length(which(entNG == 0))/nrow(entNG)
  nHVsites <- length(which(entNG > hvSiteEntCutoff))                          ####KEY CUTOFF PARAMETER
  stats <- rbind(stats,c(folder,frbl,frblng,nrow(ent),nHVsites))  
  
  ##Save results of entropy calculation
  Entropy[[i]] <- ent
  EntropyNG[[i]] <- entNG
}

#EntropyNG

##Format stats tibble-------
colnames(stats)<-c("Clade","FractionZero","FractionZeroNG","Ali_Length","N_HV_Sites")
stats<-as_tibble(stats)
#stats
stats<-mutate(stats, Clade = Clade,
              FractionZero = as.numeric(FractionZero),
              FractionZeroNG = as.numeric(FractionZeroNG),
              N_HV_Sites = as.numeric(N_HV_Sites),
              Ali_Length = as.numeric(Ali_Length)
)
stats %>% arrange(FractionZeroNG) %>% print(n=500)
stats %>% filter(FractionZeroNG<0.9) %>% arrange(N_HV_Sites) %>% print(n=500)
write_delim(stats, path = paste0(OutputDirectory,"/","Stats.txt"), delim = "\t", col_names = T)
#stop()
##Plot distributions of invariant sites with/without gaps and of highly variable sites---------
# ggplot(stats, aes(FractionZeroNG))+geom_histogram()
# ggplot(stats, aes(N_HV_Sites))+geom_histogram()#+xlim(-1,50)
# 
# ggplot(stats, aes(x=FractionZeroNG, y=N_HV_Sites))+geom_point()#+ylim(0,100)
# ggplot(stats, aes(x=FractionZeroNG, y=FractionZero))+geom_point()+xlim(0,1)+ylim(0,1)

##Use cutoff value of 90% invariant sites (without gaps) for clades not to be split
# DoNotSplit <- stats %>% filter(FractionZeroNG > 0.9) %>% arrange(Clade) %>% print(n=nrow(stats))

# ToiTol <- stats %>% filter(FractionZeroNG < 0.9) %>% arrange(N_HV_Sites) %>% print(n=nrow(stats))
# ToiTol <- mutate(ToiTol, path = paste0(Clade,"/","RAxML_bipartitionsBranchLabels.",Clade,".First.out"))
# write_delim(ToiTol %>% select(path),"To_iTol.list",col_names = F, delim = " ")

##Plot Entropy figures for all Clades, save in ./Plots/----------------------
## Set plot hight in bits
cat("======================================\nPlotting entropy\n")
m <- max(2.2,na.omit(unlist(Entropy)))
## Make Plots directory if it is not there already
plots_dir <- paste0(OutputDirectory,"/","Plots")
if (!dir.exists(plots_dir)){dir.create(plots_dir)}
## Entropy plotting
if(1==0){for (l in seq_along(Entropy)){
  folder <- str_remove(colnames(Entropy[[l]]),"Entropy_")
  cat(folder)
  cat("\n")
  Ent <- as_tibble(cbind(1:nrow(Entropy[[l]]),Entropy[[l]],EntropyNG[[l]]))
  colnames(Ent)<-c("Position","Entropy", "EntropyNG")
  
  ggplot(Ent, aes(x = Position))+
    geom_line(aes(y = Entropy), color = "blue")+
    ylim(0,m)+
    ggtitle(paste(folder)) +
    xlab("Position") + 
    ylab("Shannon Entropy with Gaps")+
    theme_classic()
  ggsave(paste0(plots_dir,"/",folder,"_Entropy_Masked",".pdf"))
  
  ggplot(Ent, aes(x = Position))+
    geom_line(aes(y = EntropyNG), color = "red")+
    ylim(0,m)+
    ggtitle(paste(folder)) +
    xlab("Position") + 
    ylab("Shannon Entropy")+
    theme_classic()
  ggsave(paste0(plots_dir,"/",folder,"_Entropy_MaskedNG",".pdf"))
}}

#stop()

########################################
###Analyse clade trees------------------
###Import trees, calculate tree statistics, store all trees in BigTable object
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

trees <- snakemake@input[["tree"]]
Annotation<-read_csv(snakemake@input[["annotation"]], col_names = c("Assembly","Name"))
(Annotation<-mutate(Annotation, label = toupper(Name)))
adress <- list.files("~/Dropbox/NLRomes/Maize_NLRome/RAxML_tree_pbNB-ARC/Initial_Clades/",pattern = ".out")[[1]]

BigTable<-vector("list",length = length(trees))
for (i in seq_along(trees)){
  #i <- 1
  adress <- trees[[i]]  
  folder <- str_split(basename(adress),pattern = "[.]")[[1]][2]
  cat(folder)
  cat("\n")
  tree <- read.raxml(file = paste(adress))
  (table <- as_tibble(tree))
  (table<-mutate(table,Clade = folder, mrca_id =  mclapply(table$node, get_mrca_name,tree  = table,mc.cores = ncores)%>%unlist))
  (table <- left_join(table,Annotation, by = "label"))
  ###Get number of unique ecotypes left and right of every node
  ###Get number of duplicated ecotypes left and right of every node 
  x<-table
  alltips<-tree@phylo$tip.label
  CountTable<-vector()
  for (nd in x$node){
    #nd<-434
    (R_tips <- offspring(x,nd, tiponly=T, self_include=T))
    (L_tips <- x[which(alltips %ni% R_tips$label),])
    nrow(L_tips)
    if (nrow(L_tips) == 0){L_NE <- 0; L_ND <-0}else{
      (L_ecotypes <- L_tips$Assembly)
      L_NE <- length(unique(L_ecotypes))
      L_ND <- length(unique(L_ecotypes[which(duplicated(L_ecotypes))]))
    }
    (R_ecotypes <- R_tips$Assembly)
    R_NE <- length(unique(R_ecotypes))
    R_ND <- length(unique(R_ecotypes[which(duplicated(R_ecotypes))]))
    OE <- length(intersect(R_ecotypes,L_ecotypes))
    vec <- c(nd,R_NE,R_ND,nrow(R_tips),L_NE,L_ND,nrow(L_tips),OE)
    vec
    CountTable <- rbind(CountTable, vec)
  }
  CountTable
  ###Merge with tree data
  colnames(CountTable) <- c("node","R_Eco","R_DuplEco","R_tips","L_Eco","L_DuplEco","L_tips","N_OverlapEco")
  x <- left_join(x,as_tibble(CountTable), by = "node")
  ###Add current tree with node statistics to the BigTable object
  BigTable[[i]] <- x      
}

## Flatten the list of tables into one table 
Full_table<-vector()
for (i in seq_along(BigTable)){Full_table <- rbind(Full_table,BigTable[[i]])} ##This way runs faster than doing rbind within loop above
BigTable <- Full_table

## Check that the table is populated correctly
BigTable %>% filter(!is.na(label)) %>% select(label) %>% distinct() ###Number of genes in the clades under refinement
BigTable %>% select(label, bootstrap, Clade, Assembly, Name, R_tips, N_OverlapEco,mrca_id)

# ggplot(BigTable,aes(x=branch.length))+geom_histogram()#+ylim(-1,201)
# ggplot(BigTable,aes(x=R_Eco))+geom_histogram()#+ylim(-1,201)

export <- BigTable %>% filter(branch.length >0.1) %>% print(n=200) ##These branches should be checked
write_delim(export %>% select(Clade,parent,node,branch.length,label,bootstrap,mrca_id,N_OverlapEco) %>% arrange(Clade,branch.length),snakemake@output[["tsv"]])


#### Use below to export additional cut nodes:
# setwd("~/Dropbox/NLRomes/Maize_NLRome/RAxML_tree_pbNB-ARC/Initial_Clades/")
# export <- BigTable %>% filter(branch.length >0.06,bootstrap == 100,Clade == "Int5286_45") %>% print(n=200) ##These branches should be checked
# write_delim(export %>% select(Clade,parent,node,branch.length,label,bootstrap,mrca_id,N_OverlapEco) %>% arrange(Clade,branch.length),"Extra_nodes.txt")


ii<-5
sink(cladestar, append = F)
cat("DATASET_SYMBOL
SEPARATOR SPACE

#label is used in the legend table (can be changed later)
DATASET_LABEL CutNodes")
#cat(paste0("_",min_seq,"_",max_seq,"\n"))
cat("

#dataset color (can be changed later)
COLOR #ffff00
MAXIMUM_SIZE 10
BORDER_WIDTH 2
BORDER_COLOR #000000
DATA
")
for (ii in 1:length(export$label)){
  if(!is.na(export[ii,7])){
  cat(export[ii,7] %>% unlist())}else{cat(export[ii,4] %>% unlist())}
  if (export[ii,3]>0.3){cat(" 3 10 #00FF00 1 0.25\n")}else
    if (export[ii,3]>0.1){cat(" 3 10 #FFFF00 1 0.25\n")}else
    {cat(" 3 10 #ff0000 1 0.25\n")}
}
sink()




save.image(file = "SM.RData")
# (a <- BigTable %>% filter(branch.length >0.1,N_OverlapEco>13,bootstrap>98))
# BigTable %>% filter(branch.length >0.1,N_OverlapEco>1,bootstrap>98) %>% arrange(-N_OverlapEco)
# 
# Cut_Nodes_10 <- vector()
# Cut_Nodes_10 <- unique(rbind(Cut_Nodes_10,a)) %>% print(n=100)
# 
# Cut_Nodes_10 %>% filter(grepl("_219",Clade)) %>% arrange(-branch.length)%>% print (n=200)
# 
# Cut_Nodes_10 %>% select(node,branch.length,label,bootstrap,Clade,Name,N_OverlapEco)%>% print(n=300)
# ggplot(stats %>% filter(Clade %ni% Cut_Nodes_10$Clade), aes(x=FractionZeroNG))+geom_histogram()
# 
# 
# 
# AutoPick<-Cut_Nodes_10 %>% arrange(branch.length)
# 
# # ## Import Saved cut list: 
# # Cut_Nodes_10<-read_delim("CutListSplit1_70.txt",delim = " ",col_names = T)
# # Cut_Nodes_10<-read_delim("../Autoclades_70_Refinement_1/CutList_Refinement_1.txt",delim = " ",col_names = T)
# # Cut_Nodes_10<-read_delim("../Autoclades_70_Refinement_2/CutList_Refinement_2.txt",delim = " ",col_names = T)
# # Cut_Nodes_10<-read_delim("../Autoclades_70_Refinement_3/CutList_Refinement_3.txt",delim = " ",col_names = T)
# # Cut_Nodes_10 <- read_delim(paste0(OutputDirectory,"/CutList_Refinement_4_Int21623.txt"),delim = ' ')
# 
# 
# Cut_Nodes_10
# 
# ### Write new cut list table. Do Not Use below without redoing subsequent refinements
#  write_delim(Cut_Nodes_10, paste0(OutputDirectory,"/CutList_Refinement_1.txt"),quote_escape = "double")
# # write_delim(Cut_Nodes_10 ,paste0(OutputDirectory,"/CutList_Refinement_2.txt"),quote_escape = "double")
# # write_delim(Cut_Nodes_10 ,paste0(OutputDirectory,"/CutList_Refinement_3.txt"),quote_escape = "double")
# # write_delim(Cut_Nodes_10 ,paste0(OutputDirectory,"/CutList_Refinement_4.txt"),quote_escape = "double")
# # write_delim(Cut_Nodes_10, paste0(OutputDirectory,"/CutList_Refinement_1_Int21623.txt"),quote_escape = "double")
# # write_delim(Cut_Nodes_10, paste0(OutputDirectory,"/CutList_Refinement_2_Int21623.txt"),quote_escape = "double")
# # write_delim(Cut_Nodes_10, paste0(OutputDirectory,"/CutList_Refinement_3_Int21623.txt"),quote_escape = "double")
# # write_delim(Cut_Nodes_10, paste0(OutputDirectory,"/CutList_Refinement_4_Int21623.txt"),quote_escape = "double")
# # 
# ## The trees from raxml are unrooted, therefore the program traverses the trees, choosing the ""biggest"" node to split on
# ## Takes everything behind it as a clade, then takes individual leaves and walks them up to the first Cut_Nodes_10 node
# ## A new column "Split_1" with the name of that node appended to the original Clade name works as new ubclade name.
# ## Clade lists can then be exported based on that column. 
# ## Run over all the nodes in Big Table (the concatenated tree collection), and for each tip assign a new Split_1 clade
# ## For non-tip lines, give NA. For tips trees with no nodes in the cut list, give _NS_N to signify lack of further splitting.
# ## For other tips figure out if the current Tip is in the offspring of the Top Clade (one with the most tips)
# ## If no, assing the name CurClade_TopClade_L
# ## If yes, find the smallest ancestor in the SplitNodes
# ## Assign the name CurClade_SmallestAncestor_R
# 
# ## Append a column to BigTable that contains the TRUE values for nodes in the cut list.
# Cut_Nodes_10 <- mutate(Cut_Nodes_10, Split_Node_1 = T)
# BigTable <- Full_table
# BigTable <- left_join(BigTable, Cut_Nodes_10 %>% select(label,node,Clade,Split_Node_1)%>%distinct(), by = c("label","node","Clade"))
# BigTable %>% filter(Split_Node_1)
# Cut_Nodes_10 %>% filter(Split_Node_1)%>%select(Clade)
# 
# ## Generate Split_1 column with new clade assignment for every tip of every tree
# Split_1 <-vector(length = nrow(BigTable))
# for (i in 1:nrow(BigTable)){if (!is.na(BigTable[i,]$label)){                    #works
#   (Tip <- BigTable[i,])
#   (CurClade <- BigTable[i,]$Clade)
#   (TreeData <- filter(BigTable,Clade == CurClade))
#   (SplitNodes <- filter(TreeData,Split_Node_1))
#   if (nrow(SplitNodes)>0){
#     (TopClade <- SplitNodes %>% filter(R_tips == max(SplitNodes$R_tips)))
#     if (is.na(Tip$Split_Node_1)){
#       if (nrow(dplyr::intersect(ancestor(TreeData,Tip$node),SplitNodes))>0){
#         SplitAncestors <- dplyr::intersect(ancestor(TreeData,Tip$node),SplitNodes)
#         SmallestAncestor <- SplitAncestors %>% filter(R_tips == min(SplitAncestors$R_tips))
#         BestNode <- paste0(CurClade,'_',SmallestAncestor$node,'_R')
#       }else{BestNode <- paste0(CurClade,'_',TopClade$node,'_L')}
#     }else{BestNode <- paste0(CurClade,'_',Tip$node,'_R')}
#   }else{BestNode <- paste0(CurClade,"_NS_N")}
#   Split_1[[i]] <- BestNode 
# }else{Split_1[[i]] <- NA}}
# 
# as_tibble(Split_1) %>% print(n=3000)
# (BigTable_1 <- mutate(BigTable, Split_1 = as.factor(Split_1)))
# BigTable_1 %>% select(label, Split_1) %>% print(n=1000)
# ## Count numbers of Maize tips in each new clade
# CladeStat <- BigTable_1 %>% filter(!is.na(Split_1),!grepl('SORBI',label)) %>% group_by(Clade,Split_1) %>% summarise(n=n()) %>% print(n=300)
# CladeStat %>% ungroup() %>% select(Clade) %>% distinct()
# CladeStat %>% ungroup() %>% select(Split_1) %>% distinct()
# BigTable  %>% select(label) %>% distinct()
# write_delim(CladeStat, paste0(OutputDirectory,"/CladeStat.txt"), delim = "\t", append = F, col_names = T)
# 
# ggplot(CladeStat,aes(x=n))+geom_density()
# ggplot(CladeStat,aes(x=n))+geom_density()+xlim(0,100)
# ggplot(CladeStat,aes(x=n))+geom_histogram()
# 
# # ## First Refinement. Print out new clade assignment for every tip ###This needs to change to get names identical to clade names below i.e. end in number of genes in the clade
#  FirstAssignment <- BigTable_1 %>% filter(!is.na(Split_1)) %>% select(label,Split_1)
#  write_delim(x = FirstAssignment, path = paste0(OutputDirectory,"/firstassignment.txt") , delim = "\t",quote_escape = "double",append = F,col_names = T)
# 
# # ## Second Refinement. Print out new clade assignment for every tip
#  # SecondAssignment <- BigTable_1 %>% filter(!is.na(Split_1)) %>% select(label,Split_1)
#  # write_delim(x = SecondAssignment, path =  paste0(OutputDirectory,"/secondassignment.txt") , delim = "\t",quote_escape = "double",append = F,col_names = T)
# 
# # ## Third Refinement. Print out new clade assignment for every tip
# # ThirdAssignment <- BigTable_1 %>% filter(!is.na(Split_1)) %>% select(label,Split_1)
# # write_delim(x = ThirdAssignment, path =  paste0(OutputDirectory,"/thirdassignment.txt") , delim = "\t",quote_escape = "double",append = F,col_names = T)
# 
# # ## Fourth Refinement. Print out new clade assignment for every tip
#  # Assignment <- BigTable_1 %>% filter(!is.na(Split_1)) %>% select(label,Split_1)
#  # write_delim(x = Assignment, path =  paste0(OutputDirectory,"/fourthassignment.txt") , delim = "\t",quote_escape = "double",append = F,col_names = T)
# 
# ###Find duplicate tips
# Duplicates <- BigTable_1 %>% group_by(label) %>%summarise(n=n()) %>% filter(n>1)
# Duplicates <- Duplicates %>% filter(!is.na(label))
# Duplicates
# BigTable_1 %>% filter(label %in% Duplicates$label) %>% group_by(Split_1) %>% summarise(n=n())
# BigTable_1 %>% filter(label %in% Duplicates$label) %>% arrange(label) %>% print(n=50) #%>% group_by(Clade) %>% summarise(n=n())
# 
# BigTable %>% filter(!is.na(Clade)) %>% select(Clade) %>% distinct()
# 
# Duplicates
# 
# ###Write Clade Lists - DO NOT RUN without repeating refinement
# for (n in 1:(nrow(CladeStat))) {
#   clade <- CladeStat[n,]$Split_1
#   tips <- BigTable_1 %>% filter(Split_1 == clade)
#   tipnames <- tips$label
#   write_delim(x = as.data.frame(tipnames), path = paste0(OutputDirectory,"/",clade, "_",length(tips$label),".txt"), delim = "\t",quote_escape = "double",append = F,col_names = F)
# }
