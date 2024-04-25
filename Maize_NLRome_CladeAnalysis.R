## ---------------------------
##
## Script name: Maize_NLRome_CladeAnalysis.R
##
## Purpose of script: Synthesize info on clade assignment, ecotypes, alleles, etc. in one "Common" table.
##
## Author: Daniil Prigozhin
##
## Date Created: 2019-08-20
##
## Copyright (c) Daniil Prigozhin, 2020
## Email: daniilprigozhin@lbl.gov
##
## ---------------------------
##
## Notes: Last tested under R version 3.6.3
##   
##
## ---------------------------

setwd("~/Dropbox/NLRomes/Maize_NLRome/")

library(tidyverse)
########################################################################
### Build a Common table of clade assignments for every NLR ------------
########################################################################

'%ni%' <- Negate('%in%')
(zeroth <-read_delim("./RAxML_tree_pbNB-ARC/Initial_Clades.pbNB-ARC.list",col_names = c("Clade_0","Gene"), delim = "\t"))
(first <-read_delim("./RAxML_tree_pbNB-ARC/Refinement_1.pbNB-ARC.Clades.list",col_names =T, delim = "\t"))
first <- first %>% mutate(Clade_1=Split,Gene=label)%>%select(-Split,-label)
(second <-read_delim("./RAxML_tree_pbNB-ARC/Refinement_2.pbNB-ARC.Clades.list",col_names = T, delim = "\t"))
second <- second %>% mutate(Clade_2=Split,Gene=label)%>%select(-Split,-label)
(third <-read_delim("./RAxML_tree_pbNB-ARC/Refinement_3.pbNB-ARC.Clades.list",col_names = T, delim = "\t"))
third <- third %>% mutate(Clade_3=Split,Gene=label)%>%select(-Split,-label)

####Remove duplicate assignment
#(first <- first %>% filter(Clade_1 != "Int9836_17_15_L_16"))

###JOIN DATA TOGETHER############
Common<-left_join(zeroth,first,by = "Gene")
Common<-left_join(Common,second,by = "Gene")
Common<-left_join(Common,third,by = "Gene")
Common %>% print(n=444444)
Common<-Common[c(2,1,3,4,5)] ###Rearrange Columns

###FIND Multiple assignments for individual hits###########################
Common %>% select(Gene) %>% unique()

Common  %>% group_by(Gene) %>% filter(n()>1)  %>% summarise(N=n())%>% arrange(Gene) %>% print(n=300)
Common %>% group_by(Gene) %>% filter(n()>1)  %>% arrange(Gene) %>% print(n=300)


######MAKE a Column that uniquely establishes final assigned clade######
CladeA<-vector()
for (l in c(1:nrow(Common))){
  b<-Common[l,]
  a<-paste(b[2], sep ="")
  if (!is.na(b[3])){a<-paste(b[3], sep ="")}
  if (!is.na(b[4])){a<-paste(b[4], sep ="")}
  if (!is.na(b[5])){a<-paste(b[5], sep ="")}
  CladeA<-append(CladeA, a, after = length(CladeA))
}


Common<-mutate(Common, Clade = CladeA)
Common
Common %>% select(Gene, Clade) %>% unique() %>% group_by(Gene) %>% filter(n()>1) %>% print(n=39) 
Common %>% select(Gene, Clade) %>% unique() 

####Find problem Clades#########
Common %>% select(Gene, Clade) %>% unique() %>% group_by(Gene) %>% filter(n()>1) %>% ungroup() %>% select(Clade) %>% unique()
Common <-Common %>% filter(Clade != "Int4916_25_31_L_7") ## removes duplicate clade


#####Make a column that contains Ecotype ID only######



meta <- read_delim("NLR_Map.csv", delim = ",", col_names = c("Ecotype","Gene"))
meta <- meta %>% mutate(Gene = Gene %>% toupper())
Common <- left_join(Common,meta)
Common
###GET allele number##############
(Common <-mutate(Common, Allele = 1))

# #####Find non-unique gene assignments#########
# FixList <- Common %>% select(Gene,Clade) %>% distinct() %>% group_by(Gene) %>% filter(n()>1) %>% select(Gene) %>% distinct()
# FixList
# Common %>% filter(Gene %in% FixList$Gene) %>% arrange(Gene)
# ####Need to think about how to fix these. Just 1 at the moment########

###FIND Singletons#########
SingletonClades <- Common %>% group_by(Clade) %>% summarise(n=n()) %>% filter(n<2) %>% select(Clade)
Common %>% filter(Gene %in% (Common %>% filter(Clade %in% SingletonClades$Clade))$Gene) %>% arrange(Gene) %>% print(n=111111)
##Currently 26 singletons.

#### Create a Pointer table to link individual genes to assemblies. Add column to Common.
dirs <- dir("RAxML_tree_pbNB-ARC", pattern = "_", full.names = T, recursive = F)
files <- list.files(path = dirs, pattern = ".afa$",full.names = T)
files
(clades <- unlist(strsplit(files, "/"))[3*1:length(files)] %>% str_remove(".afa"))
Pointer <- tibble(Clade = clades, File = files)
Pointer %>% select(Clade) %>% distinct()
Common <- left_join(Common,Pointer, by = "Clade")
Common %>% filter(is.na(File))%>% select(Clade) %>% unique() %>% print(n=111111) #27 singleton clades

###############################################
### Look at clade properties ------------------
###############################################

###CLADE Sizes####
Frequency_table <- Common %>% filter(Ecotype != "Athaliana", Allele ==1)%>% group_by(Clade) %>% summarise(n=n()) %>% arrange(n)
Frequency_table <- Common %>% group_by(Clade) %>% summarise(n=n()) %>% arrange(n) %>% print(n=nrow(Common))
ggplot(Frequency_table, aes (x=n))+geom_density()+theme_classic()
ggplot(Frequency_table, aes (x=n))+geom_density()+theme_classic()+xlim(0,85)
ggplot(Frequency_table, aes (x=n))+geom_bar()+theme_classic()#+xlim(0,85)
Frequency_table %>% select(n) %>% sum() ###3496
###Currently 190 clades based on the 3496-tip tree
Frequency_table <- Common %>% group_by(Clade) %>% summarise(n=n()) %>% arrange(n)
Frequency_table %>% filter(n<20) %>% select(n) %>% sum()
739/3496
###Clade Characteristics###
AllClades<-levels(as.factor(Common$Clade))
Common$Clade
AllClades
AllEcotypes<-levels(Common$Ecotype)
N_Ecotype<-vector()
NDupl_Ecotype<-vector()
for (l in AllClades){
  #l <- "2_AT4G16920.1"
  (a <- Common %>% filter(Ecotype != "Athaliana" & Allele == 1)%>%filter(Clade == l) %>% select (Gene,Clade,Ecotype) %>% distinct() %>% group_by(Ecotype) %>% filter(n()>1) %>% distinct(Ecotype) %>% nrow())
  (b <- Common %>% filter(Ecotype != "Athaliana" & Allele == 1)%>%filter(Clade == l) %>% select(Ecotype) %>% distinct()%>% nrow())
  N_Ecotype<-rbind(N_Ecotype, c(l,b))
  NDupl_Ecotype<-rbind(NDupl_Ecotype, c(l,a))
}
N_Ecotype
colnames(N_Ecotype)<-c("Clade","N_Ecotype")
colnames(NDupl_Ecotype)<-c("Clade","NDupl_Ecotype")
(N_Ecotype<-as_tibble(N_Ecotype))
(NDupl_Ecotype<-as_tibble(NDupl_Ecotype))

EcoTibble<-left_join(N_Ecotype,NDupl_Ecotype,by = "Clade")
EcoTibble$N_Ecotype <- as.numeric(EcoTibble$N_Ecotype)
EcoTibble$NDupl_Ecotype <- as.numeric(EcoTibble$NDupl_Ecotype)
EcoTibble

EcoTibble %>% arrange(-NDupl_Ecotype)%>% print(n=20000)

pd <- position_dodge(width = 1)
ggplot(EcoTibble, aes(x=N_Ecotype, y=NDupl_Ecotype))+geom_point(position = pd)

EcoTibble %>% filter(NDupl_Ecotype >0) %>%arrange(N_Ecotype) %>% print(n=200)


##################################################
### Calculate entropy of the clade alignments  ---

### Now it is possible to create a list of genes with entropy strings 
### with positions corresponding to the residues of the gene (i.e. remove gaps at the level of the gene)
### Can do this while cycling over Clades to save on reading
library(msa)
library(entropy)
hvSiteEntCutoff <- 1.5
MinGapFraction <- 1 
MinGapBlockWidth <- 1
Alph_21 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","-")
Alph_20 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
files<-levels(as.factor(Common$File))
EntropyNG <- vector("list",length(files))
Entropy <- vector("list",length(files))
stats<-vector()
## Entropy calculation
i <- 1
for (i in seq_along(files)) {
  ## Read alignment
  
  maa <- readAAMultipleAlignment(files[[i]])
  
  ## Extracting folder name
  str <- str_split(files[[i]], "/")[[1]]
  folder <- str[[length(str)]]%>%str_remove(".afa")
  print(paste("Looking at clade ",folder,"..."))
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
  (Tidy_CM_Gaps <- select(Tidy_CM,Alph_21))
  (Tidy_CM_NoGaps <- select(Tidy_CM,Alph_20))
  
  ##Entropy Calculation
  ent <- apply(Tidy_CM_Gaps, 1, entropy,unit="log2") %>% as_tibble()
  colnames(ent)<-paste0("Entropy_",folder)
  
  
  ##Entropy Calculation Ignoring Gaps
  entNG <- apply(Tidy_CM_NoGaps, 1, entropy,unit="log2") %>% as_tibble()
  colnames(entNG)<-paste0("EntropyNoGaps_",folder)
  
  
  ##Save fraction of invariant positions with/without gaps and number of highly variable positions (without gaps)
  frbl <- length(which(ent == 0))/nrow(ent)
  frblng <- length(which(entNG == 0))/nrow(entNG)
  nHVsites <- length(which(entNG > hvSiteEntCutoff))                          ####KEY CUTOFF PARAMETER
  stats <- rbind(stats,c(folder,frbl,frblng,nrow(ent),nHVsites))  
  
  ##Save results of entropy calculation
  Entropy[[i]] <- ent
  EntropyNG[[i]] <- entNG
}

##Format stats tibble-------
colnames(stats)<-c("Clade","FractionZero","FractionZeroNG","Ali_Length","N_HV_Sites")
stats<-as_tibble(stats)
stats
stats<-mutate(stats, Clade = Clade,
              FractionZero = as.numeric(FractionZero),
              FractionZeroNG = as.numeric(FractionZeroNG),
              N_HV_Sites = as.numeric(N_HV_Sites),
              Ali_Length = as.numeric(Ali_Length)
)
stats %>% arrange(-N_HV_Sites) %>% print(n=500)

N_Tips <- Common %>% group_by(Clade)%>%select(Gene)%>%summarise(N=n())
EcoTibble <- left_join(EcoTibble,N_Tips)
EcoTibble <- left_join(EcoTibble,stats)
HV_Clades <- stats %>% filter(N_HV_Sites>20)
HV_Clades
Common <- mutate(Common, HV = ifelse(Common$Clade %in% HV_Clades$Clade,1,0))
Common
#################################
### Export common gene table ----
#################################
write_delim(Common,"Maize_NLRome_GeneTable.txt",delim = "\t")
Common <- read_delim("Maize_NLRome_GeneTable.txt",delim = "\t")

### Plot highly variable gene per genotype-----
hvFreq <- Common %>% 
  filter(HV == 1, Allele == 1)%>%
  group_by(Ecotype)%>%
  summarise(n=n())%>% 
  arrange(n) %>%
  print(n=70)
ggplot(hvFreq,aes(x=n))+geom_histogram()
hvFreq %>% select(n) %>% unlist %>% median
### 4 of the original clades contain hvNLRs
Common %>% filter(HV == 1, Allele == 1) %>% select(Clade_0) %>% unique() %>% arrange()

# ### Chromosome 11 has 3 different hvNLR clades on it
# Common %>% filter(HV == 1, Allele == 1) %>% select(Chr,Clade_0) %>%distinct()%>% group_by(Chr) %>% count() %>% arrange()

## Final clade assignment has 190 clades (26 singletons)
Common %>% select(Clade) %>% unique()

# ## Export lists of genes of interest-----
# AthaHV<-Common %>% filter(Clade %in% HV_Clades$Clade, Allele == 1,Ecotype == "Athaliana")%>%select(Gene)%>%unique()%>% arrange(Gene)%>%print(n=50)
# # write_delim(AthaHV,"AthaHV.txt", col_names = F, delim = "\t", na = "")
# AthaHVPlus<-Common %>% filter(Clade %in% HV_Clades$Clade, Allele != 1,Ecotype == "Athaliana")%>%select(Gene)%>%unique()%>% arrange(Gene)%>%print(n=50)
# # write_delim(AthaHVPlus,"AthaHVPlus.txt", col_names = F, delim = "\t", na = "")

# ## Map common names to hvNLRs ----
# CommonNames <- read_delim("AthaKnownGenes.txt",delim = "\t", col_names = c("Gene","CommonName"))
# CommonNames %>% print(n=100)
# AthaHV <- left_join(AthaHV,CommonNames)
# AthaHV %>% distinct() %>% arrange(Gene) %>% print(n=50)

## Export a table of clade representation in every ecotype--------
CladeEco_long <- Common %>% select("Gene","Clade","Ecotype") %>% unique() %>% group_by(Clade,Ecotype) %>% summarise(N=n()) 
CladeEco_wide <- pivot_wider(CladeEco_long,names_from = Ecotype,values_from = N, values_fill = list (N = 0))  %>% ungroup()
write_delim(CladeEco_wide, "CladeEco.txt", delim = "\t")
getwd()

##############################################################
### Produce hvNLR iTOL annotation file---
##############################################################
Common
genreg <- read_delim("HMM_search_pbNB-ARC/all_samples.pbNB-ARC.allgene.list",delim = " ",col_names = "GeneReg")
genreg <- genreg %>% mutate(Gene = GeneReg %>% str_remove("/.*"))
col_hv <- tibble(HV = c(0,1),Color = c("#ffffff","#fdae61"))
export <- left_join(genreg, Common%>%select(Gene,HV))%>%left_join(col_hv)%>%select(GeneReg,Color,HV)

sink("Annotation/HV_String.iTOL.txt", append = F)
cat("DATASET_COLORSTRIP
SEPARATOR SPACE
DATASET_LABEL hvNLRs
COLOR #ff0000
COLOR_BRANCHES 0
DATA
")
for (ii in seq_along(export$GeneReg)){
  cat(paste(export[ii,],sep = " "))
  cat("\n")
}
sink()


Common %>% filter(HV == 1) %>% select(Gene)%>% distinct() %>% write_delim("Maize_hvNLR_List.txt",col_names = F)
