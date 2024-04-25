
package_list<-c("entropy","dplyr","msa","tidyverse")

load_pack <- function(x){
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE )
      require( i , character.only = TRUE , quietly = T)
    }
  }
}
load_pack(package_list)
#########################################################
### What % of Maize RLPs, RLKs and NLRs are highly variable ---
#########################################################
setwd("~/Dropbox/NLRomes/Maize_NLRome/")
Maize_Common <- read_delim("Maize_NLRome_GeneTable.txt",delim = "\t")
Maize_Common %>% select(Gene,HV)%>%distinct() %>%group_by(HV) %>% count
346/(346+2905)
## 10% of NLRs in hvNLR clades

##########################################################
### Find All Uniprot Genes in NLR/RLK/RLP hv Clades ------
##########################################################

Uni_Wm <- read_delim("~/Dropbox/NLRomes/Maize_NLRome/Proteomes/B73_AF.txt",delim = "\t",col_names = F)
Uni_Key <- Uni_Wm[,1:2] %>%transmute(Gene = X1, Uniprot = X2)
Maize_Common <- Maize_Common %>% left_join(Uni_Key) 
Maize_Common %>% filter(Ecotype == "b73") %>% filter(is.na(Uniprot)) %>%group_by(HV)
### there are 3 of 135 NLRs that do not have a 1-1 match in uniprot
Maize_Common %>% filter(HV ==1, Ecotype == "b73") %>% select(Clade,Uniprot) %>%arrange(Clade) %>%print(n=1000)
### All hvNLR clades have one or more 1-1 mapped uniprot representatives
Maize_Common %>% filter(HV ==1, Ecotype == "b73", !is.na(Uniprot)) %>% select(Clade,Uniprot) %>%arrange(Clade) %>%print(n=1000)
Maize_Common %>% filter(HV ==1, Ecotype == "b73", !is.na(Uniprot)) %>% select(Clade) %>% distinct() %>% nrow
### Have 9 Uniprot proteins representing 4 refined clades of 4 total refined hvNLR clades (from 4 initial clades)
Maize_Common %>% filter(HV ==1) %>% select(Clade_0) %>% distinct() %>% nrow
(NLR_models <- Maize_Common %>% filter(HV ==1, !is.na(Uniprot)) %>% select(Clade,Gene,Uniprot,File))


get_ent <- function(ali, refseq){
  gene <- refseq
  CN <- refseq
  maa <- ali
  ### Check that the protein name matches one alignment key ------------------
  lm <- length(grep(pattern = gene, x=rownames(maa)))
  if (lm == 1) {
    cat ("Found Reference Sequence\n")}else 
      if (lm ==0) {
        stop("Protein name not found in alignment", call.=FALSE)}else
          if (lm > 1) {
            stop("More than one protein matched the name provided", call.=FALSE)}else{stop("Error at protein name", call.=FALSE)}
  ## Masking columns by reference gene ----------------
  Alph_21 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","-")
  Alph_20 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  RefGene <- gene
  RefSeq <- as(maa, "AAStringSet")[grep(pattern = gene, x=rownames(maa))]
  GapMask<- NULL
  for (i in 1:width(RefSeq)){
    c<-as.vector(RefSeq[[1]][i]) %in% c("-")
    GapMask<-append(GapMask,c,length(GapMask))
  }
  colmask(maa) <- IRanges(GapMask)
  #Retrieving the non-masked subset ------------------
  RefAli <- as(maa, "AAStringSet")
  RefLen <- width(RefAli[1])
  ## Calculating Consensus Matrix -------------------
  Tidy_CM<-as_tibble(t(consensusMatrix(RefAli, baseOnly = T)))
  ## Compensating for consensus matrix not keeping full alphabet in output
  for (a in setdiff(Alph_21,colnames(Tidy_CM))){
    vec <- as_tibble(0*(1:nrow(Tidy_CM)))
    colnames(vec) <- paste(a)
    Tidy_CM <- as_tibble(cbind(Tidy_CM,vec))
  } 
  ##Selecting relevant columns
  Tidy_CM_NoGaps <- select(Tidy_CM,all_of(Alph_20))
  ##Entropy Calculation Ignoring Gaps ----------------------
  entNG <- apply(Tidy_CM_NoGaps, 1, entropy,unit="log2") %>% as_tibble()
  colnames(entNG)<-paste0("EntropyNoGaps_",gene)
  return(entNG)
}

print_ent <- function(entNG,file){
  #Output entropy results to file --------------------------
  sink(file = file,append = F)
  cat("attribute: shannonEntropy\n")
  cat("match mode: 1-to-1\n")
  cat("recipient: residues\n")
  for (ii in seq_along(entNG[[1]])){
    cat("\t")
    cat(paste0(":",ii))
    cat("\t")
    cat(sprintf("%.5f", entNG[[1]][ii]))
    cat("\n")
  }
  sink()   
}

pdf_ent <- function(entNG,file,refseq){
  Ent <- as_tibble(cbind(1:nrow(entNG),entNG))
  colnames(Ent)<-c("Position","Entropy")
  ggplot(Ent, aes(x = Position))+
    geom_line(aes(y = Entropy), color = "blue")+
    ylim(0,3)+
    ggtitle(paste(refseq)) +
    xlab("Position") + 
    ylab("Shannon Entropy")+
    theme_classic()
  ggsave(file)
}

#### Running NLR prep steps ----
setwd("~/Dropbox/NLRomes/Maize_NLRome/")
getwd()
if (!dir.exists("AF_NLR")) {dir.create("AF_NLR")}
prefix <- "AF_NLR/"    

Maize_Uniprot <- readAAStringSet("Proteomes/All_Maize_Uniprot.fasta")

for(clade in NLR_models$Clade%>%unique()){
  print(clade)
  if (!dir.exists(paste0(prefix,clade))) {dir.create(paste0(prefix,clade))}
  clade_dir <- paste0(prefix,clade,"/")  
  Uniprot_list <- NLR_models%>%filter(Clade ==clade) %>%select(Uniprot,File)%>%distinct()
  for (unip in Uniprot_list$Uniprot){
    system(paste0("wget -O ",clade_dir,unip,".pdb https://alphafold.ebi.ac.uk/files/AF-",unip,"-F1-model_v2.pdb"))
  }
  list.files(clade_dir)
    ### from a Soy Uniprot collection write a fasta file with new sequences
    writeXStringSet(Maize_Uniprot[Uniprot_list$Uniprot],   paste0(clade_dir ,"uniprot.fasta"))
    ### use system command to run mafft -add
    system(paste0(" /Users/prigozhin/miniconda/envs/snakemake/bin/mafft --add ",clade_dir,"uniprot.fasta ",Uniprot_list$File[[1]]," > ",clade_dir,clade,".uniprot.afa"))
  maa <- readAAMultipleAlignment(paste0(clade_dir,clade,".uniprot.afa")) ## change to read the modified alignment
  for (unip in Uniprot_list$Uniprot){
    unip
    maa@unmasked@ranges@NAMES
    ent_ng <- get_ent(maa,unip)
    print_ent(ent_ng,paste0(clade_dir,unip,".ChimeraEntropy.txt"))
    pdf_ent(ent_ng,paste0(clade_dir,unip,".Entropy.pdf"),unip)
  }
}

NLR_Common %>%filter(HV ==1, !is.na(Uniprot)) %>% group_by(Clade_0)%>%count 

### Preparing hv graphs for publication ---

#4787 - RppM - "#d73027"
#6329 - RppC - "#fc8d59"
#6648 - Rp1 - "#4575b4"
#3480 - NA - "#6a329f"
prefix <- "AF_NLR/" 
pretty_print <- (NLR_models[c(1,2,4,7),])%>%mutate(Color =c("#6a329f","#d73027","#fc8d59","#4575b4") )
clade <- "Int3480_75_130_L_68"
for(clade in pretty_print$Clade%>%unique()){
    print(clade)
    clade_dir <- paste0(prefix,clade,"/")  
    list.files(clade_dir)
    maa <- readAAMultipleAlignment(paste0(clade_dir,clade,".uniprot.afa")) ## change to read the modified alignment
    unip <- pretty_print %>% filter(Clade == clade)%>%select(Uniprot)%>%unlist()
    color <- pretty_print %>% filter(Clade == clade)%>%select(Color)%>%unlist()
    maa@unmasked@ranges@NAMES
    ent_ng <- get_ent(maa,unip)
    #print_ent(ent_ng,paste0(clade_dir,unip,".ChimeraEntropy.txt"))
    #pdf_ent(ent_ng,paste0(clade_dir,unip,".Entropy.pdf"),unip)
    Ent <- as_tibble(cbind(1:nrow(ent_ng),ent_ng))
    colnames(Ent)<-c("Position","Entropy")
    ggplot(Ent, aes(x = Position))+
      geom_line(aes(y = Entropy), color = color)+
      ylim(0,3)+
      xlim(0,1400)+
      ggtitle(paste(unip)) +
      xlab("Position") + 
      ylab("Shannon Entropy")+
      theme_classic()
    ggsave(paste0("Figures/Figure_3/",unip,".PP.Entropy.pdf"))
}
pretty_print$Color
