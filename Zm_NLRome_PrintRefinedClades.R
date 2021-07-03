require("parallel")
require("tidyverse")
require("ggtree")
require("treeio")
require("tidytree")


#setwd("~/Dropbox/NLRomes/Maize_NLRome/")
#setwd("~/Dropbox/NLRomes/Rice_NLRome/")

Cut_Nodes_10<-read_delim(snakemake@input[["tsv"]],delim = " ",col_names = T)
load(snakemake@input[["BigTable"]])
OutputDirectory <- snakemake@output[["dir"]]
if(!dir.exists(OutputDirectory)){dir.create(OutputDirectory)}
save.image("Test_printrefined.RData")
#load("SM.RData")
#stop()

## Append a column to BigTable that contains the TRUE values for nodes in the cut list.
Cut_Nodes_10 <- mutate(Cut_Nodes_10, Split_Node_1 = T)
BigTable <- left_join(BigTable, Cut_Nodes_10 %>% select(label,node,Clade,Split_Node_1)%>%distinct(), by = c("label","node","Clade"))
cat("==============================\nClades to refine this round:")
BigTable %>% filter(Split_Node_1) %>% select(Clade) %>% distinct() %>% print(n=1000)

## Generate Split_1 column with new clade assignment for every tip of every tree
Split_1 <-vector(length = nrow(BigTable))
i<-1
for (i in 1:nrow(BigTable)){if (!is.na(BigTable[i,]$label)){                    #works
  (Tip <- BigTable[i,])
  (CurClade <- BigTable[i,]$Clade)
  (TreeData <- filter(BigTable,Clade == CurClade))
  (SplitNodes <- filter(TreeData,Split_Node_1))
  if (nrow(SplitNodes)>0){
    (TopClade <- SplitNodes %>% filter(R_tips == max(SplitNodes$R_tips)))
    TopClade <- TopClade[1,]
    if (is.na(Tip$Split_Node_1)){
      if (nrow(dplyr::intersect(ancestor(TreeData,Tip$node),SplitNodes))>0){
        SplitAncestors <- dplyr::intersect(ancestor(TreeData,Tip$node),SplitNodes)
        SmallestAncestor <- SplitAncestors %>% filter(R_tips == min(SplitAncestors$R_tips))
        BestNode <- paste0(CurClade,'_',SmallestAncestor$node,'_R')
      }else{BestNode <- paste0(CurClade,'_',TopClade$node,'_L')}
    }else{BestNode <- paste0(CurClade,'_',Tip$node,'_R')}
  }else{BestNode <- paste0(CurClade)}
  Split_1[[i]] <- BestNode
}else{Split_1[[i]] <- NA}}

BigTable_1 <- mutate(BigTable, Split_1 = as.factor(Split_1))
BigTable_1 %>% select(label, Split_1) %>% print(n=1000)
## Count numbers of Maize tips in each new clade
CladeStat <- BigTable_1 %>% filter(!is.na(Split_1),!grepl('SORBI',label)) %>% 
                            group_by(Clade,Split_1) %>% summarise(n=n()) %>%
                            mutate(Split = ifelse(Clade == Split_1,Clade,paste0(Split_1,"_",n))) %>% print(n=300)
CladeStat %>% ungroup() %>% select(Clade) %>% distinct()
CladeStat %>% ungroup() %>% select(Split) %>% distinct()

CladeStat %>% ungroup() %>%filter(n > 3,Split != Clade) %>% select(Split) %>% distinct() %>% write_delim("pbNB-ARC_Clade.r1.tree.list",col_names = F, append = T)

CladeStat %>% ungroup() %>% group_by(Clade,Split) %>% count() %>% filter(Clade == Split)%>%print(n=1000)
BigTable_1 <- left_join(BigTable_1,CladeStat %>% select(-n)) %>% select(-Split_1)
#BigTable  %>% select(label) %>% distinct()
#write_delim(CladeStat, paste0(OutputDirectory,"/CladeStat.txt"), delim = "\t", append = F, col_names = T)

# ggplot(CladeStat,aes(x=n))+geom_density()
# ggplot(CladeStat,aes(x=n))+geom_density()+xlim(0,100)
# ggplot(CladeStat,aes(x=n))+geom_histogram()

# ## First Refinement. Print out new clade assignment for every tip ###This needs to change to get names identical to clade names below i.e. end in number of genes in the clade
 FirstAssignment <- BigTable_1 %>% filter(!is.na(Split)) %>% select(label,Split)
 write_delim(x = FirstAssignment, file = snakemake@output[["list"]] , delim = "\t",quote_escape = "double",append = F,col_names = T)

###Find duplicate tips
Duplicates <- BigTable_1 %>% group_by(label) %>%summarise(n=n()) %>% filter(n>1)
Duplicates <- Duplicates %>% filter(!is.na(label))
Duplicates %>% print(n=10000)
BigTable_1 %>% filter(label %in% Duplicates$label) %>% group_by(Split) %>% summarise(n=n())
BigTable_1 %>% filter(label %in% Duplicates$label) %>% arrange(label) %>% print(n=50) 


write_delim(BigTable_1 %>% filter(label %in% Duplicates$label) %>% group_by(Split) %>% summarise(n=n()), paste0(OutputDirectory,"/OverlappingClades.txt"), delim = "\t", append = F, col_names = T)

###Write Clade Lists - DO NOT RUN without repeating refinement
n<-1
for (n in 1:(nrow(CladeStat))) {
  clade <- CladeStat[n,]$Split
  tips <- BigTable_1 %>% filter(Split == clade)
  tipnames <- tips$label
  if(tips[1,]$Clade != tips[1,]$Split){
    write_delim(x = as.data.frame(tipnames), file = paste0(OutputDirectory,"/",clade,".txt"), delim = "\t",quote_escape = "double",append = F,col_names = F)
  }
}  

save.image("Test_printrefined.RData")
