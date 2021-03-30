require("tidyverse")

Cut_Nodes_10<-read_delim(snakemake@input[["tsv"]],delim = " ",col_names = T)
BigTable <- read_delim(snakemake@input[["bigtable"]],delim = "\t",col_names = T)
OutputDirectory <- snakemake@output[["dir"]]

## Append a column to BigTable that contains the TRUE values for nodes in the cut list.
Cut_Nodes_10 <- mutate(Cut_Nodes_10, Split_Node_1 = T)
BigTable <- left_join(BigTable, Cut_Nodes_10 %>% select(label,node,Clade,Split_Node_1)%>%distinct(), by = c("label","node","Clade"))
BigTable %>% filter(Split_Node_1)
Cut_Nodes_10 %>% filter(Split_Node_1)%>%select(Clade)

## Generate Split_1 column with new clade assignment for every tip of every tree
Split_1 <-vector(length = nrow(BigTable))
for (i in 1:nrow(BigTable)){if (!is.na(BigTable[i,]$label)){                    #works
  (Tip <- BigTable[i,])
  (CurClade <- BigTable[i,]$Clade)
  (TreeData <- filter(BigTable,Clade == CurClade))
  (SplitNodes <- filter(TreeData,Split_Node_1))
  if (nrow(SplitNodes)>0){
    (TopClade <- SplitNodes %>% filter(R_tips == max(SplitNodes$R_tips)))
    if (is.na(Tip$Split_Node_1)){
      if (nrow(dplyr::intersect(ancestor(TreeData,Tip$node),SplitNodes))>0){
        SplitAncestors <- dplyr::intersect(ancestor(TreeData,Tip$node),SplitNodes)
        SmallestAncestor <- SplitAncestors %>% filter(R_tips == min(SplitAncestors$R_tips))
        BestNode <- paste0(CurClade,'_',SmallestAncestor$node,'_R')
      }else{BestNode <- paste0(CurClade,'_',TopClade$node,'_L')}
    }else{BestNode <- paste0(CurClade,'_',Tip$node,'_R')}
  }else{BestNode <- paste0(CurClade,"_NS_N")}
  Split_1[[i]] <- BestNode
}else{Split_1[[i]] <- NA}}

as_tibble(Split_1) %>% print(n=3000)
(BigTable_1 <- mutate(BigTable, Split_1 = as.factor(Split_1)))
BigTable_1 %>% select(label, Split_1) %>% print(n=1000)
## Count numbers of Maize tips in each new clade
CladeStat <- BigTable_1 %>% filter(!is.na(Split_1),!grepl('SORBI',label)) %>% group_by(Clade,Split_1) %>% summarise(n=n()) %>% print(n=300)
CladeStat %>% ungroup() %>% select(Clade) %>% distinct()
CladeStat %>% ungroup() %>% select(Split_1) %>% distinct()
BigTable  %>% select(label) %>% distinct()
write_delim(CladeStat, paste0(OutputDirectory,"/CladeStat.txt"), delim = "\t", append = F, col_names = T)

# ggplot(CladeStat,aes(x=n))+geom_density()
# ggplot(CladeStat,aes(x=n))+geom_density()+xlim(0,100)
# ggplot(CladeStat,aes(x=n))+geom_histogram()

# ## First Refinement. Print out new clade assignment for every tip ###This needs to change to get names identical to clade names below i.e. end in number of genes in the clade
 FirstAssignment <- BigTable_1 %>% filter(!is.na(Split_1)) %>% select(label,Split_1)
 write_delim(x = FirstAssignment, path = snakemake@output[["list"]] , delim = "\t",quote_escape = "double",append = F,col_names = T)

###Find duplicate tips
Duplicates <- BigTable_1 %>% group_by(label) %>%summarise(n=n()) %>% filter(n>1)
Duplicates <- Duplicates %>% filter(!is.na(label))
Duplicates
BigTable_1 %>% filter(label %in% Duplicates$label) %>% group_by(Split_1) %>% summarise(n=n())
BigTable_1 %>% filter(label %in% Duplicates$label) %>% arrange(label) %>% print(n=50) #%>% group_by(Clade) %>% summarise(n=n())

BigTable %>% filter(!is.na(Clade)) %>% select(Clade) %>% distinct()


write_delim(Duplicates, paste0(OutputDirectory,"/Duplicates.txt"), delim = "\t", append = F, col_names = T)

###Write Clade Lists - DO NOT RUN without repeating refinement
for (n in 1:(nrow(CladeStat))) {
  clade <- CladeStat[n,]$Split_1
  tips <- BigTable_1 %>% filter(Split_1 == clade)
  tipnames <- tips$label
  write_delim(x = as.data.frame(tipnames), path = paste0(OutputDirectory,"/",clade, "_",length(tips$label),".txt"), delim = "\t",quote_escape = "double",append = F,col_names = F)
}
save.image("Test_printrefined.RData")
