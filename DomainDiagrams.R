## ---------------------------
##
## Script name: DomainDiagrams.R
##
## Purpose of script:
##
## Author: Daniil Prigozhin
##
## Date Created: 2021-03-04
##
## Copyright (c) Daniil Prigozhin, 2020
## Email: daniilprigozhin@lbl.gov
##
## ---------------------------
##
## Notes:derives from Atha_NLRome_DomainDiagrams
##       aims to work inside a snakemake flow to integrate LRR-predictor output and non-redundant Pfam scan output
##       and produce iTol domain annotation file. 
## ---------------------------
## load packages

require(tidyverse)
require(Biostrings)
'%ni%' <- Negate('%in%') 
## set working directory

setwd("~/Dropbox/NLRomes/Maize_NLRome/")  


#### Get Pfam domain definitions --------------
Pfam_domains <- read_delim("HMM_search_pbNB-ARC/all_samples.pbNB-ARC.Pfam_scan_reduced.out",delim = "\t")
min_pfam_domains <- tibble(Gene = Pfam_domains$target_name,           
                           Dom =  Pfam_domains$query_name,
                           Start = Pfam_domains$env_from,
                           Stop = Pfam_domains$env_to,
                           Eval = Pfam_domains$fullseq_Evalue)
min_pfam_domains <- min_pfam_domains %>%filter(!grepl("LRR_",Dom))
min_pfam_domains

#### Get LRR predictor domain definitions --------------
if(!file.exists("Annotation/LRRpredictor.all.tsv")){
          dirs<-dir("Annotation/LRRpredictor/",pattern = "ZM.*_P001",full.names = T)
          files<-list.files(path = dirs, pattern = ".*.predshort.txt",full.names = T)
          all_lrr<-vector()
          for (jj in seq_along(files)){
            all_lrr <- rbind(all_lrr,read_delim(files[[jj]],delim = "\t"))
          }
          write_delim(path = "Annotation/LRRpredictor.all.tsv", x = all_lrr %>%select(-X17,-X24),delim = "\t",col_names = T)
          cat("=================================\nWrote file: Annotation/LRRpredictor.all.tsv\n")
}else {
          all_lrr <- read_delim("Annotation/LRRpredictor.all.tsv",delim = "\t",col_names = T)
          cat("=================================\nRead file: Annotation/LRRpredictor.all.tsv\n")
}

min_lrrpred <- tibble(Gene = all_lrr$`#Prot`,
                       Dom = "LRR",
                       Start = all_lrr$pos-5,
                       Stop = all_lrr$pos+15,
                       Eval = all_lrr$LRRpred)
min_lrrpred

#### Combine Pfam and LRR predictor domain definitions --------------
domains <- rbind(min_pfam_domains,min_lrrpred)

domains %>% group_by(Dom) %>% summarise(n=n()) %>% arrange (n) %>% print(n=300)
#domains %>% group_by(Dom) %>% summarise(n=n()) %>% arrange (n) %>% filter (!grepl("LRR", Dom),!grepl("NB-ARC", Dom)) %>% print(n=300)
#domains %>% mutate(Length = 1+Stop-Start) %>% filter(Length <100) %>% group_by(Dom) %>% summarise(n=n()) %>% arrange(-n)%>% print(n=300)
#domains %>% group_by(Dom) %>% summarise(n=n())%>%arrange(-n)%>% print(n=300)

### iTOL export----
# need to assign a color and a shape to every domain
# this can then be joined to the domains table and subset on the leaves in the tree
# finally sink+cat can be used to export iTOL-readable annotation file
main_doms <- tibble(Dom = c("pbNB-ARC","NB-ARC", "RPW8", "Rx_N", "TIR","TIR_2", "LRR"), 
                    Shape = c("HH","HH", "EL", "EL", "EL","EL", "RE"),
                    Color = c("#ffff99","#ffff99", "#6a3d9a", "#ff7f00", "#33a02c","#33a02c","#e31a1c"))

if (!file.exists("Annotation/All_Domains_with_colors.tsv")){
        
        extra_doms <- domains %>% filter(Dom %ni% main_doms$Dom) %>% group_by(Dom) %>% summarise(n=n()) %>% arrange (n) %>% print(n=300)
        
        ### 169 excluding LRR (LRRs from LRR predictor can be added later)
        ### RPPW, Rx_N, TIR, and NB-ARC are the most popular, The first three should share shape and have different colors
        ### NB-ARC and LRR should have different colors and shapes
        ### NB-ARC - yellow hor hexagon
        ### LRR - red rectangle
        ### Post-LRR - red diamond
        ### TIR - orange ellipse
        ### RPP8 - purple ellipse
        ### Rx_n - dark green ellipse
        ### Others  - randomly assign 12 divergent colors and shapes from the list: HV, TR, TL, PL, PR, PU, PD, OC
        
        shapes <- c("HV", "TR", "TL", "PL", "PR", "PU", "PD", "OC")
        #sample(shapes, 165, replace = T)
        #
        ### 12 Colors to use from Color Brewer----
        colors <- c("#a6cee3",
                    "#1f78b4",
                    "#b2df8a",
                    "#33a02c",
                    "#fb9a99",
                    "#e31a1c",
                    "#fdbf6f",
                    "#ff7f00",
                    "#cab2d6",
                    "#6a3d9a",
                    "#ffff99",
                    "#b15928")
 
        ### Building the color tables ---------
        color_doms <- tibble(Dom = extra_doms$Dom, 
                             Shape = sample(shapes, nrow(extra_doms), replace = T),
                             Color = sample(colors, nrow(extra_doms), replace = T)
        )
        
        all_doms <- rbind(color_doms,main_doms)
        write_delim(all_doms, "Annotation/All_Domains_with_colors.tsv", delim = "\t", col_names = T)

}else{all_doms<-read_delim("Annotation/All_Domains_with_colors.tsv",delim = "\t")}


#### Domain shapes for reference ----
#RE  rectangle
#HH  horizontal hexagon
#HV  vertical hexagon
#EL  ellipse
#DI  rhombus (diamond)
#TR  right pointing triangle
#TL  left pointing triangle
#PL  left pointing pentagram
#PR  right pointing pentagram
#PU  up pointing pentagram
#PD  down pointing pentagram
#OC  octagon
#GP  rectangle (gap; black filled rectangle with 1/3 normal height)


### Export to iTOL ------------
color_doms <- left_join(domains, all_doms, by = "Dom") 
color_doms <- color_doms %>% mutate(iTOL = paste(Shape, Start, Stop, Color, Dom, sep = "|"))
### Need to leftjoin a table of protein lengths to this and will be ready for export!!!
### For the original tree will also need to bind names/start-stop to match leaf IDs
fasta_files <- list.files(path = "Proteomes", pattern = "corr.*fasta", full.names = T)
lengths<-vector("list",length = length(fasta_files))
ii<-1
for (ii in seq_along(fasta_files)){
  print(paste0("Looking at file ",ii, " of ", length(fasta_files), " (",fasta_files[[ii]],")..."))
  a<-readAAStringSet(fasta_files[[ii]])
  lengths[[ii]] <- tibble(Gene = a@ranges@NAMES, Length = a@ranges@width)
}

prot_l<-lengths[[1]]
for (jj in 2:length(lengths)){prot_l<-rbind(prot_l,lengths[[jj]])}

color_doms <- left_join(color_doms,prot_l) %>% filter(!is.na(Length))
color_doms %>%arrange(Gene)

afa_file <- list.files(path = "HMM_search_pbNB-ARC",pattern = ".afa$",full.names = T)
a<-readAAStringSet(afa_file)
genreg <- tibble(GeneRegion=a@ranges@NAMES)
genreg <- genreg %>% mutate(Gene = str_remove(GeneRegion,"\\/.*$")%>% str_replace_all(" ","_"))


genreg<-left_join(genreg,color_doms, by = "Gene")

#%>% mutate(Gene = GeneRegion) %>% select(-GeneRegion) ->genreg


sink("Annotation/BigTree_domains.txt", append = F)
cat(
  "DATASET_DOMAINS
SEPARATOR COMMA
DATASET_LABEL,Domains
COLOR,#ff0000
DATA
")
for (gene in levels(as.factor(genreg$GeneRegion))){
  tbl <- genreg %>% filter(GeneRegion == gene)
  cat(paste(gene, tbl[[1,10]], sep = ","))
  
  for (ii in seq_along(tbl$Dom)){
    cat(",")
    cat(tbl[[ii,9]])
  }
  cat("\n")
}
sink()

sink("Annotation/SmallTrees_domains.txt", append = F)
cat(
  "DATASET_DOMAINS
SEPARATOR COMMA
DATASET_LABEL,Domains
COLOR,#ff0000
DATA
")
for (gene in levels(as.factor(genreg$Gene))){
  tbl <- genreg %>% filter(Gene == gene)
  cat(paste(gene, tbl[[1,10]], sep = ","))
  
  for (ii in seq_along(tbl$Dom)){
    cat(",")
    cat(tbl[[ii,9]])
  }
  cat("\n")
}
sink()

#### Annotate interesting groups of genes to help orienting the big trees ------------------

for(name in main_doms$Dom){
#name<-"RPW8"
sink(paste0("Annotation/",name,".genes.iTol.txt"),append = F)
cat(paste0(
  "DATASET_BINARY
SEPARATOR TAB

DATASET_LABEL	",name,"
COLOR	#ff0000

FIELD_SHAPES	1
FIELD_LABELS	",name,"
FIELD_COLORS	#ff0000


DATA

")
)
sink()
genreg %>% filter(Dom == name) %>% select(GeneRegion) %>% distinct()%>%mutate(Flag = 1) %>% write_delim(paste0("Annotation/",name,".genes.iTol.txt"),col_names = F,delim = "\t",append = TRUE)
}
#genreg %>% filter(Dom == "TIR_2") %>% select(GeneRegion) %>% distinct()%>%mutate(Flag = 1) %>% write_delim("Annotation/TIR_2.genes.iTol.txt",col_names = F,delim = "\t")

sink(paste0("Annotation/ID.genes.iTol.txt"),append = F)
cat(paste0(
  "DATASET_BINARY
SEPARATOR TAB

DATASET_LABEL	ID
COLOR	#ff0000

FIELD_SHAPES	1
FIELD_LABELS	ID
FIELD_COLORS	#ff0000


DATA

")
)
sink()
genreg %>% filter(Dom %ni% main_doms$Dom) %>% select(GeneRegion) %>% distinct()%>%mutate(Flag = 1) %>% write_delim("Annotation/ID.genes.iTol.txt",col_names = F,delim = "\t", append = TRUE)

#genreg %>% filter(GeneRegion %ni% (genreg %>% filter(Dom =="LRR") %>% select(GeneRegion) %>% unlist())) %>%select(GeneRegion) %>% distinct()%>%mutate(Flag = 1) %>% write_delim("Annotation/NoLRR.genes.iTol.txt",col_names = F,delim = "\t")

#### Make pdf domain diagrams for genes we want to plot in figures ----------
color_doms <- all_color_doms %>% filter(Gene %in% c(
  "A0A1D6MVM5",
  "A0A1D6IRB2",
  "A0A1D6IRQ2",
  "A0A1D6IRV2",
  "Example seq from Int4084",
  "Example seq from Int3480_75"
))
 

doms <- rbind(
  c("A0A1D6MVM5","CC",82,222,1365),
  c("A0A1D6MVM5","NB-ARC",257,593,1365),
  c("A0A1D6MVM5","LRR",595,1307,1365),
  c("A0A1D6IRB2","CC",1,142,1043),
  c("A0A1D6IRB2","NB-ARC",160,509,1043),
  c("A0A1D6IRB2","LRR",512,921,1043),
  c("A0A1D6IRQ2","CC",1,195,1192),
  c("A0A1D6IRQ2","NB-ARC",223,561,1192),
  c("A0A1D6IRQ2","LRR",565,1131,1192),
  c("A0A1D6IRV2","CC",1,39,1184),
  c("A0A1D6IRV2","NB-ARC",75,416,1184),
  c("A0A1D6IRV2","LRR",418,1178,1184)
)
doms <- tibble(Gene=doms[,1],Dom=doms[,2],Start=doms[,3]%>%col_integer,Stop=doms[,4],Length=doms[,5],col)
doms <- doms %>% mutate(Start=as.integer(Start)) %>%
  mutate(Stop=as.integer(Stop)) %>%
  mutate(Length=as.integer(Length)) 
  
main_doms[[4,1]]<-"CC"
color_doms <- doms %>% left_join(main_doms%>%select(Dom,Color))
color_doms %>%select(Gene)%>% distinct()
floor <- 0
spacing <- 40
graph_tbl <- vector()
for (gene in levels(as.factor(color_doms$Gene)) )     {
  #gene <- "ZM00001EB015450_P001"
  gene_tbl <- color_doms %>% filter(Gene == gene) %>%select(Gene,Dom,Start,Stop,Color,Length)
  all <- tibble(Gene = gene,Dom = "",Start = 1, Stop = gene_tbl[[1,6]],Color = "#000000",Length = gene_tbl[[1,6]])
  gene_tbl <- gene_tbl %>% mutate(ymin = floor+0,ymax=floor+10)
  all <- all  %>% mutate(ymin = floor+4,ymax=floor+6)
  tbl <- rbind(all,gene_tbl)%>%mutate(y_cent = floor+5)
  graph_tbl <- rbind(graph_tbl,tbl)
  floor<- floor+spacing
}

plot <- ggplot() + 
  geom_rect(aes(xmin = -400, xmax = tbl[[1,6]], ymin = 0, ymax = tbl[[1,6]]/8), 
            fill = "white")+
  geom_rect(graph_tbl, mapping = 
              aes(xmin = Start, xmax = Stop, 
                  ymin = ymin, ymax = ymax),
            fill = graph_tbl$Color,
            #alpha = 1,
            color = "black") + 
  geom_text(aes(label = graph_tbl$Dom,x = rowMeans(graph_tbl[,3:4]),y=graph_tbl$y_cent))+
  geom_text(aes(label = graph_tbl$Gene,x = -300, y = graph_tbl$y_cent))+
  geom_text(aes(label = graph_tbl$Length,x = graph_tbl$Length + 100 , y = graph_tbl$y_cent))+
  theme_void()
plot
ggsave(x=plot,"Figures/Drafts/test_plot.pdf")
getwd()



