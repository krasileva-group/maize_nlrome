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

#setwd("~/Dropbox/NLRomes/Maize_NLRome/")  
#setwd("~/Dropbox/NLRomes/Atha_SnkMk_NLRome")  


#### Get Pfam domain definitions --------------
Pfam_domains <- read_delim(snakemake@input[["pfam"]],
                           delim = "\t",
                           col_types = cols(
                             target_name = col_character(),
                             tlen = col_double(),
                             query_name = col_character(),
                             qlen = col_double(),
                             fullseq_Evalue = col_double(),
                             dom_N = col_double(),
                             dom_of = col_double(),
                             dom_cEvalue = col_double(),
                             hmm_from = col_double(),
                             hmm_to = col_double(),
                             ali_from = col_double(),
                             ali_to = col_double(),
                             env_from = col_double(),
                             env_to = col_double(),
                             hmm_frac = col_double()
                           ))
min_pfam_domains <- tibble(Gene = Pfam_domains$target_name,           
                           Dom =  Pfam_domains$query_name,
                           Start = Pfam_domains$env_from,
                           Stop = Pfam_domains$env_to,
                           Eval = Pfam_domains$fullseq_Evalue)
min_pfam_domains <- min_pfam_domains %>%filter(!grepl("LRR_",Dom))
cat ("========================================\nRead reduced PFAM table:\n")
min_pfam_domains

#### Get LRR predictor domain definitions, save results for future use --------------
# all_lrr <- read_delim("~/Dropbox/NLRomes/Atha_SnkMk_NLRome/Annotation/all_samples.LRRpred.tsv",
#                       col_names = c("Prot","pos","clf1","clf2","clf3","clf4","clf5","clf6","clf7","clf8","LRRpred","-5","-4","-3","-2","-1", "SP1","1","2","3","4","5","6","SP2","+6","+7","+8","+9","+10"),
#                       delim = "\t",comment = "#")
all_lrr <- read_delim(snakemake@input$lrrpred,
                      col_names = c("Prot","pos","clf1","clf2","clf3","clf4","clf5","clf6","clf7","clf8","LRRpred","-5","-4","-3","-2","-1", "SP1","1","2","3","4","5","6","SP2","+6","+7","+8","+9","+10"),
                      delim = "\t",comment = "#")
#all_lrr <- read_delim(snakemake@input$lrrpred,col_names = TRUE, delim = "\t",comment = "#")


## Rice version:
# min_lrrpred <- tibble(Gene = all_lrr$Prot %>% str_replace("NC_","NC@") %>% str_replace("_",".") %>% str_replace("NC@","NC_"),
#                       Dom = "LRR",
#                       Start = all_lrr$pos-5,
#                       Stop = all_lrr$pos+15,
#                       Eval = all_lrr$LRRpred)

min_lrrpred <- tibble(Gene = all_lrr$Prot,
                       Dom = "LRR",
                       Start = all_lrr$pos-5,
                       Stop = all_lrr$pos+15,
                       Eval = all_lrr$LRRpred)
cat ("========================================\nRead LRR Predictor results:\n")
min_lrrpred

#### Combine Pfam and LRR predictor domain definitions --------------
domains <- rbind(min_pfam_domains,min_lrrpred)

cat ("========================================\nDomains Observed:\n")

domains %>% group_by(Dom) %>% summarise(n=n()) %>% arrange (n) %>% print(n=300)

write_delim(domains,path = snakemake@output[["domains"]],col_names = TRUE)

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

}else{all_doms<-read_delim("Annotation/All_Domains_with_colors.tsv",
                           delim = "\t", 
                           col_types = cols(
                                            Dom = col_character(),
                                            Shape = col_character(),
                                            Color = col_character()
                                            )
                           )
}


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
all_color_doms <- left_join(domains, all_doms, by = "Dom")
all_color_doms <- all_color_doms %>% mutate(iTOL = paste(Shape, Start, Stop, Color, Dom, sep = "|"))
### Need to leftjoin a table of protein lengths to this and will be ready for export!!!
### For the original tree will also need to bind names/start-stop to match leaf IDs

cat ("========================================\nExtracting Protein Lengths:\n")

fasta_files <- snakemake@input$fasta
lengths<-vector("list",length = length(fasta_files))
#save.image(file = "~/Dropbox/NLRomes/Soy_NLRome/Test.RData")
#setwd("~/Dropbox/NLRomes/Soy_NLRome/")

for (ii in seq_along(fasta_files)){
  print(paste0("Getting protein lengths from file ",ii, " of ", length(fasta_files), " (",fasta_files[[ii]],")..."))
  a<-readAAStringSet(fasta_files[[ii]])
  lengths[[ii]] <- tibble(Gene = str_trim(a@ranges@NAMES), Length = a@ranges@width)
}

prot_l<-lengths[[1]]
for (jj in 2:length(lengths)){prot_l<-rbind(prot_l,lengths[[jj]])}

all_color_doms <- left_join(all_color_doms,prot_l,by = "Gene") %>% filter(!is.na(Length))

afa_file <- snakemake@input[["afa"]]
a<-readAAStringSet(afa_file)
genreg <- tibble(GeneRegion=a@ranges@NAMES)
genreg <- genreg %>% mutate(Gene = str_remove(GeneRegion,"\\/.*$")%>% str_replace_all(" ","_"))


genreg<-left_join(genreg,all_color_doms, by = "Gene")

#%>% mutate(Gene = GeneRegion) %>% select(-GeneRegion) ->genreg

cat ("========================================\nExporting iTOL Domain Annotations:\n")

sink(snakemake@output[["bigtree"]], append = F)
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
cat ("Exported iTOL Domains for Big Tree (ProteinName/nnn-mmm)\n")

sink(snakemake@output[["smalltree"]], append = F)
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
cat ("Exported iTOL Domains for small Trees (ProteinName)\n")

#### Annotate interesting groups of genes to help orienting the big trees ------------------
cat ("========================================\nExporting iTOL Binary Annotations:\n")

for(name in main_doms$Dom){
  cat(paste0("Exporting ", name, " binary annotation\n"))
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
cat(paste0("Exporting non-standard domain (ID) binary annotation\n"))
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
cat ("========================================\n\n")

#genreg %>% filter(GeneRegion %ni% (genreg %>% filter(Dom =="LRR") %>% select(GeneRegion) %>% unlist())) %>%select(GeneRegion) %>% distinct()%>%mutate(Flag = 1) %>% write_delim("Annotation/NoLRR.genes.iTol.txt",col_names = F,delim = "\t")



#save.image("Test.RData")
