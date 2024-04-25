require(tidyverse)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
require(Biostrings)
B73 <- readDNAStringSet("~/Downloads/Zm-B73-REFERENCE-NAM-5.0.fa")
B73@ranges@width[1:10]
B73@ranges@NAMES[1:10]
data <- data.frame(Chr = B73@ranges@NAMES[1:10],
                   size = B73@ranges@width[1:10])
#rm(B73)
Common <- read_tsv("~/Dropbox/NLRomes/Maize_NLRome/Maize_NLRome_GeneTable.txt")

GFF <- read_tsv("~/Dropbox/NLRomes/Maize_NLRome/Zm_NAM_gff/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3",comment = "#",col_names = c("Chr","source","type","start","end","score","strand","phase","info"))
Genes <- GFF %>% filter(type == "gene") %>% mutate(Gene = info %>%str_remove_all("ID=")%>%str_remove(";.*")%>%toupper())
B73_Common <- Common %>% mutate(Gene = Gene %>% str_remove("_.*"))%>% left_join(Genes, by = "Gene")%>%filter(!is.na(Chr))
  
B73_Common %>% select(Chr)%>%distinct%>%print(n=100)
B73_Common <- B73_Common %>% filter(Chr != "scaf_487")
max <- B73_Common %>% filter(HV ==1,Chr =="chr10")%>%select(end)%>%max
min <- B73_Common %>% filter(HV ==1,Chr =="chr10")%>%select(start)%>%min

adj <- 0.21
chrom_plot <- ggplot() +
  geom_segment(data = data,
               aes(x = fct_inorder(Chr), xend = Chr, y = 0, yend = size),
               lineend = "round", color = "lightgrey", size = 5) +
  geom_segment(data = B73_Common,
               aes(x = as.integer(Chr%>%str_remove("chr")) - adj, xend = as.integer(Chr%>%str_remove("chr")) + adj,
                   y = start, yend = end, color = as.factor(HV)),
               size = 1) +
  geom_segment(data = B73_Common %>% filter(HV == 1),
               aes(x = as.integer(Chr%>%str_remove("Chr")) - adj, xend = as.integer(Chr%>%str_remove("Chr")) + adj,
                   y = start, yend = end, color = as.factor(HV)),
               size = 1) +
  theme_classic() + theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none")  +
  scale_y_reverse()
ggsave(chrom_plot,filename = "Chrom_plot.pdf",height = 6.44,width = 4,dpi = 300,units = "in")
getwd()
# 
install.packages("devtools")
require(devtools)
devtools::install_github("thackl/thacklr")
devtools::install_github("thackl/gggenomes")
install.packages("reticulate")
require("thacklr")
require("gggenomes")
require(tidyverse)
library(reticulate)


##################################################
#### Trying with maize data -----------
B73_test_gff <- filter(GFF, Chr =="chr10" )
B73_gff <- read_gff3("~/Downloads/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3")
B73_test_gff <- B73_gff[2:2944,]

B73_test_gff %>% select(type)%>%distinct()

B73_NLR <- B73_gff%>%
  filter(type %in% c("gene","CDS"))%>%
  mutate(Gene = feat_id%>%toupper%>%str_remove("_.*")) %>% filter(Gene %in% B73_Common$Gene)
B73_norm <- B73_gff%>%
  filter(type %in% c("gene","CDS"))%>%
  mutate(Gene = feat_id%>%toupper%>%str_remove("_.*")) %>% filter(Gene %ni% B73_Common$Gene)
B73_test <- B73_gff%>%
  mutate(Gene = feat_id%>%toupper%>%str_remove("_.*"))%>%left_join(B73_Common%>%select(Gene,HV,starts_with("Clade")), by = "Gene")

to_plot <- B73_NLR %>% filter(seq_id == "chr10")
to_plot <- B73_test_gff
to_plot <- B73_norm
to_plot <- B73_test %>% filter(seq_id == "chr10")

the_plot <- gggenomes(to_plot%>%
                        filter(type %in% c("gene","CDS"))%>%
                        filter(start<max|end<max)%>%
                        filter(start>min|end>min)
                        #
                      )+geom_seq()+geom_gene(aes(color=as.factor(Clade_0),fill = as.factor(Clade_0))) #+ scale_x_bp(limits=c(15e5, 30e5))#+theme(legend.position="none")
the_plot + geom_bin_label()

TE_gff <- read_gff3(file = "~/Downloads/Zm-B73-REFERENCE-NAM-5.0.TE.gff3")
TE_reg <- TE_gff %>% 
  filter(seq_id =='chr1')%>%
  filter(start<3500000|end<3500000)%>%
  filter(start>1500000|end>1500000)

TE_gff %>% 
  filter(seq_id =='chr1')%>%
  filter(start<3500000|end<3500000)%>%
  filter(start>1500000|end>1500000)%>%
  select(start,end,strand,type)%>%
  print(n=300)
TE_reg <- TE_reg %>%mutate(strand = ifelse(strand =="?",NA,strand))
TE_reg <- TE_reg %>% mutate(kind=type)%>%mutate(type = "gene")
TE_reg <- TE_reg %>%mutate(type = "CDS")

gggenomes(TE_reg)+geom_seq()+geom_gene(aes(color=as.factor(kind),fill = as.factor(kind)))
TE_gff %>% group_by(type)%>%count

TE_gff %>% filter(type =="centromeric_repeat", seq_id =="chr10")%>%select(start)%>%ggplot(aes(x=start))+geom_density()+xlim(1,152435371)
TE_gff %>% filter( seq_id =="chr10")%>%ggplot(aes(x=start,color = type))+geom_density(scaled = 1)+xlim(1,152435371)

types <- TE_gff %>% select(type)%>%distinct()
TE_gff %>% filter( seq_id =="chr10",type %in% types$type[1:19])%>%ggplot(aes(x=start,color = type))+geom_density()+xlim(1,152435371)
data
TE_gff %>% filter(type =="knob", seq_id =="chr10")%>%select(start)%>%ggplot(aes(x=start))+geom_density()+xlim(1,152435371)
TE_gff %>% filter(type =="knob", seq_id =="chr10")%>%select(start)%>%ggplot(aes(x=start))+geom_density()+xlim(1,152435371)



#############################
### NAM GFF -----------------
setwd("~/Dropbox/NLRomes/Maize_NLRome/Zm_NAM_gff/")
gff_list <- list.files("~/Dropbox/NLRomes/Maize_NLRome/Zm_NAM_gff/",pattern = "^Zm.*.gff3$")
gff_list <- gff_list[2:27]
all_plot <- vector()
for (i in seq_along(gff_list)){
          gff_name <- gff_list[i]
          gff_name
          line_name <- (gff_name%>%str_split("-"))[[1]][2]
          gff <- read_delim(gff_name,delim = "\t",comment = "#",col_names = c("Chr","source","type","start","end","score","strand","phase","info"))
          genes <- gff %>% filter(type == "gene") %>% mutate(Gene = info %>%str_remove_all("ID=gene:")%>%str_remove_all("ID=")%>%str_remove(";.*")%>%toupper())
          genes$Gene
          (local_common <- Common %>% mutate(Gene = Gene %>% str_remove("_.*"))%>% left_join(genes, by = "Gene")%>%filter(!is.na(Chr)))
          
          #Common %>% mutate(Gene = Gene %>% str_remove("[AE].*"))%>%select(Gene) %>% distinct()%>%arrange(Gene)%>%print(n=30)
          
          if (nrow(local_common)<10){next}
          local_common %>% select(Chr) %>% distinct %>% arrange(Chr) %>% print(n=100)
          local_common <- local_common %>% filter(grepl(x = local_common$Chr,pattern = "chr"))
          chr <- (local_common %>% filter(HV ==1) %>% group_by(Chr) %>%count %>%arrange(-n))[[1,1]]
          #if (chr != "chr3"){next}
          max <- local_common %>% filter(HV ==1,Chr ==chr)%>%select(end)%>%max
          min <- local_common %>% filter(HV ==1,Chr ==chr)%>%select(start)%>%min
          gff_subset_name <- paste0("subset.",gff_name)
          gff_subset <- gff %>% filter(Chr == chr,
                                       start<max|end<max,
                                      start>min|end>min)
          write_delim(file = gff_subset_name,x = gff_subset,col_names = F,delim = "\t")
          gff_subset <- read_gff3(gff_subset_name)
          to_plot <- gff_subset%>%
            mutate(Gene = feat_id%>%toupper%>%str_remove("_.*"))%>%
            left_join(local_common%>%select(Gene,HV,starts_with("Clade")), by = "Gene")%>%
            mutate(seq_id = paste0(line_name," ",chr))
          all_plot <- rbind(all_plot,to_plot)
          # the_plot <- gggenomes(to_plot%>%
          #                         filter(type %in% c("gene","CDS"))%>%
          #                         filter(start<max|end<max)%>%
          #                         filter(start>min|end>min)
          #                       #
          # )+geom_seq()+geom_gene(aes(color=as.factor(Clade_0),fill = as.factor(Clade_0))) #+ scale_x_bp(limits=c(15e5, 30e5))#+theme(legend.position="none")
          # the_plot
          # plot_name <-  paste0(gff_name,".pdf")
          # ggsave(the_plot,filename = plot_name)
}

##### NEED TO CHECK HOW MANY OF THE RELEVANT CLADES ARE ON UNPLACED SCAFFOLDS
seq_ids <- all_plot %>%select(seq_id)%>%distinct()%>%print(n=101)

seq_order <- c(1,2,16,17,18,19,20,21,14,15,11,13,23,22,3,4,8,9,5,24,25,6,7,12,10,26)
length(seq_order)
order_tbl <- tibble(seq_id = seq_ids$seq_id,seq_order)
order_tbl <- order_tbl %>% arrange(seq_order)
order_tbl <- order_tbl %>% mutate(col = c("#ffd579",rep("#02b0f0",6),rep("#bfbfbf",3),"#d883ff","#ff9300","#ff9300",rep("#02b04f",13)))
all_plot <- all_plot %>%left_join(order_tbl)%>%arrange(seq_order)
col_vec <- c("#ffd579","#02b0f0",rep("#02b04f",8),"#d883ff","#ff9300",rep("#02b04f",2),rep("#02b0f0",2),rep("#bfbfbf",2),"#02b0f0",rep("#02b04f",2),rep("#02b0f0",2),"#ff9300","#bfbfbf","#02b04f")
col_vec %>% length
the_plot <- gggenomes(all_plot%>%filter(type %in% c("gene","CDS")))+
  geom_seq()+
  geom_gene(size = 3, aes(color=as.factor(Clade_0),fill = as.factor(Clade_0)))+
  geom_bin_label(color = col_vec)+
  scale_colour_manual(values =c("#fee090","#d73027","#fc8d59","#4575b4"))+
  scale_fill_manual(values = c("#fee090","#d73027","#fc8d59","#4575b4"))

#scale_fill_brewer("Genes", palette="Dark2", na.value="darkgrey")

plot_name <-  paste0("~/Dropbox/NLRomes/Maize_NLRome/Figures/Figure_2/test",".pdf")
#ggsave(the_plot,filename = plot_name)
#the_plot  %>% ggsave(filename = plot_name)
the_plot %>% pick(order_tbl$seq_id) %>% ggsave(filename = plot_name, width = 8,height =6)

#4787 - RppM - "#d73027"
#6329 - RppC - "#fc8d59"
#6648 - Rp1 - "#4575b4"
#3480 - NA - "#6a329f"


Common %>% filter(HV==1) %>%group_by(Clade_0)%>%count%>%print(n=30)

Common %>% filter(HV==1, Clade_0 =="Int3480_75") %>%group_by(Ecotype)%>%count%>%print(n=30)
Common %>% filter(HV==1, Clade_0 =="Int4787_129") %>%group_by(Ecotype)%>%count%>%print(n=30)
Common %>% filter(HV==1, Clade_0 =="Int6329_131") %>%group_by(Ecotype)%>%count%>%print(n=30)
Common %>% filter(HV==1, Clade_0 =="Int6648_150") %>%group_by(Ecotype)%>%count%>%print(n=30)

Common

cat(order_tbl$seq_id)
install.packages("cowplot")
require(cowplot)
plot_grid(chrom_plot,the_plot)

chr10_cluster <- all_plot

Common %>% filter(Clade_0 =="Int3480_75") %>% group_by(HV,Clade)%>%count
Common %>% filter(Clade_0 =="Int4787_129") %>% group_by(HV,Clade)%>%count
Common %>% filter(Clade_0 =="Int6329_131") %>% group_by(HV,Clade)%>%count
Common %>% filter(Clade_0 =="Int6648_150") %>% group_by(HV,Clade)%>%count
  

all_plot%>%filter(seq_id =="B73 chr10")%>%select(start)%>%filter(start >1)%>%unlist()%>%min
all_plot%>%filter(seq_id =="B73 chr10")%>%select(start)%>%filter(start >1)%>%unlist()%>%max
