require(tidyverse)
setwd("~/Dropbox/NLRomes/Maize_NLRome/")
list.files("Annotation")
chrom_table <- read_delim("Annotation/chrom_genes.list",delim = "\t",col_names = c("Chrom","Gene"))
target_table <- read_delim("Annotation/taget.list",delim = "\t",col_names = c("GeneReg"))
target_table <- target_table %>% mutate(ShortCaps = str_remove(GeneReg,"/.*$"))
chrom_table <- chrom_table %>% mutate(ShortCaps = str_to_upper(Gene))
left_join(target_table,chrom_table)
colors <- c(
"#a50026",
"#d73027",
"#f46d43",
"#fdae61",
"#fee090",
"#e0f3f8",
"#abd9e9",
"#74add1",
"#4575b4",
"#313695")
chroms <- (paste0("chr",1:10))
col_chro <- tibble(Colors = colors,Chrom = chroms)
export <- left_join(target_table,chrom_table) %>% left_join(col_chro) %>%select(GeneReg,Colors,Chrom) %>%filter(!is.na(Colors))

sink("Annotation/ChromStrip.iTOL.txt", append = F)
cat("DATASET_COLORSTRIP
SEPARATOR SPACE
DATASET_LABEL Chromosomes
COLOR #ff0000
COLOR_BRANCHES 1
DATA
")
for (ii in seq_along(export$GeneReg)){
  cat(paste(export[ii,],sep = " "))
  cat("\n")
}



sink()
