require(tidyverse)
current_dir <- "/Volumes/DMP_Data_Drive/Soy_LRR_Pred/LRR_Exchange/"
write_file <- "~/Dropbox/NLRomes/Soy_NLRome/LRR_pred.tsv"
current_dir <- snakemake@input[["lrrpred"]]
write_dir <- snakemake@output[["lrrpred"]]


dirs<-dir(current_dir,pattern = "",full.names = T)
files<-list.files(path = dirs, pattern = ".*.predshort.txt",full.names = T)
all_lrr<-vector()
for (jj in seq_along(files)){
  all_lrr <- rbind(all_lrr,read_delim(files[[jj]],delim = "\t"))
}
write_delim(path = write_file, x = all_lrr %>%select(-X17,-X24),delim = "\t",col_names = T)
cat(paste0("=================================\nWrote file: ",write_file,"\n"))
