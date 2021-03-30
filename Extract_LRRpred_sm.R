require(tidyverse)

dirs<-dir(snakemake@input[["lrrpred"]],pattern = "",full.names = T)
files<-list.files(path = dirs, pattern = ".*.predshort.txt",full.names = T)
all_lrr<-vector()
for (jj in seq_along(files)){
  all_lrr <- rbind(all_lrr,read_delim(files[[jj]],delim = "\t"))
}
write_delim(path = snakemake@output[["lrrpred"]], x = all_lrr %>%select(-X17,-X24),delim = "\t",col_names = T)
cat(paste0("=================================\nWrote file: ",snakemake@output[["lrrpred"]],"\n"))