## EMseq processing 

The methylation pipeline is implemented by `methylation_calling_commands.sh`, with associated individual job scripts in this directory. `per_exon_methylation.py` is a custom python script to calculate per gene % methylation, and `methylation_combination.ipynb` combines the % methylation files into a single document. 

Figures were generated in R studio. The `.Rmd` file contains all the code necessary to recreate the figures from intermediate data files available on Zenodo. The `.md` file is more user friendly to view in github and displays the plots alongside the code. The plots themselves are saved in the `./figure-gfm` subdirectories.
