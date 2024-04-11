## RNAseq processing 

The pipeline is outlined in `expression_calling_commands_per_tissue.bash`. Briefly, raw sequence files are downloaded using the custom `rna_download.py` script, then uncompressed using `unpigz.sh`. STAR reference genomes are created using `make_STAR_genome_dir.bash`. Reads are aligned to the respective reference genome using `STAR.bash`. TPM is calculated from the outputed counts file by `tpm_calc.py`. QC of the coverage of each NLR is determined by `samtools_coverage_runner.bash`. 

Parallel pipelines developed for RNA reads simulated by polyester (`simulate_experiment_runner.bash`) are also in this directory. 

Figures were generated in R studio. The `.Rmd` file contains all the code necessary to recreate the figures from intermediate data files available on Zenodo. The `.md` file is more user friendly to view in github and displays the plots alongside the code. The plots themselves are saved in the ./figure-gfm subdirectories.
