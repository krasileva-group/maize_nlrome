import os 
import pandas as pd
import sys
import glob
import numpy as np 

#untested 
#takes in sample name, and coagulates all the counts files into a single csv with TPM and log2TPM 
sample=sys.argv[1]

#load bed file 
bed_path=glob.glob('/global/scratch/users/chandlersutherland/e16/'+sample+'/genome/*all_gene.bed')
bed=pd.read_csv(bed_path[0], sep='\t', lineterminator='\n')
bed['gene_length'] = 1 + bed['chromEnd'] - bed['chromStart']

#get all counts files 
count_files=glob.glob('/global/scratch/users/chandlersutherland/e16/'+sample+'/rna_*/STAR/*ReadsPerGene.out.tab')

#populate a dataframe 
counts=pd.DataFrame()
for i in range(0, len(count_files)):
    accession=count_files[i].split('/')[6]
    rep=count_files[i].split('/')[-1].split('_')[0]
    tissue=count_files[i].split('/')[7].split('_')[1]
    
    #read in file 
    count=pd.read_csv(count_files[i], sep='\t', names=['name', 'unstranded', 'stranded_1', 'stranded_2'])
    count=count.drop(labels=['unstranded', 'stranded_2'], axis=1)
    
    #filter out summary stats and keep the column with 
    count=count[~count['name'].str.contains('N_')]
    
    #add accession and sample columns 
    count['accession']=accession
    count['rep']=rep
    count['tissue']=tissue
    
    counts=pd.concat([counts,count])
    
counts

#Merge counts output and gene length by gene name 
merged=pd.merge(bed, counts, on=['name'])
merged
#Calculate RPK by dividing readcount by length of gene in kb
merged['RPK']=merged['stranded_1']*1000/merged['gene_length']

#Calculate per million scaling factor for the sample 
#Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
grouped=merged.groupby(['rep'])

tpm=pd.DataFrame()
for name, group in grouped: 
    scaling=group['RPK'].sum()/1000000
    group['TPM']=group['RPK']/scaling 
    tpm=pd.concat([tpm, group])

tpm['log2(TPM)']=np.log2(tpm['TPM']+1)
tpm=tpm[['accession', 'rep', 'tissue', 'name', 'chrom', 'chromStart', 'chromEnd', 'strand', 'gene_length', 'stranded_1', 'TPM', 'log2(TPM)']]
tpm.to_csv('/global/scratch/users/chandlersutherland/e16/cs_reports/'+sample+'_all_tissue.tsv', sep='\t')
