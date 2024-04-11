#download rna reads from ENA 
#using the middle of the 10th leaf 

import os 
import pandas as pd 
import time 

em_info=pd.read_csv('/global/home/users/chandlersutherland/e16/methylation_download_info.txt', sep='\t', header=0)
em_info[['sample', 'rep', 'corrected', 'S', 'R', '001']]=em_info['submitted_ftp'].str.split(pat=';',expand=True)[1].str.split('/',expand=True)[5].str.replace('-','_').str.split(pat='.', expand=True)[0].str.split('_',
expand=True)

#filter out the root and tassel samples 
leaf=em_info[(em_info['rep']=='rep1') | (em_info['rep']=='rep2') | (em_info['rep']=='rep3')]

#pretty messy metadata, so download everything (corrected, 3 reps, resequencing)

for i in range(38, len(leaf)):
    start_time = time.time()
    #set sample with for loop
    accession_name=leaf.iloc[i,9]
    
    #make dir 
    os.system("mkdir -p /global/scratch/users/chandlersutherland/e16/"+accession_name+"/em")
    os.chdir("/global/scratch/users/chandlersutherland/e16/"+accession_name+"/em")
    
    #want the fastq version of the files 
    r1=leaf.iloc[i,6].split(';')[0]
    r2=leaf.iloc[i,6].split(';')[1]
    
    os.system("wget ftp://"+r1)
    os.system("wget ftp://"+r2)
    end_time=time.time()
    print("finished rna download for accession", accession_name, ". Sample ", leaf.iloc[i,3], "Total time taken: ", end_time - start_time)
