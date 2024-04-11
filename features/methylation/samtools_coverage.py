import os 
import sys
import pandas as pd 

#read in NLR bed file, and take bam file from bash input 
NLR_position=pd.read_csv(sys.argv[1], sep='\t', index_col=False)
NLR_position=NLR_position.drop(NLR_position.columns[[0]], axis=1)
print('NLR_position file loaded')

bismark_file=str(sys.argv[2])
print('starting bismark file '+bismark_file)
basename=bismark_file.split('/')[-1].replace('.bam', '')
print(basename) 

#create header 
os.system("echo 'rname  startpos        endpos  numreads        covbases        coverage        meandepth       meanbaseq       meanmapq' >> "+basename+"_coverage.tsv")

print('starting for loop')

#loop the coverage function over the chromosome coordinates given in the NLR file 
for i in range(1, len(NLR_position)):
    chromosome=NLR_position.iloc[i,0]
    start=NLR_position.iloc[i,1]
    end=NLR_position.iloc[i,2]
    call=str(chromosome) +':' + str(int(start)) + '-' + str(int(end))
    os.system('samtools coverage -r '+call+' '+bismark_file+'  >> '+basename+'_coverage.tsv')

os.system("sed '/^#/d' "+basename+"_coverage.tsv > "+basename+"_clean_coverage.tsv") 

#add gene names
cov=pd.read_csv(basename+"_clean_coverage.tsv", sep='\t', skiprows=1, names=['rname', 'startpos', 'endpos', 'numreads', 'covbases', 'coverage', 'meandepth', 'meanbaseq', 'meanmapq'])
merged=pd.merge(cov, NLR_position, how='left', left_on=['rname', 'startpos', 'endpos'], right_on=['chrom', 'chromStart', 'chromEnd'])
merged=merged.drop(merged.columns[[9, 10, 11, 13]], axis=1)
cols=['name', 'rname', 'startpos', 'endpos', 'numreads', 'covbases', 'coverage', 'meandepth', 'meanbaseq', 'meanmapq']
merged=merged[cols]

merged.to_csv(basename+"_clean_coverage.tsv", sep='\t') 

