import os 
import sys
import pandas as pd 
import glob 
import numpy as np 

pd.options.mode.chained_assignment = None  # default='warn'

#this python file does the per gene processing after bismark methylation extractor, returning a per gene average methylation
#In the case of CG context, the symmetrical CGs are averaged as they are not statistically independent  
#The output is a tsv file with the average methylation % per gene, and the number of high coverage cytosines (or CG sites for CG context) per gene 

#input is an unzipped .cov file 
cov_file=str(sys.argv[1])
methylation = pd.read_csv(cov_file, skiprows=1, sep = '\t', names=['Chrom', 'start_pos', 'end_pos', 'meth_percentage', 'count_methylated', 'count_unmethylated']).sort_values('Chrom')
filename=cov_file.split('/')[-1].split('.')[0]
context=filename.split('_')[0]
sample=cov_file.split('/')[6]
print(filename+' loaded successfully')

#initialize output directory (same as input directory)
f=len(cov_file.split('/'))
b=cov_file.split('/')[0:f-1]
out_dir='/'.join(b)
print('Output directory is '+out_dir)

#input necessary exon coordiate name files variables (bed files, chromosomes)
#second passed argument should be the bed file for the entire genome
chroms=methylation['Chrom'].unique()
pos_file=str(sys.argv[2])
all_positions = pd.read_csv(pos_file, sep = '\t',header=0,  names=['index', 'Chrom', 'start', 'end', 'gene', 'exon', 'strand'], index_col=False).drop(['index'], axis=1)
all_positions=all_positions[['Chrom', 'start', 'end', 'strand', 'gene', 'exon']].sort_values('Chrom')

#step 1: filter out ChrC and ChrM, since we don't care about them 
print('Filtering out ChrC and ChrM...')
met=methylation.loc[(methylation['Chrom']!='ChrC') & (methylation['Chrom']!='ChrM')] 

#step 2: Average the symmetrical CG pairs since they are not independent 
def cg_average(df):
    cpg_cov=df.sort_values(['Chrom','start_pos'])
    
    avg_cpg = []
    
    for i in range(0,len(cpg_cov)-1):
    #find and average the symmetrical CpG pairs 
    #this relies on sorting by chromosome and then by position, searching for neighbors
        if cpg_cov.iloc[i, 1] == cpg_cov.iloc[i+1,1]-1:
            obj = {'Chrom':cpg_cov.iloc[i,0],
                   'start_pos':cpg_cov.iloc[i,1],
                   'end_pos':cpg_cov.iloc[i+1,1],
                   'meth_percentage':np.mean([cpg_cov.iloc[i,3],
                                              cpg_cov.iloc[i+1,3]])
            }
            avg_cpg.append(obj)
        else:
            continue
    avg_cpg_df = pd.DataFrame.from_dict(avg_cpg)
    #avg_cpg_df['count'] == 'NA'
    return avg_cpg_df

if context=='CpG':
    print('CpG context detected, averaging symmetrical positions')
    met=cg_average(met)
    print('averaging complete')
else: 
    print('CHH or CHG context, no averaging necessary')
        

#step 3: assign cytosines to genes. Define a few functions first 

#takes in a dataframe and a list of chromosomes, and returns a list of dataframes broken by chromosome.
#Run on positions to create pos 
def break_chroms(df, chroms):
  df_chroms = []
  for i in chroms: 
    chrom_met = df[df['Chrom'] == i]
    df_chroms.append(chrom_met)
  return df_chroms

#takes in each individual chunk by chromosome, and "names" it, assigning it to the correct exon  
def exon_namer(pos, met):
  for j in range(0, len(pos)):
    start = pos.iloc[j, 1]
    end = pos.iloc[j, 2]
    exon = pos.iloc[j, 5]
    met_named=met #as written now, this overwrites met. 
    met_named.loc[
      (met_named['start_pos'].between(start,end)), 'Exon'] = exon
  return met_named.dropna()

pos_chroms = break_chroms(all_positions, chroms)  
met_chroms= break_chroms(met, chroms)
print('Chromosomes broken')

def combine(pos_chroms, met_chroms, chrom_number):
  print('chrom number: '+str(chrom_number))

  #check if there are exons on this chromosome 
  if pos_chroms[chrom_number].empty== True:
    print('no exons on this chromosome')
    final=pd.DataFrame()
  
  #check if there are methylated cytosines on this chromosome
  elif met_chroms[chrom_number].empty== True:
    print('no methylation on this chromosome')
    final=pd.DataFrame()

  #if there are, run exon namer 
  else: 
    chrom_name=pos_chroms[chrom_number]['Chrom'].unique()[0]
    print('chrom name: '+chrom_name)
    print('starting '+pos_chroms[chrom_number].iloc[0,0])

    final=exon_namer(pos_chroms[chrom_number], met_chroms[chrom_number])
    
    #check if there are methylated chromosomes inside exons 
    if len(final) == 0: #case where there are no methylated cytosines in positions 
      final=pd.DataFrame()

    else: 
      #print('some cytosines found')
      final[['Gene', 'transcript_exon']]=final['Exon'].str.split('_', expand=True)
      final[['transcript', 'del', 'exon_num']]=final['transcript_exon'].str.split('.', expand=True)
      final=final[final['transcript']=='T001'] #filter to just T1 to not overrepresent any exons 

      print('aggregate '+chrom_name+' over exons')
      final=final.groupby(['Chrom', 'Gene']).agg({'meth_percentage':['mean', 'count']}).reset_index()
  return final 

all_gene_count = pd.concat([combine(pos_chroms, met_chroms, i) for i in range(0, len(chroms))])

print('Gene names assigned')

#step 4: save file to a tsv 
all_gene_count.to_csv(out_dir+'/'+filename+'_per_exon_met_'+sample+'.tsv', sep='\t', header=True, index=False)
print('Gene counts saved to '+out_dir+'/'+filename+'_per_exon_met_'+sample+'.tsv')

print('Finished with '+sample+' in the '+context+' context!!')