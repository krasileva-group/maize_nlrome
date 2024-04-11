#download rna reads from ENA by tissue type 
import os 
import pandas as pd 
import time 
import sys
import glob
pd.options.mode.chained_assignment = None  # default='warn'

#read in metadata file 
rna_info=pd.read_csv('/global/home/users/chandlersutherland/e16/rna_metadata.csv', index_col=0)
sample=sys.argv[1]

#define function that takes the name of a sample, searches to see what has been downloaded, and returns a df of what's left 
def count_left_s(sample):
    paths=glob.glob(os.path.join('/global/scratch/users/chandlersutherland/e16/'+sample+'/rna_*/STAR/*_ReadsPerGene.out.tab'))
    subset=rna_info[rna_info['Sample']==sample]
    df=pd.DataFrame(paths)
    done=df.iloc[:,0].str.split('/', expand=True)
    done=done.iloc[:,9].str.strip('_ReadsPerGene.out.tab')

    fastq=subset['fastq_ftp'].str.split(';', expand=True)#.str.strip('.gz')
    #fastq.iloc[:,0]
    subset['fastq_r1']=fastq.iloc[:,0].str.split('/', expand=True).iloc[:,6].str.strip('.gz')
    subset['fastq_r2']=fastq.iloc[:,1].str.split('/', expand=True).iloc[:,6].str.strip('.gz')
    subset

    subset_long=pd.melt(subset, id_vars=['study_accession', 'sample_accession', 'experiment_accession', 'run_accession', 'tax_id', 'scientific_name', 'fastq_ftp', 'submitted_ftp', 'sra_ftp', 'Sample', 'tissue'], value_vars=['fastq_r1', 'fastq_r2'], var_name='fastq_path')
    subset_long
    left=subset_long[-subset_long['run_accession'].isin(done)]
    left=left.reset_index()[['study_accession', 'sample_accession', 'experiment_accession', 'run_accession', 'tax_id', 'scientific_name', 'fastq_ftp','submitted_ftp', 'sra_ftp', 'Sample', 'tissue', 'fastq_path', 'value']]
    to_do=len(left)
    print(sample+' '+str(to_do))
    return(left)

#function that takes in that df, then downloads using ftp 
def downloader(left):
    sample=left.iloc[0,9] # have to assign to get glob to work? 
    for run_accession in left['run_accession'].unique():
        corrupt= glob.glob("/global/scratch/users/chandlersutherland/e16/"+sample+"/rna_*/"+ run_accession+'*')
        print('found' + str(len(corrupt)) + 'corrupt files, removing...')
        for i in corrupt: 
                print("removing file "+ i)
                os.remove(i) 
                
    for i in range(0, len(left)):
        start_time = time.time()
        sample=left.iloc[i, 9]
        tissue=left.iloc[i, 10]
        os.system("mkdir -p /global/scratch/users/chandlersutherland/e16/"+sample+"/rna_"+tissue)
        os.chdir("/global/scratch/users/chandlersutherland/e16/"+sample+"/rna_"+tissue)
        
        list=left.iloc[i,6].split(';')[0].split('/')[0:6]
        base='/'.join(list)
        ex=left.iloc[i,12]
        
        os.system("wget ftp://"+base+'/'+ex+'.gz')
        file=str(left.iloc[i,12])+'.gz'
        if os.path.isfile(file) == False:
            print("download failed for", sample, tissue, file)
        else:
            end_time=time.time()
            print("finished rna download for ", sample, tissue,". Total time taken: ", end_time - start_time)
            
def STAR_run(left):
    sample=left.iloc[0,9]
    for tissue in left['tissue'].unique():
        STAR_command="sbatch --job-name="+sample+".STAR --export=base=/global/scratch/users/chandlersutherland/e16,sample="+sample+",tissue="+tissue+" -A co_minium --qos minium_htc4_normal --partition savio4_htc /global/home/users/chandlersutherland/e16/expression/STAR.bash"" -A co_minium --partition savio4_htc expression/STAR.bash"
        os.system(STAR_command) 
        print("running star for "+sample+" "+tissue)
        

left=count_left_s(sample)
downloader(left)

#unpigz_command="sbatch --job-name="+sample+".unpigz --export=base=/global/scratch/users/chandlersutherland/e16,sample="+sample+" -A co_minium -p savio4_htc --qos minium_htc4_normal /global/home/users/chandlersutherland/e16/unpigz.sh"
#os.system(unpigz_command)
#STAR_run(left)

print('\n\n\n')
#print('the following files for sample', sample, 'failed:')
#count_left_s(sample) 
