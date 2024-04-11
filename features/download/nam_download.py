#download .fa genome file, the TE annotation file, and the gene annotation file for each NAM line 
#get assembly name from nam_genome_info.txt 
import os 
import pandas as pd 
import time 

nam_genome_info=pd.read_csv('/global/home/users/chandlersutherland/e16/nam_genome_info.txt', sep='\t\t', header=0)

for i in range(0, len(nam_genome_info)):
    start_time = time.time()
    #set sample with for loop
    sample=nam_genome_info.iloc[i,0]
    accession_name=sample.split('-')[1]

    #os.system("mkdir -p /global/scratch/users/chandlersutherland/e16/"+accession_name+"/genome")
    os.chdir("/global/scratch/users/chandlersutherland/e16/"+accession_name+"/genome")

    #get readme
    #os.system("wget https://download.maizegdb.org/"+sample+"/"+sample+".README")

    #get annotations 
    #os.system("wget https://download.maizegdb.org/"+sample+"/"+sample+"TE.gff3.gz")

    #get gene annotation
    #get cross reference name first, different for each genome
    annotation_name=nam_genome_info.iloc[i,3].split('->')[1]
    #os.system("wget https://download.maizegdb.org/"+sample+"/"+sample+"_"+annotation_name+".gff3.gz")

    #finally download genome 
    #os.system("wget https://download.maizegdb.org/"+sample+"/"+sample+".fa.gz")
    
    #download the cds file for polyester 
    #os.system("wget   https://download.maizegdb.org/"+sample+"/"+sample+"_"+annotation_name+".canonical.cds.fa.gz")
    #unpigz
    #os.system("unpigz *")
    
    #get gene info txt file 
    os.system("wget https://download.maizegdb.org/" +sample+"/"+annotation_name+".fulldata.txt")
    end_time=time.time()
    print("finished genome download for assembly ", accession_name, ". Total time taken: ", end_time - start_time)
