import argparse
import subprocess
import os 
import time
import zipfile
import glob 
import shutil
import pandas as pd 

print("################################################################### ")
print("# Welcome to the divergent transcription promoter analysis tool!   #")
print("# For usage information please type --help in the command terminal #")
print("################################################################### ")

parser = argparse.ArgumentParser()
parser.add_argument("--accession", action="store", required=True, help="the accession is the study number from Gene Expression Omnibus (e.g GSE87821)")
#Function to download fastq files from SRA 
def SRA_download(args): 
    os.makedirs("SRA", exist_ok=True)
    print('A Directory called SRA was created.')
    print('-----------------')
    os.chdir("SRA")
    print('Attempting to download fastq files from GEO')
    print('----------------------------------------------------')
    download = subprocess.run("esearch -db sra -query " +  args.accession + " | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | xargs fastq-dump  --skip-technical  --readids --read-filter pass --dumpbase --split-3", shell=True)
    if download.returncode !=0:
        print('Error occured during download of fastQ files')
        exit(1)
    else: 
        print('Download completed')
        print('------------------')
        time.sleep(0.5)
        
def main():
    args = parser.parse_args()
    SRA_download(args)

if __name__ == '__main__':
    main()