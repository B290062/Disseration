import argparse
import subprocess
import os 
import time
import zipfile
import glob 
import shutil
import pandas as pd 

print("################################################################### ")
print("# Welcome to the divergent transcription promoter analysis tool!  #")
print("# For usage information please type --help in the command terminal#")
print("################################################################### ")

parser = argparse.ArgumentParser()
parser.add_argument("--accession", action="store", required=True, help="The accession is the study number from Gene Expression Omnibus (e.g GSE87821)")
#store_true used for optional values, when the value is present it will be True so this can be used to write functions.
#https://docs.python.org/3/library/argparse.html#action
parser.add_argument("--trim", action="store_true", required=False, help = "The data can be optionally trimmed if the user requires")
parser.add_argument("--adapter1",help="The first adapter required for trimming")
parser.add_argument("--adapter2", help="The second adapter required for trimming") 
#Function to download fastq files from SRA 
def SRA_download(args): 
    #replaced with makedirs instead of os.mkdirs as it has the exist_ok function which prevents crashing
    #when the folder is already created.
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

def Quality_control(args):
    print('Beginning performing quality control...')
    print('---------------------------------')
    time.sleep(0.5)
    print('Creating new directory for FASTQC files...')
    print('----------------------------------')
    time.sleep(0.5)
    os.makedirs("FastQC", exist_ok = True)
    print("FastQC directory created")
    print('--------')


#running FastQC

    print('Running FastQC....')
    time.sleep(0.5)
    print('---------------')
    fastqc_run = subprocess.run('fastqc * -o FastQC', shell=True)
    if fastqc_run.returncode !=0:
        print('FastQC were not able to be produced')
        exit(1)
    else:
        print('Finished FastQC analysis')
        print('-------------------------')
        time.sleep(0.5)
        print('All the files are in the fastqc directory')
        print('---------------------------------------')


def Trimming(args):

#trimming opton Note- this is for simple trimming only, not linked adapter trimming
    
    if args.trim is False:
        print('Proceeding without trimming...')
        return
    
    if args.trim is True:
        os.makedirs("Trimmed_data", exist_ok=True)
        if not os.path.exists('SRA'):
            print('The sra directory does not exist. Exiting.')
            exit(1)
            
        for file in os.listdir('SRA'):
                if file.endswith('_1.fastq') or file.endswith('_1.fastq.gz'):
                    input_file = os.path.join('SRA', file)
                    input_file_2 = os.path.join('SRA', file.replace('_1.fastq', '_2.fastq').replace('_1.fastq.gz', '_2.fastq.gz'))
                    output_file = os.path.join('Trimmed_data', file.replace('.fastq', '_trimmed.fastq').replace('.fastq.gz', '_trimmed.fastq.gz'))
                    output_file_2 = os.path.join('Trimmed_data', file.replace('_1.fastq', '_2_trimmed.fastq').replace('_1.fastq.gz', '_2_trimmed.fastq.gz'))

                    cutadapt = subprocess.run('cutadapt -a ' + args.adapter1 + ' -A ' + args.adapter2 + ' -o ' + output_file + ' -p ' + output_file_2 + '' + input_file + '' + input_file_2, shell=True)
                    if cutadapt.returncode !=0:
                        print('Error occured')
                    else: 
                        print('The trimming is now done..')
                        print('-----------------')
                        
                    

def main():
    args = parser.parse_args()
    SRA_download(args)
    Quality_control(args)
    Trimming(args)

if __name__ == '__main__':
    main()