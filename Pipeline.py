import argparse
import subprocess
import os 
import time
import zipfile
import glob 
import shutil
import pandas as pd 
import configparser

print("################################################################### ")
print("# Welcome to the Divergent Transcription Promoter Analysis Tool   #")
print("# For usage information please type --help in the command terminal#")
print("################################################################### ")

config = configparser.ConfigParser()
config.read('config.ini')

parser = argparse.ArgumentParser()
parser.add_argument("--sra", action="store", required=True, help="This is the SRA number from Gene Expression eg for GSE87821 it would be SRP091444)")
#store_true used for optional values, when the value is present it will be True so this can be used to write functions.
#https://docs.python.org/3/library/argparse.html#action
parser.add_argument("--trim", action="store_true", required=False, help = "The data can be optionally trimmed if the user requires")
parser.add_argument("--adapter1",required=False,help="The first adapter required for trimming")
parser.add_argument("--adapter2",required=False, help="The second adapter required for trimming") 
parser.add_argument("--multiqc", action="store_true", required=False, help = "This argument enables the use of MutliQC")
#Function to download fastq files from SRA 
parser.add_argument("--fasta", default=config["urls"]["fasta_url"], help = "This is the link/path to the FASTA file from Ensembl it can be" \
"changed in the command line or via pasting a new link/path into the config file")
parser.add_argument("--gtf", default =config["urls"]["gtf_url"], help = "Link to the GTF file from Ensembl, can be overwriten in " \
"the config or changed via the command line.")

def SRA_download(args): 
    #replaced with makedirs instead of os.mkdirs as it has the exist_ok function which prevents crashing
    #when the folder is already created.
    os.makedirs("SRA", exist_ok=True)
    print('A Directory called SRA was created.')
    print('-----------------')
    os.chdir("SRA")
    print('Attempting to download fastq files from GEO')
    print('----------------------------------------------------')
    #if statement adapted from https://www.geeksforgeeks.org/python/check-if-directory-contains-file-using-python/
    #this checks if the current directory "SRA" contains one or more files and skips the download process if true
    #essential as the SRA files for any given study are very large.
    if len(os.listdir(".")) > 0:
        print("SRA files were detected, inside the SRA directory, therefore the download will be skipped")
        return
    download = subprocess.run("esearch -db sra -query " +  args.sra + " | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | xargs fastq-dump  --skip-technical  --readids --read-filter pass --dumpbase --split-3", shell=True)
    if download.returncode !=0:
        print('Error occured during download of fastQ files')
        exit(1)
    else: 
        print('Download completed')
        print('------------------')
        time.sleep(0.5)

def Quality_control(args):
    #produces FastQC files which allow the user to examine to determine the quality of the data.
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
    #trimming option Note- this is for simple trimming only, not linked adapter trimming
    
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
                        
def Multiqc(args):
    #this method takes all of the FastQC files and combines them into a MultiQC for easier interpretation
    if args.multiqc is False:
        print('MultiQC was not performed')
        return
    if args.multiqc is True:
        print('Performing MultiQC....')
        time.sleep(0.5)
        print('--------------------------')
        multiqc_run = subprocess.run('multiqc .', shell = True)
        if multiqc_run.returncode !=0:
            print('Error occured while running MultiQC')
            exit(1)
        else:
            print('-----------------------------')
            print('Finished MultiQC')
            time.sleep(0.5)
            print('MultiQC directory was produced containing the results')
            time.sleep(0.5)
            print('---------------------------')                    

def STAR_files_fasta(args):
    #this function enables the download of the STAR FASTA genome files
    #checks if .gz files already exist in the directory to skip the download
    #this was adapted from the previous code from the unzip function.
    gz_files = glob.glob('*.gz')
    if len(gz_files) > 0:
        print("Fasta and GFT files detected. skipping the download")
        return
    print(f'Downloading assembly {args.fasta}')
    print('--------------------------')
    fatsa_download = subprocess.run('wget ' +  args.fasta, shell = True)
    if fatsa_download.returncode !=0:
        print('Error occured while attempting to download FASTA file...')
        exit(1)
    else:
        print('FASTA download successful')
        print('-----------------------------')
        time.sleep(0.5)
        
def STAR_files_GTF(args):
    # function provides the annotation coordinates
    print(f'Downloading GTF {args.gtf}')
    print('--------------------------')
    GTF_download = subprocess.run(' wget ' + args.gtf, shell = True)
    if GTF_download.returncode !=0:
        print('Error occured while attempting to download GTF file...')
        exit(1)
    else: 
        print('GTF download successful')
        print('-----------------------------')
        time.sleep(0.5)

def Unzip(args):
    # STAR files and compressed when downloaded and therefore need to be unzipped. 
    time.sleep(0.5)
    gz_files = glob.glob('*.gz')
    if not gz_files:
        print('Unzipping of STAR files failed, please try again.')
        exit(1)
    else:
        print('Begin to unzip files')
        print('--------------------')
        unzip = subprocess.run(['gunzip'] + gz_files)
        if unzip.returncode !=0:
            print('The zipped files were unable to be unzipped.')
            exit(1)
        else:
            print('Both files unziped')
            print('-------------------')
            time.sleep(0.5)
            files = os.listdir('.')
            print('Files in the current directory:')
            print('-------------------------------')
            for file in files:
                print(file)

def Indexing():
    print('Now starting to build index')
    time.sleep(0.5)
    while True:
        print('Please provide full pass to the files from your diretory')
        print('Example: /home/s2614505/Diss/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa')
        time.sleep(0.5)
        fa_path = input('Please provide full path for fa file: ')
        if os.path.isfile(fa_path):
            print('The file exists', fa_path)
            time.sleep(0.5)
            break
        else: 
            print('File not found, please try again')

    while True:  
        print('Example for GTF: /home/s2614505/Diss/Mus_musculus.GRCm39.112.gtf ')
        GTF_path = input('Please provide path for GTF file: ')
        if os.path.isfile(GTF_path):
            print('The file exists', GTF_path)
            time.sleep(0.5)
            break
        else: 
            print('File not found, please try again')
    
    while True:
        print('Please provide a full path to ypur home directory where sra and ref directories can be found')
        print('Example: /home/s2614505/Diss/')
        home_path = input('Path: ')
        if os.path.isdir(home_path):
            print('Path is valid proceeding...')
            print('Proceding...')
            time.sleep(0.5)
            os.chdir(home_path)
            break 
        else: 
            print('Not valid path, please try again')
    
    directory_ref = ('ref')
    try:
        os.mkdir(directory_ref)
    except OSError as error:
        print('the error:', error)
        time.sleep(0.5)
        print('Directory called ' + directory_ref + ' craeted' )
        print('Results of indexing will be stored there')
        time.sleep(0.5)
        print('--------')

        
    index = subprocess.run('STAR --runMode genomeGenerate --genomeDir ' + directory_ref + '/ --genomeFastaFiles ' + fa_path + ' --sjdbGTFfile ' + GTF_path + ' --runThreadN 10', shell=True)
    if index.returncode !=0:
        print('Error occured, please try again')
    else:
        print('Indexing done!')
        print('---------------')
        time.sleep(0.5)    

def main():
    args = parser.parse_args()
    SRA_download(args)
    Quality_control(args)
    Trimming(args)
    Multiqc(args)
    STAR_files_fasta(args)
    STAR_files_GTF(args)
    Unzip(args)

if __name__ == '__main__':
    main()