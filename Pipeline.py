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
parser.add_argument("--multiqc", action="store_true", required=False, help = "This argument enables the use of MultiQC")
#Function to download fastq files from SRA 
parser.add_argument("--fasta", required = False, default =config["dir"]["fasta"], help = "This is the file path for the fasta file")
parser.add_argument("--gtf",required = False, default=config["dir"]["gtf"] , help = "This is the file path for the gtf file")
parser.add_argument("--mode", choices= ["rnaseq", "groseq"], required = True, help= "The type of data the user wants to analyse, either " \
"RNA-Seq or GRO-seq data")
parser.add_argument("--mask", required = True, help="This is the BED file that has coordinates of the region of interest e.g promoters, enhancers.")
parser.add_argument("--window", type= int, default = 500, help="Number of bases flanking the promoter/enchancer region")

def SRA_download(args): 
    #replaced with makedirs instead of os.mkdirs as it has the exist_ok function which prevents crashing
    #when the folder is already created.
    os.makedirs("SRA", exist_ok=True)
    print('A Directory called SRA was created.')
    print('-----------------')
    print('Attempting to download fastq files from GEO')
    print('----------------------------------------------------')
    #if statement adapted from https://www.geeksforgeeks.org/python/check-if-directory-contains-file-using-python/
    #this checks if the current directory "SRA" contains one or more files and skips the download process if true
    #essential as the SRA files for any given study are very large.
    #this could also be replaced with a download skip argument that runs main without the SRA, FASTA, and GTF download steps. implement later.
    if len(os.listdir("SRA")) > 0:
        print("SRA files were detected, inside the SRA directory, therefore the download will be skipped")
        return
    #this search function looks at at the SRA database, however the search command could be adapted for GEO number...
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
    #runs fastqc on all the files in the SRA folder
    fastqc_run = subprocess.run('fastqc * -o FastQC', shell=True)
    if fastqc_run.returncode !=0:
        print('FastQC were not able to be produced')
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
            print('The SRA directory does not exist. Exiting.')
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
        print("Fasta detected. skipping the download")
        return
    print(f'Downloading assembly {config["urls"]["fasta_url"]}')
    print('--------------------------')
    fatsa_download = subprocess.run('wget ' +  config['urls']["fasta_url"], shell = True)
    if fatsa_download.returncode !=0:
        print('Error occured while attempting to download FASTA file...')
        exit(1)
    else:
        print('FASTA download successful')
        print('-----------------------------')
        time.sleep(0.5)
        
def STAR_files_GTF(args):
    # function provides the annotation coordinates
    gz_files = glob.glob('*.gz')
    if len(gz_files) > 0:
        print("GFT files detected. skipping the download")
        return
    print(f'Downloading GTF {config["urls"]["gtf_url"]}')
    print('--------------------------')
    GTF_download = subprocess.run(' wget ' + config["urls"]["gtf_url"], shell = True)
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
            print('Both files unzipped')
            print('-------------------')
            time.sleep(0.5)
            files = os.listdir('.')
            print('Files in the current directory:')
            print('-------------------------------')
            for file in files:
                print(file)

def Indexing(args):
    print('Index is being built')
    time.sleep(0.5)
    if os.path.isfile(args.fasta):
        print('The file exists', args.fasta)
        time.sleep(0.5)
        
    else: 
        print('Fasta file not found')
        exit(1)
    
    if os.path.isfile(args.gtf):
        print('The file exists', args.gtf)
        time.sleep(0.5)
        
    else:
        print('GTF file not found, please try again')
        exit(1)
    
    if os.path.exists('ref') and len(os.listdir("ref")) > 1:
        print("the index is already present so STAR is skipped")
        return
    os.makedirs("ref", exist_ok=True)
    #runs STAR alignment with the fasta and gtf file to produce a reference
    index = subprocess.run('STAR --runMode genomeGenerate --genomeDir ' + "ref" + '/ --genomeFastaFiles ' + args.fasta + ' --sjdbGTFfile ' + args.gtf + ' --runThreadN 10', shell=True)
    if index.returncode !=0:
        print('Index error occured')
        exit(1)
    else:
        print('Indexing done!')
        print('---------------')
        time.sleep(0.5)    

def STAR_map(args):
     #main adaptation of this function involves mapping GRO-seq data too, I investigated how other pipelines do this
    # The difference in processing is that GRO-seq is single end sequencing while RNA-Seq is paired end.
    # https://github.com/Danko-Lab/proseq2.0/blob/master/proseq2.0.bsh
    print('The next step is to perform mapping')
    #this function was especially cluttered in the previous code, the issue of the FASTQC file directory being a subdirectory of
    # SRA was fixed earlier so many of the inputs previously present here are not relevent
    #  also many of the inputs appeared irrelevent so they were re   
    os.makedirs("Alignment", exist_ok=True)
    print("Alignment directory was created.")
    
    find_base = subprocess.run("find . -name 'SRR*' -print| sort | uniq", shell=True, capture_output= True, text=True)
    if find_base.returncode !=0:
        print('Error occured')
    else:
        time.sleep(0.5)
        print('-----------------')
        print('Reads found:')
        time.sleep(0.5)
        names = set()
        for line in find_base.stdout.splitlines():
            name = line.strip().split('/')[-1]
            name = name.split('_')[0]
            names.add(name)
        for name in sorted(names):
            print(name)
    time.sleep(0.7)
    print('Beginning aligment')
    print('-----------------')

    #defining bases
    if args.mode == "rnaseq":
        for name in sorted(names):
            fq1 = os.path.join('SRA', name + '_pass_1.fastq')
            fq2 = os.path.join('SRA', name + '_pass_2.fastq')
            aligned_read = os.path.join("Alignment", name)
            time.sleep(0.5)
            map = subprocess.run("STAR --runThreadN 10 --genomeDir " + "ref" +  " --readFilesIn " + fq1 + " " + fq2 + " " "--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix " + aligned_read + "_", shell=True)
            if map.returncode !=0: 
                print('Error occured during mapping')
                exit(1)
            else: 
                print('Finished mapping')  
    
    #if the mode is not RNA-Seq it must be GRO-Seq which, is single ended, therefore only Fq1 is needed.
    else:
        for name in sorted(names):
            fq1 = os.path.join("SRA", name + "_pass_1.fastq")
            aligned_read = os.path.join("Alignment", name)
            time.sleep(0.5)
            map = subprocess.run("STAR --runThreadN 10 --genomeDir " + "ref" +  " --readFilesIn " + fq1 + "  --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix " + aligned_read + "_", shell=True)
            if map.returncode !=0: 
                print('Error occured during mapping')
                exit(1)
            else: 
                print('Finished mapping')  
            
def Bed_file_making(args):
    #this is the bed file, that is used as coordinates
    if os.path.isfile(args.mask):
        print('File found: ', args.mask)
    else:
        print('File not found')
        time.sleep(0.5)
        print('Please try again')
        exit(1)
    
    print('To modify bed file for further analysis, flankbed tool will be used')
    time.sleep(0.5)
    print('Need to generate genome.sizes from fasta file')
    time.sleep(0.5)
    fasta_p = args.fasta
    if os.path.isfile(fasta_p):
        print('The file exists', fasta_p)
        print('------------------------')
        time.sleep(0.5)
    else: 
        print('File not found, please try again')
        exit(1)
    
    print('Modifying BED file for further analyis')
    print('--------------------------------------')
    time.sleep(0.5)
    fasta_to_fai = subprocess.run('samtools faidx ' + fasta_p, shell=True)
    time.sleep(0.5)
    if fasta_to_fai.returncode !=0:
        print('Error occured')
    else: 
        print('Fai file generated: ', fasta_p + '.fai')
        print('------------------') 
        new_fai = fasta_p + '.fai'
    
    print('By using fai can now create a genome.sizes file that is required for fblank option')
    gen_size = 'genome.size_now'
    time.sleep(0.5)
    print('A new file generated: ', gen_size)
    time.sleep(0.5)
    fai_to_gen = subprocess.run("awk '{{print $1\"\t\"$2}}' " + new_fai + " > " + gen_size, shell=True)
    time.sleep(0.5)
    if fai_to_gen.returncode !=0:
        print('Error occured')
    else: 
        print('Proceeding to customizing genome.size_now for flankbed')
        print('-------------------------')
        time.sleep(0.5)
    
    gen_size_new = 'genome.size'

    with open(gen_size, 'r') as infile, open(gen_size_new, 'w') as uotfile:
        for line in infile:
            new_line = line.strip()
            if not new_line.startswith('chr'):
                new_line = 'chr' + new_line
                uotfile.write(new_line +'\n')

    print('New file created: ', gen_size_new)

    #remove genome.size without chr 
    os.remove(gen_size)
    print('Old file removed')
    time.sleep(0.5)

    #using flank to modify bed 
    bed_file_new = args.mask.replace(".bed", "_flanked.bed")
        
    # Step 7: Read the input BED file and create a new BED file with midpoints
    print('Reading input BED file and calculating midpoints')
    print('------------------------------------------------------------')

    columns_bed = ['chr', 'start', 'end', 'gene', 'score', 'strand', 'pos1', 'pos2']
    df = pd.read_csv(args.mask, sep='\t', names=columns_bed, header=None)

    # Calculate midpoints
    df['midpoint'] = (df['start'] + df['end']) // 2

    # Create a new DataFrame for the midpoint BED file
    df_midpoints = df[['chr', 'midpoint', 'midpoint', 'gene', 'score', 'strand', 'pos1', 'pos2']]

    # Save the midpoint BED file
    midpoint_bed = 'midpoints.bed'
    df_midpoints.to_csv(midpoint_bed, sep='\t', header=False, index=False)

    print('Midpoint BED file created:', midpoint_bed)

    # Step 8: Use flankBed to generate flanking intervals around the midpoints
    print('Generating flanking intervals using flankBed')
    flank_bed = bed_file_new  # The final BED file name as specified by the user

    flank_command = f'flankBed -i {midpoint_bed} -g {gen_size_new} -b {args.window} > {flank_bed}'
    bed_to_new = subprocess.run(flank_command, shell=True)
    time.sleep(0.5)
    if bed_to_new.returncode != 0:
        print('Error occurred while running flankBed.')
        exit(1)
    else:
        print('flankBed performed successfully.')
        print('Flanking intervals BED file created:', flank_bed)
        print('-----------------')
        time.sleep(0.5)

    # Step 9: Modify the BED file to contain both strands for divergent transcription analysis
    print('Now modifying BED file further to contain both strands to identify divergent transcription')
    print('---------------------------------------------------------------------------------------')

    df = pd.read_csv(flank_bed, sep='\t', names=columns_bed, header=None)

    # Create unique identifiers for each region
    gene_counter = {}
    def get_unique_id(gene):
        if gene not in gene_counter:
            gene_counter[gene] = 1
        else:
            gene_counter[gene] += 1
        return f"{gene}_{gene_counter[gene]}"
    df['gene'] = df['gene'].apply(get_unique_id)

    df_plus = df.copy()
    df_plus['strand'] = '+'

    df_minus = df.copy()
    df_minus['strand'] = '-'

    df_strands = pd.concat([df_plus, df_minus])

    df_strands.to_csv(flank_bed, sep='\t', header=False, index=False)
    print('BED file now fully modified')
    print('--------------------------------------------')
    time.sleep(0.5)

    # Step 10: Convert the BED file to GTF format for featureCounts
    print('To perform featureCounts ' + flank_bed + ' needs to be converted to GTF file')
    time.sleep(0.5)
    print('Will create genepred file')
    genepred = 'bed.genepred'
    bed_to_genepred = subprocess.run(f'./bedToGenePred {flank_bed} {genepred}', shell=True)
    if bed_to_genepred.returncode != 0:
        print('Error occurred.')
        exit(1)
    else:
        print('Created genepred file')
        print('---------------------')
        time.sleep(0.5)

        
    gtf_file_name = args.mask.replace(".bed", ".gtf")

    genepred_to_gtf = subprocess.run(f'./genePredToGtf file {genepred} {gtf_file_name}', shell=True)
    if genepred_to_gtf.returncode != 0:
        print('Error occurred.')
        exit(1)
    else:
        print('GTF file created for featureCounts')
        print('----------------------------------')

def featureCounts(args):
    alignment_path = "Alignment"
    if not os.path.isdir(alignment_path):
        print('Alignment path not detected')
        exit(1)
    
    gtf_file_name = args.mask.replace(".bed",".gtf")
    if not os.path.isfile(gtf_file_name):
        print("Bed file not found")
        exit(1)
    
    print('Beginning to perform featureCounts')
    print('-----------------------------------')
    if args.mode == "rnaseq":
        feature = subprocess.run('featureCounts -s 2 -a ' + gtf_file_name + ' -o counts.txt -T 10 -p ' + alignment_path + '/*bam', shell=True)
        if feature.returncode !=0:
            print('Error occured')
            exit(1)
        else:
            print('Created counts.txt and counts.txt.summary')
            print('-------------------')
            time.sleep(0.5)
    #Feature counts has the -p flag which corresponds to counting paired end data, this is removed for the Gro-seq as its single ended
    else:
        feature = subprocess.run('featureCounts -s 2 -a ' + gtf_file_name + ' -o counts.txt -T 10 ' + alignment_path + '/*bam', shell=True)
        if feature.returncode !=0:
            print('Error occured')
            exit(1)
        else:
            print('Created counts.txt and counts.txt.summary')
            print('-------------------')
            time.sleep(0.5)

def Final(args):
    print('Creating bedgraph files for vizualization with IGV')
    time.sleep(0.5)
    print('------------------------------------------------')
    counts_file = 'counts.txt'
    df = pd.read_csv(counts_file, sep='\t', comment = '#')

    #now need to separate counts by strands so can have two sep files for both strands  
    df_counts_pos = df[df['Strand'] == '+']
    df_counts_neg = df[df['Strand'] == '-']

    # creating dataframe for + bedgraph 
    df_bed_pos = pd.DataFrame()
    df_bed_pos['chr'] = df_counts_pos['Chr']
    df_bed_pos['start'] = df_counts_pos['Start'] -1 
    df_bed_pos['end'] = df_counts_pos['End']
    df_bed_pos['score'] = df_counts_pos.iloc[:, -1]

    df_bed_pos.to_csv('pos_for_mid.bedgraph', sep='\t', header=False, index=False)
    print('First bedgraph craeted: ', df_bed_pos)

    #creating dataframe for - bedgraph 
    df_bed_neg = pd.DataFrame()
    df_bed_neg['chr'] = df_counts_neg['Chr']
    df_bed_neg['start'] = df_counts_neg['Start'] -1 
    df_bed_neg['end'] = df_counts_neg['End']
    df_bed_neg['score'] = df_counts_neg.iloc[:, -1]

    df_bed_neg.to_csv('neg_for_mid.bedgraph', sep='\t', header=False, index=False )
    print('Second bedgraph created: ', df_bed_neg)

    print('Bedgraph files are created')
    print('--------------------------')
    time.sleep(0.5)
    print('Now you can proceed with R analysis')
    time.sleep(0.5)
    print('Using Analysis.Rmakdown script')
    time.sleep(0.5)
    print('Exiting now......')
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
    Indexing(args)
    STAR_map(args)
    Bed_file_making(args)
    featureCounts(args)

if __name__ == '__main__':
    main()