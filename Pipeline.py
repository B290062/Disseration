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
config.read("config2.ini")

#The following are command line arguments that can be specified in the console by the user, each of which denotes its use.
parser = argparse.ArgumentParser()
parser.add_argument("--sra", action="store", required=True, help="This is the SRA number from Gene Expression eg for GSE87821 it would be SRP091444)")
#store_true used for optional values, when the value is present it will be True so this can be used to write functions.
#https://docs.python.org/3/library/argparse.html#action
parser.add_argument("--trim", action="store_true", required=False, help = "The data can be optionally trimmed if the user requires")
parser.add_argument("--adapter1",required=False,help="The first adapter required for trimming. GRO-seq only requires one adapter.")
parser.add_argument("--adapter2",required=False, help="The second adapter required for trimming") 
parser.add_argument("--multiqc", action="store_true", required=False, help = "This argument enables the use of MultiQC")
#Function to download fastq files from SRA 
parser.add_argument("--fasta", required = False, default =config["dir"]["fasta"], help = "This is the file path for the fasta file")
parser.add_argument("--gtf",required = False, default=config["dir"]["gtf"] , help = "This is the file path for the gtf file")
parser.add_argument("--mode", choices= ["rnaseq", "groseq"], required = True, help= "The type of data the user wants to analyse, either " \
"RNA-Seq or GRO-seq data")
parser.add_argument("--mask", required = True, help="This is the BED file that has coordinates of the region of interest e.g promoters, enhancers.")
#--window works as the bed file (mm39_refseq.bed) was downloaded as upstream by one base, so the window can be configured by the user as desired.
parser.add_argument("--window", type= int, default = 500, help="Number of bases flanking the promoter/enchancer region")

def SRA_download(args): 
    #replaced with makedirs instead of os.mkdirs as it has the exist_ok function which prevents crashing
    #when the folder is already created.
    os.makedirs("SRA", exist_ok=True)
    #moves into the created/existing SRA directory
    os.chdir("SRA")
    print('Attempting to download fastq files from GEO')
    #if statement adapted from https://www.geeksforgeeks.org/python/check-if-directory-contains-file-using-python/
    #this checks if the current directory "SRA" contains one or more files and skips the download process if true
    #essential as the SRA files for any given study are very large.
    if len(os.listdir(".")) > 0:
        print("SRA files were detected, inside the SRA directory, therefore the download will be skipped")
        os.chdir("..")
        return
    #this search function looks at at the SRA database, however the search command could be adapted for GEO number...
    download = subprocess.run("esearch -db sra -query " +  args.sra + " | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | xargs fastq-dump  --skip-technical  --readids --read-filter pass --dumpbase --split-3", shell=True)
    if download.returncode !=0:
        print('Error occured during download of fastQ files')
        exit(1)
    else: 
        print('Fastq file download successful, storing in the SRA directory.')
        #leaves the SRA directory backs into the main directory.
        os.chdir("..")

def Quality_control(args):
    #produces FastQC files which allow the user to examine to determine the quality of the data.
    print('Beginning performing quality control...')
    os.makedirs("FastQC", exist_ok = True)
    print("FastQC directory created")
    #this code appears in all of the functions and allows for time to be saved, if files are detected.
    #only portions of the pipeline that need to be ran are ran.
    if len(os.listdir("FastQC")) > 0:
        print("FastQC files detected in the folder, skipping this step")
        return
    #runs fastqc on all the files in the SRA folder
    fastqc_run = subprocess.run('fastqc SRA/*.fastq -o FastQC', shell=True)
    if fastqc_run.returncode !=0:
        print('FastQC were not able to be produced')
        exit(1)
    else:
        print('Finished FastQC analysis')
        print('All the files are in the fastqc directory')

def Trimming(args):
    #trimming option Note- this is for simple trimming only, not linked adapter trimming
    
    #--trim needs to be specified in the command line. Otherwise it skips this step entirely
    if args.trim is False:
        print('Proceeding without trimming...')
        return
    #with mode gro-seq this code checks if one adapter was specified by the user, if not the code exits.
    #note - if the user enters two adapters with gro-seq, no error will occur but only the first adapter will be used.
    if args.trim is True and args.mode =="groseq" and args.adapter1 is None:
        print("Please provide an adapter for trimming")
        exit(1)
    #with mode rna-seq two adapters are needed so the code checks that both have been inputted, if not it exits.
    if args.trim is True and args.mode =="rnaseq" and (args.adapter1 is None or args.adapter2 is None):
        print("Please provide both adapters for trimming")
        exit(1)

    if args.trim is True:
        os.makedirs("Trimmed_data", exist_ok=True)
        if not os.path.exists('SRA'):
            print('The SRA directory does not exist. Exiting.')
            exit(1)
        #The structure of SRA files are different between RNA-Seq and GRO-Seq. RNAS has two passes, while gro-seq only has one.
        #so the files need to be trimmed differently to account for this.    
        if args.mode == "rnaseq":
            for file in os.listdir('SRA'):
                #removed the or file.endswith("_pass_1.fastq.gz as none of the files in SRA are processed like this")
                if file.endswith("_pass_1.fastq"):
                    #the path of the SRA, folder and the file are joined together. E.G /home/diss/SRA + SRR441391_pass_1.fastq
                    input_file = os.path.join('SRA', file)
                    #this replaces the SRA file with the complimentory SRA file with _pass_2
                    input_file_2 = os.path.join('SRA', file.replace('_pass_1.fastq', '_pass_2.fastq'))
                    # The output file names originally had _trimmed in the name, but as they go into the /Trimmed_data folder
                    # this was removed, as it map the star map function work correctly with trimming.
                    output_file = os.path.join('Trimmed_data', file)
                    output_file_2 = os.path.join('Trimmed_data', file.replace('_pass_1.fastq', '_pass_2.fastq'))

                    #cutadapt command that takes the two inputs and adapters and performs trimming
                    cutadapt = subprocess.run('cutadapt -a ' + args.adapter1 + ' -A ' + args.adapter2 + ' -o ' + output_file + ' -p ' + output_file_2 + '' + input_file + '' + input_file_2, shell=True)
                    if cutadapt.returncode !=0:
                        print('Error occured')
                        exit(1)
                    else: 
                        print('The trimming is now done..')
        else:
            for file in os.listdir("SRA"):
                if file.endswith("_pass.fastq"):
                    input_file = os.path.join('SRA', file)
                    output_file = os.path.join('Trimmed_data', file)
                    #the cutadapt command is shortened for gro-seq as it only requires one input, output and adapter. The -A for second adapter
                    # and -p "paired" flags are removed.
                    cutadapt = subprocess.run('cutadapt -a ' + args.adapter1 + ' -o ' + output_file + ' ' + input_file , shell=True)
                    if cutadapt.returncode !=0:
                        print('Error occured')
                        exit(1)
                    else: 
                        print('The trimming is now done..')
                        
def Multiqc(args):
    #this method takes all of the FastQC files and combines them into a MultiQC for easier interpretation

    #if argument not specified, it returns and doesn't run
    if args.multiqc is False:
        print('MultiQC was not performed')
        return
    #if argument is specified, it runs and results are put into multiqc folder
    if args.multiqc is True:
        print('Performing MultiQC....')
        multiqc_run = subprocess.run('multiqc .', shell = True)
        if multiqc_run.returncode !=0:
            print('Error occured while running MultiQC')
            exit(1)
        else:
            print('MultiQC directory was produced containing the results')                  

def STAR_files_fasta(args):
    #this function enables the download of the STAR FASTA genome files
    #checks if .gz files already exist in the directory to skip the download
    #this was adapted from the previous code from the unzip function.
    gz_files = glob.glob('*.fa.gz')
    if len(gz_files) > 0:
        print("Fasta detected. skipping the download")
        return
    #downloads the fasta file based on the link specified in the config file
    print(f'Downloading assembly {config["urls"]["fasta_url"]}')
    fasta_download = subprocess.run('wget ' +  config['urls']["fasta_url"], shell = True)
    if fasta_download.returncode !=0:
        print('Error occured while attempting to download FASTA file...')
        exit(1)
    else:
        print('FASTA download successful')
                
def STAR_files_GTF(args):
    # function provides the annotation coordinates

    #skips download if the gtf files are already downloaded.
    gz_files = glob.glob('*.gtf.gz')
    if len(gz_files) > 0:
        print("GFT files detected. skipping the download")
        return
    #downloads gtf based on the link that is specified in the config file
    print(f'Downloading GTF {config["urls"]["gtf_url"]}')
    GTF_download = subprocess.run(' wget ' + config["urls"]["gtf_url"], shell = True)
    if GTF_download.returncode !=0:
        print('Error occured while attempting to download GTF file...')
        exit(1)
    else: 
        print('GTF download successful')
        
def Unzip(args):
    # STAR files and compressed when downloaded and therefore need to be unzipped. 

    #these variables are used to check if unzipped files are present allowing this step to be skipped.
    gtf_files = glob.glob('*.gtf')
    fa_files = glob.glob('*.fa')
    
    #skip if uncompressed file are detected
    if len(gtf_files) > 0 and len(fa_files) > 0:
        print("Unzipped file already detected, skipping.")
        return
    #the downloaded fasta and gtf files end with .gz and if they arent present unzipping cannot occur
    gz_files = glob.glob('*.gz')
    if not gz_files:
        print('No files detected to unzip.')
        exit(1)
    else:
        print('Begin to unzip files')
        #unzips the files ending with .gz
        unzip = subprocess.run(['gunzip'] + gz_files)
        if unzip.returncode !=0:
            print('The zipped files were unable to be unzipped.')
            exit(1)
        else:
            print('Both files unzipped')
            
def Indexing(args):
    print('Index is being built')
    #this is the path to the fasta file, if the file is already downloaded if can be specified in the command line.
    #the default path has been set up for this project 
    if os.path.isfile(args.fasta):
        print('The file exists', args.fasta)
        
    else: 
        print('Fasta file not found')
        exit(1)
    #the same applies for args.gtf.
    if os.path.isfile(args.gtf):
        print('The file exists', args.gtf)
        
    else:
        print('GTF file not found, please try again')
        exit(1)
    
    #skip indexing if files are detected in the ref folder.
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

def STAR_map(args):
     #main adaptation of this function involves mapping GRO-seq data too, I investigated how other pipelines do this
    # The difference in processing is that GRO-seq is single end sequencing while RNA-Seq is paired end.
    # https://github.com/Danko-Lab/proseq2.0/blob/master/proseq2.0.bsh
    print('Performing alignment')
    #this function was especially cluttered in the previous code, the issue of the FASTQC file directory being a subdirectory of
    # SRA was fixed earlier so many of the inputs previously present here are not relevent
    #  also many of the inputs appeared irrelevent so they were removed.

    #checks if alignment has already be done and skips if so.   
    if os.path.exists('Alignment') and len(os.listdir("Alignment")) > 1:
        print("Alignment has already been produced... Skipping.")
        return
    os.makedirs("Alignment", exist_ok=True)
    print("Alignment directory was created.")
    #if the data was trimmed the directory is set to trimmed data to use this data in the alignment
    if args.trim == True:
        search_directory = "./Trimmed_data"
    #otherwise it will use the SRA data
    else:
        search_directory = "./SRA"
    #creates a list of all the SRA files and sorts them in the variable find_base
    find_base = subprocess.run(f"find {search_directory} -name 'SRR*' -print| sort | uniq", shell=True, capture_output= True, text=True)
    if find_base.returncode !=0:
        print('Error occured')
    else:
        print('Reads found:')
        names = set()
        #this process extracts only the SRA name from the filepath and the _pass
        for line in find_base.stdout.splitlines():
            #splits on the file path
            name = line.strip().split('/')[-1]
            #splits on the _ to remove the _pass.fastq
            name = name.split('_')[0]
            names.add(name)
        for name in sorted(names):
            print(name)

    if args.mode == "rnaseq":
        for name in sorted(names):
            #adds the passes back onto the the SRA names to allow them to be distinguished from eachother in the mapping
            fq1 = os.path.join(search_directory, name + '_pass_1.fastq')
            fq2 = os.path.join(search_directory, name + '_pass_2.fastq')
            aligned_read = os.path.join("Alignment", name)
            map = subprocess.run("STAR --runThreadN 10 --genomeDir " + "ref" +  " --readFilesIn " + fq1 + " " + fq2 + " " "--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix " + aligned_read + "_", shell=True)
            if map.returncode !=0: 
                print('Error occured during mapping')
                exit(1)
            else: 
                print('Finished mapping')  
    
    #if the mode is not RNA-Seq it must be GRO-Seq which, is single ended, therefore only Fq1 is needed.
    else:
        for name in sorted(names):
            fq1 = os.path.join(search_directory, name + "_pass.fastq")
            aligned_read = os.path.join("Alignment", name)
            map = subprocess.run("STAR --runThreadN 10 --genomeDir " + "ref" +  " --readFilesIn " + fq1 + "  --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix " + aligned_read + "_", shell=True)
            if map.returncode !=0: 
                print('Error occured during mapping')
                exit(1)
            else: 
                print('Finished mapping')  
            
def Bed_file_making(args):
    #the final output of bed file modification is the bed file in gtf format used in featurecounts
    #therefore the pipeline can check if this exists and skip this step.
    gtf_file_check = args.mask.replace(".bed", ".gtf")
    if os.path.isfile(gtf_file_check):
        print(f'{gtf_file_check} Already exists, skipping the BED file modification')
        return

    if os.path.isfile(args.mask):
        print('File found: ', args.mask)
    else:
        print('BED file not found')
        exit(1)
    
    print('To modify bed file for further analysis, flankbed tool will be used')
    print('Need to generate genome.sizes from fasta file')
    fasta_p = args.fasta
    if os.path.isfile(fasta_p):
        print('The file exists', fasta_p)
    else: 
        print('File not found, please try again')
        exit(1)
    
    print('Modifying BED file for further analysis')
    #this line converts the FASTA format to fai, essentially it extracts the info that was generated with the FASTA file
    # it is used to get precise file coordinates. https://pmc.ncbi.nlm.nih.gov/articles/PMC8558547/
    fasta_to_fai = subprocess.run('samtools faidx ' + fasta_p, shell=True)
    if fasta_to_fai.returncode !=0:
        print('Error occured')
    else: 
        print('Fai file generated: ', fasta_p + '.fai')
        new_fai = fasta_p + '.fai'
    
    print('By using fai can now create a genome.sizes file that is required for fblank option')
    gen_size = 'genome.size_now'
    print('A new file generated: ', gen_size)
    #takes only the necessary columns from the fai file it has 5 columns but it filters to keep the chromosome and the length.
    fai_to_gen = subprocess.run("awk '{{print $1\"\t\"$2}}' " + new_fai + " > " + gen_size, shell=True)
    if fai_to_gen.returncode !=0:
        print('Error occured')
    else: 
        print('Proceeding to customizing genome.size_now for flankbed')
    
    #genomesize contains only the chromosome and the length
    gen_size_new = 'genome.size'

    #this loop puts the chr prefix infront of the chromosome numbers in the genomesize file
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

    #using flank to modify bed 
    bed_file_new = args.mask.replace(".bed", "_flanked.bed")
        
    # Step 7: Read the input BED file and create a new BED file with midpoints
    print('Reading input BED file and calculating midpoints')
    #bed6 file format downloaded from UCSC
    columns_bed = ['Chromosome', 'Start', 'End', 'Gene', 'Score', 'Strand']
    df = pd.read_csv(args.mask, sep='\t', names=columns_bed, header=None)
    #this code below makes a list of valid chromosomes that work with flankbed, it removed the errors of chrM not found and other chromosomes.
    target_chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
                          "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
                          "chr18", "chr19", "chrX", "chrY" ]
    df= df[df["Chromosome"].isin(target_chromosomes)]

    # Calculate midpoints
    df['midpoint'] = (df['Start'] + df['End']) // 2

    # Create a new DataFrame for the midpoint BED file
    df_midpoints = df[['Chromosome', 'midpoint', 'midpoint', 'Gene', 'Score' ,'Strand']]

    # Save the midpoint BED file
    midpoint_bed = 'midpoints.bed'
    df_midpoints.to_csv(midpoint_bed, sep='\t', header=False, index=False)

    print('Midpoint BED file created:', midpoint_bed)

    # Step 8: Use flankBed to generate flanking intervals around the midpoints
    print('Generating flanking intervals using flankBed')
    flank_bed = bed_file_new  # The final BED file name as specified by the user

    #use the window argument to extend the midpoints by the specified window
    flank_command = f'flankBed -i {midpoint_bed} -g {gen_size_new} -b {args.window} > {flank_bed}'
    bed_to_new = subprocess.run(flank_command, shell=True)
    
    if bed_to_new.returncode != 0:
        print('Error occurred while running flankBed.')
        exit(1)
    else:
        print('flankBed performed successfully.')
        print('Flanking intervals BED file created:', flank_bed)


    # Step 9: Modify the BED file to contain both strands for divergent transcription analysis
    print('Now modifying BED file further to contain both strands to identify divergent transcription')

    df = pd.read_csv(flank_bed, sep='\t', names=columns_bed, header=None)

    # Create unique identifiers for each region
    gene_counter = {}
    def get_unique_id(gene):
        if gene not in gene_counter:
            gene_counter[gene] = 1
        else:
            gene_counter[gene] += 1
        return f"{gene}_{gene_counter[gene]}"
    df['Gene'] = df['Gene'].apply(get_unique_id)

    df_plus = df.copy()
    df_plus['Strand'] = '+'

    df_minus = df.copy()
    df_minus['Strand'] = '-'

    df_strands = pd.concat([df_plus, df_minus])

    df_strands.to_csv(flank_bed, sep='\t', header=False, index=False)
    print('BED file now fully modified')
    
    # Step 10: Convert the BED file to GTF format for featureCounts
    print('To perform featureCounts ' + flank_bed + ' needs to be converted to GTF file')
    print('Will create genepred file')
    genepred = 'bed.genepred'
    bed_to_genepred = subprocess.run(f'./bedToGenePred {flank_bed} {genepred}', shell=True)
    if bed_to_genepred.returncode != 0:
        print('Error occurred.')
        exit(1)
    else:
        print('Created genepred file')
    

        
    gtf_file_name = args.mask.replace(".bed", ".gtf")

    genepred_to_gtf = subprocess.run(f'./genePredToGtf file {genepred} {gtf_file_name}', shell=True)
    if genepred_to_gtf.returncode != 0:
        print('Error occurred.')
        exit(1)
    else:
        print('GTF file created for featureCounts')

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
    if args.mode == "rnaseq":
        feature = subprocess.run('featureCounts -s 2 -a ' + gtf_file_name + ' -o counts.txt -T 10 -p ' + alignment_path + '/*.bam', shell=True)
        if feature.returncode !=0:
            print('Error occured')
            exit(1)
        else:
            print('Created counts.txt and counts.txt.summary')

    #Feature counts has the -p flag which corresponds to counting paired end data, this is removed for the Gro-seq as its single ended
    else:
        feature = subprocess.run('featureCounts -s 2 -a ' + gtf_file_name + ' -o counts.txt -T 10 ' + alignment_path + '/*.bam', shell=True)
        if feature.returncode !=0:
            print('Error occured')
            exit(1)
        else:
            print('Created counts.txt and counts.txt.summary')
        
           

def Final(args):
    print('Creating bedgraph files for vizualization with IGV')
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
    print('Now you can proceed with R analysis')
    print('Exiting now......')
    
    

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
    Final(args)

if __name__ == '__main__':
    main()