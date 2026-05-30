#!/bin/user/python3
import subprocess
import os 
import time
import zipfile
import glob 
import shutil
import pandas as pd 



#Begigning of teh pipline 
#This is pipline was done to make a data prepartion for R analysis of bidirectional promoters using RNA-seq data 
print('Begining the data preparation for R analysis')
print('--------------------------------------------')
time.sleep(0.5)
print('Please insure to do a configuration for SRA before begining')
time.sleep(0.5)
#Craeting directory 
directory = 'sra'
print('Creating a directory called:', directory)
print('---------------------------------------')
time.sleep(0.5)
print('Files from SRA would be sored in here')
print('--------------------')
time.sleep(0.5)


#Function to download fastq files from SRA 
def SRA_download(): 
    try:
        os.mkdir(directory)
    except OSError as error:
        print('the error:', error)
        time.sleep(0.5)
    print('Directory named sra created')
    print('-----------------')
    os.chdir(directory)
    while True:
        question1 = input('Please provide SRR number of te study: ')
    
        print('Begining to download runs from sra in sra directory')
        print('----------------------------------------------------')
        download = subprocess.run("esearch -db sra -query " +  question1 + " | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | xargs fastq-dump  --skip-technical  --readids --read-filter pass --dumpbase --split-3", shell=True)
        if download.returncode !=0:
            print('Error occured')
            question2 = input('Would you like to try again? yes/no: ')
            if question2 == 'no':
                print('Exiting pipline....')
                exit()
            if question2 == 'yes':
                print('Okay, retrying....')
        else: 
            print('Download completed')
            print('------------------')
            time.sleep(0.5)
            break



#SRA_download()

#Functiojn to perform qualuity control 
def Quality_control():

    print('Begigning performing quality control...')
    print('---------------------------------')
    time.sleep(0.5)
    print('Creating new directory for FASTQC files...')
    print('----------------------------------')
    time.sleep(0.5)

    directory2 = 'fastqc'
    print('Directory created for fastqc files:', directory2)
    

    try:
        os.mkdir(directory2)
    except OSError as error:
        print('the error:', error)
    time.sleep(0.5)
    print('Directory ' +  fastqc + ' created')
    print('--------')


#running FastQC

    print('Running FastQC....')
    time.sleep(0.5)
    print('---------------')
    fastqc_run = subprocess.run('fastqc * -o fastqc', shell=True)
    if fastqc_run.returncode !=0:
        print('Errror occured')
    else:
        print('Finished FastQC analysis')
        print('-------------------------')
        time(0.5)
        print('All the files are in the fastqc directory')
        print('---------------------------------------')

#Quality_control()

def Trimming():

#trimming opton 
    trim_q = input('Would you like to do trimming? yes/no: ')
    if trim_q == 'no':
        print('Proceeding without trimming...')
        return
    

    if trim_q == 'yes':
        directoryt = 'trim'
        try:
            os.mkdir(directoryt)
        except OSError as error:
            print('the error:', error)
            time.sleep(0.5)
            print('Directory ' +  directoryt + ' created')
            print('--------')
            os.chdir(directoryt)
            print('To trim the cutadapt is used in this pipline')
            time.sleep(0.5)
            print('Please note here the simple trimming was done')
            print('If want to perform trimming with linked adapetr will have to chenge the script manually')
            time.sleep(0.5)
        
            adapter1 = input('Please provide adapter sequence 1: ')
            adapter2 = input('Please provide adapter sequence 2: ')

            if not os.path.exists('sra'):
                print('The sra directory does not exist. Exiting.')
                return
            
            
            os.chdir(directoryt)
            for file in os.listdir('sra'):
                if file.endswith('_1.fastq') or file.endswith('_1.fastq.gz'):
                    input_file = os.path.join('sra', file)
                    input_file_2 = os.path.join('sra', file.replace('_1.fastq', '_2.fastq').replace('_1.fastq.gz', '_2.fastq.gz'))
                    output_file = os.path.join('trim', file.replace('.fastq', '_trimmed.fastq').replace('.fastq.gz', '_trimmed.fastq.gz'))
                    output_file_2 = os.path.join('trim', file.replace('_1.fastq', '_2_trimmed.fastq').replace('_1.fastq.gz', '_2_trimmed.fastq.gz'))

                    cutadapt = subprocess.run('cutadapt -a ' + adapter1 + ' -A ' + adapter2 + ' -o ' + output_file + ' -p ' + output_file_2 + '' + input_file + '' + input_file_2, shell=True)
                    if cutadapt.returncode !=0:
                        print('Error occured')
                    else: 
                        print('The trimming is now done..')
                        print('-----------------')
                        os.chdir('..')
                    
#Trimming()

def Multiqc():



    mutiqc_q = input('Would you like to perfrom MultiQC? yes/no: ')
    if mutiqc_q == 'no':
        print('Ok, proceeding further without MULTIQC...')
    if mutiqc_q == 'yes':
        print('Performing MultiQC....')
        print('Begining MultiQC analysis')
        time.sleep(0.5)
    print('--------------------------')
 
    multiqc_run = subprocess.run('multiqc .', shell = True)
    if multiqc_run.returncode !=0:
        print('Error ocuured...')
    else:
        print('-----------------------------')
        print('Finished MultiQC')
        time.sleep(0.5)
        print('Multiqc direcotry created')
        print('Can find results of analysis in there')
        time.sleep(0.5)
        print('---------------------------')

#Multiqc()

def STAR_files_fasta():
    os.chdir('..')
    print('Prior to aligment need to get the fasta and GTF files from needed genome')
    print('e.g. Mus musculus')
    time.sleep(0.5)
    print('To download those files please provide links to them from Ensembl website')
    while True:
        time.sleep(0.5)
        print('Example of link for fa file: https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz ')
        fasta = input('Please provide a link to a fasta file: ')
        print('Downloading assembly for Mus musculus')
        print('--------------------------')
        fatsa_download = subprocess.run('wget ' +  fasta, shell = True)
        if fatsa_download.returncode !=0:
            print('Error occured, please try again')
        else:
            print('Download succesfull!')
            print('-----------------------------')
            time.sleep(0.5)
            break 
#STAR_files_fasta()

def STAR_files_GTF():
    while True: 
        time.sleep(0.5)
        print('Example of link for fa file: https://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz ')
        GTF = input('Please provide a link to GTF file: ')
        print('Downloading GTF for Mus musculus')
        print('--------------------------')
        GTF_download = subprocess.run(' wget ' + GTF, shell = True)
        if GTF_download.returncode !=0:
            print('Error occured, please try again ')
        else: 
            print('Download succesfull!')
            print('-----------------------------')
            time.sleep(0.5)
            break 
#STAR_files_GTF()

def Unzip(): 
    time.sleep(0.5)
    gz_files = glob.glob('*.gz')
    if not gz_files:
        print('Error no files not find in directory')
        print('Please try again')
        #STAR_files_fasta()
        #STAR_files_GTF()
    else:
        print('Begin to unzip files')
        print('--------------------')
        unzip = subprocess.run(['gunzip'] + gz_files, check=True)
        if unzip.returncode !=0:
            print('Error occured')
        else:
            print('Both files unziped')
            print('-------------------')
            time.sleep(0.5)
            files = os.listdir('.')
            print('Files in the current directory:')
            print('-------------------------------')
            for file in files:
                print(file)

#Unzip()

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
#Indexing()

def STAR_map():
    
    print('The next step is to perform mapping')
    time.sleep(0.5)
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

    print('Please provide the full path for directories sra and ref (was created after indexing)')
    time.sleep(0.5)

    while True: 
        print('Please provide path for sra directory below')
        print('Example: /home/s2614505/PROJECT/sra')
        time.sleep(0.5)
        sra_path = input('Path : ')
        if os.path.isdir(sra_path):
            print("Directory exists")
            print('Procedding...') 
            print('-----------')
            time.sleep(0.5)
            break
        else:
            print('Diretory not found, please try again')
    while True: 
        print('Please provide path for ref directory below')
        print('Example: /home/s2614505/PROJECT/ref')
        time.sleep(0.5)
        ref_path = input('Path: ')
        if os.path.isdir(ref_path):
            print('Directory exists')
            print('Roceeding..')
            print('-----------')
            time.sleep(0.5)
            break 
        else:
            print('Directory not found, please try again')
    
   
    print('To perform mapping a output directory will be cretaed: aligment')
    directory_al = 'aligmnet'
    try: 
        os.mkdir(directory_al)
    except OSError as error:
        print('the error:', error)
    time.sleep(0.5)
    print('Directory ' +  directory_al +  ' created')
    print('----------------------------------')
    time.sleep(0.5)
    
    while True: 
        os.chdir(home_path)
        print('Please provide path for aligmnet directory below')
        print('Example: /home/s2614505/PROJECT/aligmnet')
        time.sleep(0.5)
        al_path = input('Path: ')
        if os.path.isdir(al_path):
            print("Directory exists")
            print('Procedding...') 
            print('-----------')
            time.sleep(0.5)
            break
        else:
            print('Diretory not found, please try again')

    time.sleep(0.5)
    print('To perform mapping fastqc directory has to be moved from sra')
    while True: 
        print('Below please provide a path to fastqc directory')
        print('Example: /home/s2614505/PROJECT/fastqc')
        time.sleep(0.5)
        fastqc_path = input('Path: ')
        if os.path.isdir(fastqc_path):
            print('Directory exists')
            print('Roceeding..')
            print('-----------')
            time.sleep(0.5)
            break 
        else:
            print('Directory not found, please try again')

    
    try: 
        shutil.move(fastqc_path, home_path)
    except OSError as error:
        print('Directory already moved')

    while True:
        print('Last directory path')
        print('Please provide below new path for fastqc directory')
        print('Example: /home/s2614505/PROJECT/fastqc')
        time.sleep(0.5)
        fastqc_path_new = input('Path: ')
        if os.path.isdir(fastqc_path_new):
            print('Directory there')
            print('Can proceed with mapping now...')
            time.sleep(0.5)
            print('-----------------------')
            break
        else: 
            print('Directory not found please try again')
        
    print('Now getting files from fatqc directory')
    os.chdir(fastqc_path_new)
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
    print('Begining aligment')
    print('-----------------')
    os.chdir(home_path)


    #defining bases 
    for name in sorted(names):
        fq1 = os.path.join('sra', name + '_pass_1.fastq')
        fq2 = os.path.join('sra', name + '_pass_2.fastq')
        aligned_read = os.path.join(al_path, name)
        time.sleep(0.5)
        map = subprocess.run("STAR --runThreadN 10 --genomeDir " + ref_path +  " --readFilesIn " + fq1 + " " + fq2 + " " "--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix " + aligned_read + "_", shell=True)
        if map.returncode !=0: 
            print('Error occured')
        else: 
            print('Finised')  

#STAR_map()

def Bed_file_making():

    print('the aligment and mapping are now done')
    time.sleep(0.5)
    print('Below please provide a path to the bed file with the positios of promoters')
    time.sleep(0.5)
    print('The bed file can be taken from UCSC database')
    time.sleep(0.5)
    
    while True: 
        print('Exmaple for path: /home/s2614505/Diss/ann_to_m39.bed')
        path_bed = input('Please provide path for bed file here: ')
        if os.path.isfile(path_bed):
            print('File found: ', path_bed)
            path_correct = input('Is this file correct? yes/no: ')
            if path_correct == 'yes':
                print('Proceeding...')
                print('-------------')
                time.sleep(0.5)
                break 
            if path_correct =='no':
                print('Please try again ')
        else:
            print('File not found')
            time.sleep(0.5)
            print('Please try again')
    
    print('To modify bed file for firther anlysis, flankbed tool will be used')
    time.sleep(0.5)
    print('Need tp generate genome.sizes from fasta file')
    time.sleep(0.5)
    while True:
        print('Please provide path to the Mus_musculus.GRCm39.dna_sm.primary_assembly.fa')
        print('Example: /home/s2614505/Diss/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa')
        time.sleep(0.5)
        fasta_p = input('Path: ')
        if os.path.isfile(fasta_p):
            print('The file exists', fasta_p)
            print('------------------------')
            time.sleep(0.5)
            break
        else: 
            print('File not found, please try again')
    
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

    print('New file craeted: ', gen_size_new)

    #remove genome.size without chr 
    os.remove(gen_size)
    print('Old file removed')
    time.sleep(0.5)

    #using flank to modify bed 

    while True:
        bed_file_new = input('Please provide a name for a new BED file, must end with .bed (new.bed): ')
        if bed_file_new.endswith('.bed'):
            print('Proceeding....')
            break
        else: 
            print('Please try again.')
        
    # Step 6: Get the number of base pairs for the flanking region
    while True:
        number_ofbp = input('What bp region do you want to take? Please input number (e.g., 500): ')
        if number_ofbp.isdigit() and int(number_ofbp) > 0:
            number_ofbp = int(number_ofbp)  # Convert to integer
            print('Proceeding...')
            break
        else:
            print('Not a valid number, try again')

    # Step 7: Read the input BED file and create a new BED file with midpoints
   

    print('Reading input BED file and calculating midpoints')
    print('------------------------------------------------------------')

    columns_bed = ['chr', 'start', 'end', 'gene', 'score', 'strand', 'pos1', 'pos2']
    df = pd.read_csv(path_bed, sep='\t', names=columns_bed, header=None)

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

    flank_command = f'flankBed -i {midpoint_bed} -g {gen_size_new} -b {number_ofbp} > {flank_bed}'
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

        while True:
            gtf_file_name = input('Please provide a name for the GTF file (must end with .gtf): ')
            if gtf_file_name.endswith('.gtf'):
                print('Proceeding...')
                break
            else:
                print('Please try again.')

        genepred_to_gtf = subprocess.run(f'./genePredToGtf file {genepred} {gtf_file_name}', shell=True)
        if genepred_to_gtf.returncode != 0:
            print('Error occurred.')
            exit(1)
        else:
            print('GTF file created for featureCounts')
            print('----------------------------------')


        
Bed_file_making()

def featureCounts():


    print('Please provide path for aligmnet directory: /home/s2614505/Diss/aligmnet')
    aligmnet_path = input('Path: ')
    if os.path.isdir(aligmnet_path):
        print('Path exists!')
        print('Proceeding...')
    else: 
        print('Soemthing went wrong')
    
    print('Please provide path for bed file : /home/s2614505/Diss/my.gtf')
    gtffile_path = input('Path: ')
    if os.path.isfile(gtffile_path):
        print('Path exists!')
        print('Proceeding...')
    else: 
        print('Soemthing went wrong')
    
    print('Begining to perform featureCounts')
    print('-----------------------------------')

    feature = subprocess.run('featureCounts -s 2 -a ' + gtffile_path + ' -o counts.txt -T 10 -p ' + aligmnet_path + '/*bam', shell=True)
    if feature.returncode !=0:
        print('Error occured')
    else:
        print('Created counts.txt and counts.txt.summary')
        print('-------------------')
        time.sleep(0.5)

featureCounts()

def Final():
    final_q1 = input('Do you want to vizualise results in IGV? yes/no: ')
    if final_q1 == 'yes':
        print('Okay..')
        print('Proceeding....')
        time.sleep(0.5)
        print('Creating bedgraph files for vizualization at IGV')
        time.sleep(0.5)
        print('------------------------------------------------')

        counts_file = 'counts.txt'
        df = pd.read_csv(counts_file, sep='\t', comment = '#')

        #now need to separate counts by strands so can have two sep files for both strands  
        df_counts_pos = df[df['Strand'] == '+']
        df_counts_neg = df[df['Strand'] == '-']

        # craeting dataframe for + bedgraph 

        
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
        exit()
    else: 
        time.sleep(0.5)
        print('Now you can proceed with R analysis')
        time.sleep(0.5)
        print('Using Analysis.Rmakdown script')
        time.sleep(0.5)
        print('Exiting now......')
        time.sleep(0.5)
        exit()


Final()


