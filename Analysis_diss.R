library(biomaRt)
library(dplyr)
library(data.table)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(viridis)
library(GO.db)
library(VennDiagram)
library(grid)
library(EnhancedVolcano)
library(ggrepel)

library(org.Mm.eg.db)


#note the names used are ref_ids not the gene names for the catalogue so they need to be manually converted when showing the data in tables
#Download counts.txt
gene_data <- read.table("countsRNASEQ.txt", header=TRUE, sep="\t", comment.char="#")


# Set the Geneid column as row names
rownames(gene_data) <- gene_data$Geneid


# Sort data - create separate for counts and metadata: 
metdata <- gene_data[c(1:5)]
counts <- gene_data[-c(1:6)]


#apply a scale factor, which allows reads per million to be generated
scale_factor <- colSums(counts)
RPM <- sweep(counts, 2, scale_factor, FUN='/')*1e6

#this was changed to extract _1 and _2 rather than _2 or _3 like the previous study
pattern <- '.*_[12]$'
counts_filtered <- subset(RPM, grepl(pattern, rownames(RPM)))
#This table is prepared for DESeq2 analysis that does not require normalization: 
counts_new <- subset(counts, grepl(pattern, rownames(counts)))

#Please here divide the groups according to your study 
#Here 1:15 - samples columns with counts 
groups <- list('OCT4_untreated' = 1:15,
               'OCT4_treated ' = 16:30,
               'BRG1_untreated' = 31:45,
               'BRG1_treated' = 46:60
)
#After groups are defined manually now will go into function that will produce rankings automatically 
#Extracting indices - columns for each group 
groups <- lapply(groups, function(indices) counts_filtered[, indices])

#Calculating row sums for each group
row_sums <- lapply(groups, function(group) rowSums(group, na.rm = FALSE))

#Getting everything into a data frame
RPM_normal <- data.frame(Geneid_RPM = rownames(counts_filtered))
for (i in 1:length(row_sums)) {
  RPM_normal <- cbind(RPM_normal, row_sums[[i]])
  
}

#Renamimg the column names according to groups 
colnames(RPM_normal)[-1] <- names(groups)

#Checking the data 
head(RPM_normal)

#Now need to get Geneid_base
RPM_normal <- RPM_normal %>% 
  mutate(Geneid_base = gsub('_[fr]_[12]$','', Geneid_RPM)) 

# Creating list for rankings of each group 
ranking_tables <- list()

#Counts will be set to not be less then 10, change manually if need more/less
number <- 20 

for (group_name in names(groups)) {
  ratio_group <- RPM_normal %>%
    group_by(Geneid_base) %>%
    summarize(
      smallest_value = min(.data[[group_name]], na.rm = TRUE),
      largest_value = max(.data[[group_name]], na.rm = TRUE),
      largest_ratio = ifelse(smallest_value > 0, largest_value / smallest_value, NA),
      rank = ifelse(smallest_value >= number, smallest_value * (1 / largest_ratio), NA)
    ) %>% 
    filter(smallest_value >= number)
  
  
  #Order by rank in descending order
  newdata_group <- ratio_group[order(ratio_group$rank, decreasing = TRUE), ]
  
  #Create ranking table
  ranking_table <- data.frame(
    Geneid = newdata_group$Geneid_base,
    #added the additional values that were calculated to the ranking table, these make the data look like the tables in the prev diss
    rank = newdata_group$rank,
    upstream_region = paste0(newdata_group$Geneid_base, "_f_1"),
    upsteam_counts = newdata_group$smallest_value,
    downstream_region = paste0(newdata_group$Geneid_base, "_f_2"),
    downstream_counts = newdata_group$largest_value
    
    
  )
  
  ranking_tables[[group_name]] <- ranking_table
}
#takes the appropriate column and rounds to 2dp
ranking_OCT4_untreated <- ranking_tables$OCT4_untreated %>% 
  mutate(across(where(is.numeric), ~ round(., 2)))

ranking_OCT4_treated <- ranking_tables$OCT4_treated %>% 
  mutate(across(where(is.numeric), ~ round(., 2)))

ranking_BRG1_untreated <- ranking_tables$BRG1_untreated %>% 
  mutate(across(where(is.numeric), ~ round(., 2)))

ranking_BRG1_treated <- ranking_tables$BRG1_treated %>%
  mutate(across(where(is.numeric), ~ round(., 2)))


#the following adds the gene names by creating a new column. it takes only the first portion of the name and uses two regexes to remove the excess.
#adapted from https://www.youtube.com/watch?v=NjiCr_iiPjA 
ranking_OCT4_untreated$gene <- mapIds(org.Mm.eg.db, keys = gsub("\\.\\d+$", "", gsub("_up_1_.*", "", ranking_OCT4_untreated$Geneid)), keytype = "REFSEQ", column ="SYMBOL")
ranking_BRG1_untreated$gene <- mapIds(org.Mm.eg.db, keys = gsub("\\.\\d+$", "", gsub("_up_1_.*", "", ranking_BRG1_untreated$Geneid)), keytype = "REFSEQ", column ="SYMBOL")
ranking_OCT4_treated$gene <- mapIds(org.Mm.eg.db, keys = gsub("\\.\\d+$", "", gsub("_up_1_.*", "", ranking_OCT4_treated$Geneid)), keytype = "REFSEQ", column ="SYMBOL")
ranking_BRG1_treated$gene <- mapIds(org.Mm.eg.db, keys = gsub("\\.\\d+$", "", gsub("_up_1_.*", "", ranking_BRG1_treated$Geneid)), keytype = "REFSEQ", column ="SYMBOL")
#If want to save to directory 
write.csv(ranking_OCT4_untreated, "1000ranking_OCT4_untreated.csv", row.names = FALSE)
write.csv(ranking_OCT4_treated, "1000ranking_OCT4_treated.csv", row.names = FALSE)
write.csv(ranking_BRG1_untreated, "1000ranking_BRG1_untreated.csv", row.names = FALSE)
write.csv(ranking_BRG1_treated, "1000ranking_BRG1_treated.csv", row.names = FALSE)

#Making an annotation for OCT4, this must be done manually looking at how many samples are assigned to condition
#prepare condition data frame for the analysis 
condition_1 <- c(rep('OCT4_untreated', 15), rep('OCT4_treated', 15))
data_samples_1 <- data.frame(
  row.names = colnames(counts_new)[1:30],  
  condition = condition_1
)

#SRA table needed for deseq2 comparisons between the conditions.
SRA_table <- read.table("SraRunTable.csv", header=TRUE, sep= ',', comment.char="#")
SRA_table <- SRA_table[1:30,]

#checking the data 
all(colnames(counts_new)[1:30] %in% rownames(data_samples_1))
all(colnames(counts_new)[1:30] == rownames(data_samples_1))

#getting only OCT4_untreated and OCT4_treated from counts_new data 
counts_subset <- counts_new[, 1:30]

#do dds1 
dds_1 <- DESeqDataSetFromMatrix(countData = counts_subset,
                                colData = data_samples_1,
                                design = ~ condition)

dds_1 <- collapseReplicates(dds_1, groupby = SRA_table$Run, run = SRA_table$BioSample)



#get the factor scale 
dds_1$condition <- relevel(dds_1$condition, ref = 'OCT4_untreated')

#do deseq2 
dds_1 <- DESeq(dds_1)
res_1 <- results(dds_1, contrast=c("condition","OCT4_treated","OCT4_untreated"))
res1_new <- results(dds_1, contrast=c("condition","OCT4_treated","OCT4_untreated"))

#p-values 
resOrdered <- res_1[order(res_1$pvalue),]

normalized_counts <- counts(dds_1, normalized=TRUE)

# Filter genes with at least 10 counts in total across all samples
filter <- rowSums(normalized_counts) > 10

# Apply the filter to the results
res_filtered <- res_1[filter, ]


#Now can use the result of analysis to get the top25 genes, to investigate whether both _1 and _2 up-regulated or down-regulated or whether both variants are regulated differently. Here the data is sorted: by getting the genes adter DESeq2 analysis, the expression is shown by looking at log2FoldChange parameter. Log2FoldChange is transferred into log and if log is >1 that would mean the gene is up-regulated and if it is 1>, that would mean gene is down-regulated. If log is 1 the gene is marked as unchanged. This way can tell that expression of those genes increased/decreased in OCT4_untreated. The next parameter taht is considered is p adjacent value. If p_adj < 0.05, that suggests that the change in expression is statistically significant. The genes that were choosen: would have statistically significant change, so both pairs must be increase/decrease (if unchanged taht would mean there is no divergent transcription) and their p_adj<0.05.  
#getting the the data frame of deseq2 results for OCT4 
expressed_genes_OCT4 <- as.data.frame(res_filtered)

#getting filtered from lowest to highest p_adjacent value 
expressed_genes_OCT4 <- expressed_genes_OCT4[order(expressed_genes_OCT4$padj),]

#getting geneid from row.nmaes into column so can see both variants _1 and _2
expressed_genes_OCT4 <- rownames_to_column(expressed_genes_OCT4, var = "gene_id")
#sort in increasing order
expressed_genes_OCT4 <- expressed_genes_OCT4[order(expressed_genes_OCT4$padj), ]
#making row names numbers 
row.names(expressed_genes_OCT4) <- expressed_genes_OCT4$row_number

#getting Gneid_base
expressed_genes_OCT4 <- expressed_genes_OCT4 %>%
  # regex was adapted, for the refid UCSC dataset which had a forward and reverse naming suffix with 1 and 2 rather than the previous _2 _3 variants
  mutate(Geneid_base = gsub('_[fr]_[12]$','', gene_id))

expressed_genes_OCT4 <- relocate(expressed_genes_OCT4, Geneid_base, .after=gene_id)

#getting log value and looking whether gene pairs change: increase or decrease in expression. 
expressed_genes_OCT4 <- expressed_genes_OCT4 %>%
  mutate(
    Log = 2^log2FoldChange,
    expression = case_when(
      Log > 1 ~ 'upregulated',
      Log < 1 ~ 'downregulated',
      TRUE ~ 'unchanged'
    )
  ) 

#by sorting with p adjacent < 0.05 can see whether genes actually statistically significantly change or not, if p_adj>0.05 would suggest that change is not significant 
expressed_genes_OCT4 <- expressed_genes_OCT4 %>%
  mutate(
    expression_2 = case_when(
      is.na(padj) ~ 'not significantly changed',
      padj > 0.05 ~ 'not significantly changed',
      TRUE ~ 'significantly changed'
    )
  )





#Now need to get rid of genes that are nit significantly cahnged as would mean that not divergently transcribed (counts too low if NA)
changed_genes_OCT4 <- expressed_genes_OCT4 %>%
  filter(expression_2 == 'significantly changed')

#get both gene_ids based on Geneid_base 
both_geneid <- changed_genes_OCT4 %>%
  group_by(Geneid_base) %>%
  filter(
    (all(c(paste0(Geneid_base, "_f_1"), paste0(Geneid_base, "_f_2")) %in% gene_id)) |
      (all(c(paste0(Geneid_base, "_r_1"), paste0(Geneid_base, "_r_2")) %in% gene_id)) 
  ) %>%
  pull(Geneid_base) %>%
  unique()

#Getting all pairs that are padj < 0.05 and are either upregulated or downregulated 
OCT4_genes <- expressed_genes_OCT4 %>%
  filter(
    Geneid_base %in% both_geneid,
    expression_2 == 'significantly changed'
  )

OCT4_genes_2 <- expressed_genes_OCT4 %>%
  filter(Geneid_base %in% both_geneid, 
         expression_2 == 'significantly changed',
         baseMean > 10) %>%
  group_by(Geneid_base) %>%
  filter(n() == 2) 

OCT4_genes_2 <- OCT4_genes_2[order(OCT4_genes_2$gene_id),]

#sorting by alphabetical order 
OCT4_genes <- OCT4_genes[order(OCT4_genes$gene_id),]

#saving results 
OCT4_genes_2$gene <- mapIds(org.Mm.eg.db, keys = gsub("\\.\\d+$", "", gsub("_up_1_.*", "", OCT4_genes_2$gene_id)), keytype = "REFSEQ", column ="SYMBOL") 
write.csv(OCT4_genes_2, '/home/s2089123/diss/OCTbasemean.csv')

OCT4_genes_up_down_2 <- OCT4_genes_2 %>%
  group_by(Geneid_base) %>%
  filter(any(expression == 'upregulated') & any(expression == 'downregulated')) %>%
  filter(baseMean > 10) %>%
  filter(n() == 2)

#write.csv(OCT4_genes_up_down, '/home/s2089123/diss/OCT4_genes_diff_reg.csv')

OCT4_one_sig_one_not_sig <- expressed_genes_OCT4 %>%
  group_by(Geneid_base) %>%
  filter(any(expression_2 == 'not significantly changed') & any(expression_2 == 'significantly changed')) %>%
  filter(!is.na(padj)) %>%
  filter(baseMean > 10) %>%
  filter(n() == 2)

OCT4_one_sig_one_not_sig <- OCT4_one_sig_one_not_sig[order(OCT4_one_sig_one_not_sig$gene_id), ]
write_csv(OCT4_one_sig_one_not_sig, '/home/s2089123/diss/OCT$_one_sig_one_not_sig.csv')

## Graphs can be plot to vizualize the data quality
plotMA(res_1, ylim=c(-2,2))

## PCA plot:
vsdata <- vst(dds_1, blind = FALSE)
plotPCA(vsdata, intgroup = 'condition')

#The below analysis was done for BRG1 same as for OCT4:
condition_2 <- c(rep('BRG1_untreated', 15), rep('BRG1_treated', 15))
data_samples_2 <- data.frame(
  row.names = colnames(counts_new)[31:60],  
  condition = condition_2
)
  
  
  #checking the data 
  all(colnames(counts_new)[31:60] %in% rownames(data_samples_2))
  all(colnames(counts_new)[31:60] == rownames(data_samples_2))
  #getting only BRG1_untretaed and BRG1_treated from counts_new data 
  counts_subset <- counts_new[, 31:60]
  
  # do dds1 
  dds_2 <- DESeqDataSetFromMatrix(countData = counts_subset,
                                  colData = data_samples_2,
                                  design = ~ condition)
  SRA_table <- read.table("SraRunTable.csv", header=TRUE, sep= ',', comment.char="#")
  SRA_table <- SRA_table[31:60,]
  
  dds_2 <- collapseReplicates(dds_2, groupby = SRA_table$Run, run = SRA_table$BioSample)
  
  #get the factor scale 
  dds_2$condition <- relevel(dds_2$condition, ref = 'BRG1_untreated')
  
  #do deseq2 
  dds_2 <- DESeq(dds_2)
  res_2 <- results(dds_2, contrast=c("condition","BRG1_treated","BRG1_untreated"))
  
  #p-values 
  resOrdered <- res_2[order(res_2$pvalue),]
  
#Now getting the data frame: 
#gettingthe the data frame of deseq2 results for OCT4 
expressed_genes_BRG1 <- as.data.frame(res_2)
#getting filtered frm lowest tp highest p_adjacent value 
  
  expressed_genes_BRG1 <- expressed_genes_BRG1[order(expressed_genes_BRG1$padj),]
  
  expressed_genes_BRG1 <- rownames_to_column(expressed_genes_BRG1, var = "gene_id")
  expressed_genes_BRG1 <- expressed_genes_BRG1[order(expressed_genes_BRG1$padj), ]
  row.names(expressed_genes_BRG1) <- expressed_genes_BRG1$row_number
  
  
  #getting Gene_base from expressed genes to find pairs 
  expressed_genes_BRG1 <- expressed_genes_BRG1 %>%
    mutate(Geneid_base = gsub('_[fr]_[12]$','', gene_id))
  
  expressed_genes_BRG1 <- relocate(expressed_genes_BRG1, Geneid_base, .after=gene_id)
  
  
  #getting log value and looking whether gene pairs change: incarese or decraese in expression 
  expressed_genes_BRG1 <- expressed_genes_BRG1 %>%
    mutate(
      Log = 2^log2FoldChange,
      expression = case_when(
        Log > 1 ~ 'upregulated',
        Log < 1 ~ 'downregulated',
        TRUE ~ 'unchanged'
      )
    ) 
  
  #by sorting with adkustong p value < 0.05 can see whether genes actually change or not, if p_adj>0.05 tehy would suggest that change is nit significant 
  expressed_genes_BRG1 <- expressed_genes_BRG1 %>%
    mutate(
      expression_2 = case_when(
        is.na(padj) ~ 'not significantly changed',
        padj > 0.05 ~ 'not significantly changed',
        TRUE ~ 'significantly changed'
      )
    )
  
  
  
  
  
  #Now need to get rid of genes that are not significantly changed as would mean that not divergently transcribed (counts too low if NA)
  changed_genes_BRG1 <- expressed_genes_BRG1 %>%
    filter(expression_2 == 'significantly changed')
  
  # Identify Geneid_base that have both _2 and _3 variants significantly expressed
  both_ids <- changed_genes_BRG1 %>%
    group_by(Geneid_base) %>%
    filter(
  (all(c(paste0(Geneid_base, "_f_1"), paste0(Geneid_base, "_f_2")) %in% gene_id)) |
  (all(c(paste0(Geneid_base, "_r_1"), paste0(Geneid_base, "_r_2")) %in% gene_id)) 
  ) %>%
    pull(Geneid_base) %>%
    unique()
  
  # Filter the original dataframe to include only genes with valid Geneid_base
  BRG1_genes <- expressed_genes_BRG1 %>%
    filter(
      Geneid_base %in% both_ids,
      expression_2 == 'significantly changed',
    )
  
  BRG1_genes <- BRG1_genes[order(BRG1_genes$gene_id),]
  
  #write.csv(BRG1_genes, '/home/s2614505/Diss/BRG1_genes_padj<0.05.csv')
  
  BRG1_genes_up_down <- BRG1_genes %>%
    group_by(Geneid_base) %>%
    filter(any(expression == 'upregulated') & any(expression == 'downregulated'))
  
  
  BRG1_one_sig_one_not_sig <- expressed_genes_BRG1 %>%
    group_by(Geneid_base) %>%
    filter(any(expression_2 == 'not significantly changed') & any(expression_2 == 'significantly changed')) %>%
    filter(!is.na(padj)) %>%
    filter(baseMean > 10) %>%
    filter(n() == 2) 
  
  
  BRG1_one_sig_one_not_sig <- BRG1_one_sig_one_not_sig[order(BRG1_one_sig_one_not_sig$gene_id), ]
  
  write.csv(BRG1_one_sig_one_not_sig, '/home/s2089123/diss/BRG1_ONE_NOT.csv')
  
  
  
  BRG1_genes_2 <- expressed_genes_BRG1%>%
    filter(Geneid_base %in% both_geneid, 
           expression_2 == 'significantly changed',
           baseMean > 2) %>%
    group_by(Geneid_base) %>%
    filter(n() == 2) 
  
  BRG1_3 <- BRG1_genes %>%
    group_by(Geneid_base) %>%
    filter(baseMean > 10) %>%
    filter(n() == 2) 
  
  BRG1_genes_2 <- BRG1_genes_2[order(BRG1_genes_2$gene_id),]
  BRG1_3 <- BRG1_3[order(BRG1_3$gene_id),]
  
 BRG1_3$gene <- mapIds(org.Mm.eg.db, keys = gsub("\\.\\d+$", "", gsub("_up_1_.*", "", BRG1_3$gene_id)), keytype = "REFSEQ", column ="SYMBOL")
  #sorting by alphabetical order
  #saving results 
  #write.csv(OCT4_genes_2, '/home/s2614505/Diss/BRG1_basemeansort.csv')
  write.csv(BRG1_3, '/home/s2089123/diss/BRG1_basemeansort3.csv')
  
  BRG1_genes_up_down_2 <- BRG1_3 %>%
    group_by(Geneid_base) %>%
    filter(any(expression == 'upregulated') & any(expression == 'downregulated')) %>%
    filter(n() == 2)
  
  #write.csv(OCT4_genes_up_down_2, '/home/s2614505/Diss/OCT4_basemeansortregdif.csv')
  write.csv(BRG1_genes_up_down_2, '/home/s2089123/diss/BRG1_basemeansortregdif.csv')
  
  
  ## Plots for BRG1 data: 
  
  # MA plot: 
  plotMA(res_2, ylim=c(-2,2))
  
  # PCA plot for: 
  vsdata_2 <- vst(dds_2, blind = FALSE)
  plotPCA(vsdata_2, intgroup = 'condition')
  

# Venn Diagramm of gene pairs that differentially expressed

#doing venn diagram for overlaping genes
 genes_venn_OCT4 <- OCT4_genes_2 %>% filter(padj < 0.05) %>% pull(gene_id)
 genes_venn_BRG1 <- BRG1_3 %>% filter(padj < 0.05) %>% pull(gene_id)

genes_one_venn_OCT4 <- OCT4_one_sig_one_not_sig %>% pull(gene_id)
genes_one_venn_BRG1 <- BRG1_one_sig_one_not_sig %>%  pull(gene_id)

overlapping_genes <- intersect(genes_venn_OCT4, genes_venn_BRG1)

overlapping_one_genes <- intersect(genes_one_venn_BRG1, genes_one_venn_OCT4)

 #Generate Venn Diagram with improved styling
 venn.plot <- venn.diagram(
  x = list(OCT4 = genes_venn_OCT4, BRG1 = genes_venn_BRG1),
 category.names = c("OCT4 genes", "BRG1 genes"),
 filename = NULL,
 output = TRUE,
  fill = c("#66C2A5", "#FC8D62"),  # Distinct colors
 alpha = 0.6,  # Adjust transparency
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-30, 30),  # Adjust position of category labels
  cat.dist = c(0.05, 0.05),
  cat.fontfamily = "sans",
  lwd = 2,  # Line width of circles
  lty = 'dashed',  # Line type of circles
  col = c("darkgreen", "darkred"),  # Outline colors
  main = "Divergent pairs where both regions are significantly expressed",
  main.cex = 2,  # Main title size
  main.fontface = "bold",
  main.fontfamily = "sans"
)

 #Draw Venn Diagram
  grid.newpage()
 grid.draw(venn.plot)

 #Display overlapping genes in a prettier way
 overlapping_genes <- intersect(genes_venn_OCT4, genes_venn_BRG1)  # Find overlapping genes

 
 #This section was removed from the figure as they are not gene names, and there is much more significance than in the previous study
 # crowds the figure.
# grid.text(
#   paste(overlapping_genes, collapse = "\n"),
#   x = unit(0.5, "npc"),
#   y = unit(0.2, "npc"),  # Adjust position to avoid overlap with Venn Diagram
#   gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "sans")
#  )


 venn.plot <- venn.diagram(
   x = list(OCT4 = genes_one_venn_OCT4, BRG1 = genes_one_venn_BRG1),
   category.names = c("OCT4 genes", "BRG1 genes"),
  filename = NULL,
   output = TRUE,
   fill = c("#A6D854", "#FFD92F"),
   alpha = 0.6,  # Adjust transparency
   cex = 1.5,
   fontface = "bold",
   fontfamily = "sans",
   cat.cex = 1.5,
   cat.fontface = "bold",
   cat.default.pos = "text",
   cat.pos = c(-30, 30),  # Adjust position of category labels
  cat.dist = c(0.05, 0.05),
   cat.fontfamily = "sans",
   lwd = 2,  # Line width of circles
   lty = 'dashed',  # Line type of circles
   col = c("darkgreen", "darkred"),  # Outline colors
   main = "Divergent pairs where one region is significantly changed",
   main.cex = 2,  # Main title size
   main.fontface = "bold",
   main.fontfamily = "sans"
 )

 #Draw Venn Diagram
 grid.newpage()
 grid.draw(venn.plot)


 # grid.text(
 #   paste(overlapping_one_genes, collapse = "\n"),
 #   x = unit(0.5, "npc"),
 #   y = unit(0.2, "npc"),  # Adjust position to avoid overlap with Venn Diagram
 #   gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "sans")
 # )
 
 #load in the counts data for gro-seq
 groseq_data <- read.table("countsGROSEQ.txt", header=TRUE, sep="\t", comment.char="#")
 rownames(groseq_data) <- groseq_data$Geneid
 grometdata <-groseq_data[c(1:5)]
 untreatedgro_seq <- groseq_data[, 7:9]
 differentiated_gro_seq <- groseq_data[, 10:12]
 
 scale_untreated <- colSums(untreatedgro_seq)
 scale_diff <- colSums(differentiated_gro_seq)
 #normalise each condition
 untreated_groseq_rpm <- sweep(untreatedgro_seq, 2, scale_untreated, FUN='/')*1e6
 diff_groseq_rpm <- sweep(differentiated_gro_seq, 2, scale_diff, FUN='/')*1e6
 
 #combine the two conditions for filtering
 groseq_combined <- cbind(untreated_groseq_rpm, diff_groseq_rpm)
 rownames(groseq_combined) <- rownames(groseq_data)
 
 #this was changed to extract _1 and _2 rather than _2 or _3 like the previous study
 pattern <- '.*_[12]$'
 groseq_filtered <- subset(groseq_combined, grepl(pattern, rownames(groseq_combined)))
 
 #adjusted the groups for the gro-seq data
 gro_groups <- list("groseq_untreated" = 1:3,
                "groseq_differentiated" = 4:6) 
 
 #After groups are difined manually now will go into function that will produce rankings automatically 
 
 #Extracting indices - columns for each group 
 groups <- lapply(gro_groups, function(indices) groseq_filtered[, indices])
 
 #Calculating row sums for each group
 gro_row_sums <- lapply(groups, function(group) rowSums(group, na.rm = FALSE))
 
 #Getting everything into a data frame
 gro_RPM_normal <- data.frame(Geneid_RPM = rownames(groseq_filtered))
 for (i in 1:length(gro_row_sums)) {
   gro_RPM_normal <- cbind(gro_RPM_normal, gro_row_sums[[i]])
   
 }
 #Reanmimg the column names according to groups 
 colnames(gro_RPM_normal)[-1] <- names(gro_groups)
 
 #Checking the data 
 head(gro_RPM_normal)
 
 #Now need to get Geneid_base
 gro_RPM_normal <- gro_RPM_normal %>% 
   mutate(Geneid_base = gsub('_[fr]_[12]$','', Geneid_RPM)) 
 
 # Craeting kist for rankings of each group 
 gro_ranking_tables <- list()
 
 #Counts will be set to not be less then 10, change manually if need more/less
 number <- 20 
 
 for (group_name in names(gro_groups)) {
   ratio_group <- gro_RPM_normal %>%
     group_by(Geneid_base) %>%
     summarize(
       smallest_value = min(.data[[group_name]], na.rm = TRUE),
       largest_value = max(.data[[group_name]], na.rm = TRUE),
       largest_ratio = ifelse(smallest_value > 0, largest_value / smallest_value, NA),
       rank = ifelse(smallest_value >= number, smallest_value * (1 / largest_ratio), NA)
     ) %>% 
     filter(smallest_value >= number)
   
   
   #Order by rank in descending order
   newdata_group <- ratio_group[order(ratio_group$rank, decreasing = TRUE), ]
   
   #Create ranking table
   gro_ranking_table <- data.frame(
     Geneid = newdata_group$Geneid_base,
     #added the additional values that were calculated to the ranking table, these make the data look like the tables in the prev diss
     rank = newdata_group$rank,
     upstream_region = paste0(newdata_group$Geneid_base, "_f_1"),
     upsteam_counts = newdata_group$smallest_value,
     downstream_region = paste0(newdata_group$Geneid_base, "_f_2"),
     downstream_counts = newdata_group$largest_value
     
     
   )
   
   gro_ranking_tables[[group_name]] <- gro_ranking_table
 }
 ranking_groseq_untreated <- gro_ranking_tables$groseq_untreated %>% 
   mutate(across(where(is.numeric), ~ round(., 2)))
 
 ranking_groseq_differentiated <- gro_ranking_tables$groseq_differentiated %>% 
   mutate(across(where(is.numeric), ~ round(., 2)))

 
#map the refseq id to gene name so both can be displayed
 ranking_groseq_untreated$gene <- mapIds(org.Mm.eg.db, keys = gsub("\\.\\d+$", "", gsub("_up_1_.*", "", ranking_groseq_untreated$Geneid)), keytype = "REFSEQ", column ="SYMBOL")
 ranking_groseq_differentiated$gene <- mapIds(org.Mm.eg.db, keys = gsub("\\.\\d+$", "", gsub("_up_1_.*", "", ranking_groseq_differentiated$Geneid)), keytype = "REFSEQ", column ="SYMBOL")
 write.csv(ranking_groseq_untreated, "1000ranking_groseq_untreated.csv", row.names = FALSE)
 write.csv(ranking_groseq_differentiated, "1000ranking_groseq_differentiated.csv", row.names = FALSE)

 rnaseq_expressed <- ranking_OCT4_untreated$Geneid
 groseq_untreated_expressed <- ranking_groseq_untreated$Geneid
 both <- intersect(rnaseq_expressed, groseq_untreated_expressed)
 length(both)
 #checks the genes that are not in the RNA seq untreated from the gro-seq and makes a list. aka the 263 genes that are only in gro-seq
 groseq_specific <- ranking_groseq_untreated[!ranking_groseq_untreated$Geneid %in% ranking_OCT4_untreated$Geneid & !ranking_groseq_untreated$Geneid %in% ranking_tt_untreated$Geneid,]
 write.csv(groseq_specific, "groseq_only.csv", row.names=FALSE)
 
 #comparsion between rnaseq and groseq venn diagram
 venn.plot <- venn.diagram(
   x = list(RNASeq = rnaseq_expressed, GROSeq = groseq_untreated_expressed),
   category.names = c("RNA-Seq Untreated", "GRO-Seq Untreated"),
   filename = NULL,
   output = TRUE,
   fill = c("#D81B60", "#1E88E5"),  # Distinct colors
   alpha = 0.6,  # Adjust transparency
   cex = 1.5,
   fontface = "bold",
   fontfamily = "sans",
   cat.cex = 1.5,
   cat.fontface = "bold",
   cat.default.pos = "text",
   cat.pos = c(-30, 30),  # Adjust position of category labels
   cat.dist = c(0.05, 0.05),
   cat.fontfamily = "sans",
   lwd = 2,  # Line width of circles
   lty = 'dashed',  # Line type of circles
   col = c("#D81B60", "#1E88E5"),  # Outline colors
   main = "Divergent pairs where both regions are significantly expressed",
   main.cex = 2,  # Main title size
   main.fontface = "bold",
   main.fontfamily = "sans",
   scaled = FALSE
 )
 
 grid.newpage()
 grid.draw(venn.plot)
 
 tt_seq_data <- read.table("countsttseq.txt", header=TRUE, sep="\t", comment.char="#")
 rownames(tt_seq_data) <- tt_seq_data$Geneid
 tt_untreated <-tt_seq_data[, 7:8]
 
 scale_tt <- colSums(tt_untreated)
 untreated_tt_rmp <- sweep(tt_untreated, 2, scale_tt, FUN='/')*1e6

 #this was changed to extract _1 and _2 rather than _2 or _3 like the previous study
 pattern <- '.*_[12]$'
 tt_filtered <- subset(untreated_tt_rmp, grepl(pattern, rownames(untreated_tt_rmp)))
 
 #adjusted the groups for the gro-seq data
 tt_groups <- list("ttseq_untreated" = 1:2) 
 
 #After groups are difined manually now will go into function that will produce rankings automatically 
 
 #Extracting indices - columns for each group 
 groups <- lapply(tt_groups, function(indices) tt_filtered[, indices])
 
 #Calculating row sums for each group
 tt_row_sums <- lapply(groups, function(group) rowSums(group, na.rm = FALSE))
 
 #Getting everything into a data frame
 tt_RPM_normal <- data.frame(Geneid_RPM = rownames(tt_filtered))
 for (i in 1:length(tt_row_sums)) {
   tt_RPM_normal <- cbind(tt_RPM_normal, tt_row_sums[[i]])
   
 }
 #Reanmimg the column names according to groups 
 colnames(tt_RPM_normal)[-1] <- names(tt_groups)
 
 #Checking the data 
 head(tt_RPM_normal)
 
 #Now need to get Geneid_base
 tt_RPM_normal <- tt_RPM_normal %>% 
   mutate(Geneid_base = gsub('_[fr]_[12]$','', Geneid_RPM)) 
 
 # Craeting kist for rankings of each group 
 tt_ranking_tables <- list()
 
 #Counts will be set to not be less then 10, change manually if need more/less
 number <- 20
 
 for (group_name in names(tt_groups)) {
   ratio_group <- tt_RPM_normal %>%
     group_by(Geneid_base) %>%
     summarize(
       smallest_value = min(.data[[group_name]], na.rm = TRUE),
       largest_value = max(.data[[group_name]], na.rm = TRUE),
       largest_ratio = ifelse(smallest_value > 0, largest_value / smallest_value, NA),
       rank = ifelse(smallest_value >= number, smallest_value * (1 / largest_ratio), NA)
     ) %>% 
     filter(smallest_value >= number)
   
   
   #Order by rank in descending order
   newdata_group <- ratio_group[order(ratio_group$rank, decreasing = TRUE), ]
   
   #Create ranking table
   tt_ranking_table <- data.frame(
     Geneid = newdata_group$Geneid_base,
     #added the additional values that were calculated to the ranking table, these make the data look like the tables in the prev diss
     rank = newdata_group$rank,
     upstream_region = paste0(newdata_group$Geneid_base, "_f_1"),
     upsteam_counts = newdata_group$smallest_value,
     downstream_region = paste0(newdata_group$Geneid_base, "_f_2"),
     downstream_counts = newdata_group$largest_value
     
     
   )
   
   tt_ranking_tables[[group_name]] <- tt_ranking_table
 }
 ranking_tt_untreated <- tt_ranking_tables$ttseq_untreated %>% 
   mutate(across(where(is.numeric), ~ round(., 2)))
 ranking_tt_untreated$gene <- mapIds(org.Mm.eg.db, keys = gsub("\\.\\d+$", "", gsub("_up_1_.*", "", ranking_tt_untreated$Geneid)), keytype = "REFSEQ", column ="SYMBOL")
 write.csv(ranking_tt_untreated, "1000ranking_tt_untreated.csv", row.names = FALSE)
 
 tt_expressed <- ranking_tt_untreated$Geneid
 
 #min expression of 20 rpm ranking between the two conditions
 venn.plot <- venn.diagram(
   x = list(RNASeq = rnaseq_expressed, GROSeq = groseq_untreated_expressed, TTSeq = tt_expressed),
   category.names = c("RNA-Seq Untreated", "GRO-Seq Untreated", "TT-Seq Untreated"),
   filename = NULL,
   output = TRUE,
   fill = c("#D81B60", "#1E88E5", "#FFC107"),
   alpha = 0.6,  # Adjust transparency
   cex = 1.5,
   fontface = "bold",
   fontfamily = "sans",
   cat.cex = 1.5,
   cat.fontface = "bold",
   cat.default.pos = "text",
   cat.pos = c(-30, 30, 180),  # Adjust position of category labels
   cat.dist = c(0.05, 0.05, 0.05),
   cat.fontfamily = "sans",
   lwd = 2,  # Line width of circles
   lty = 'dashed',  # Line type of circles
   col = c("#D81B60", "#1E88E5", "#FFC107"),  # Outline colors
   main = "Divergent expression in RNA-seq Vs Gro-seq Vs TT-Seq",
   main.cex = 2,  # Main title size
   main.fontface = "bold",
   main.fontfamily = "sans",
   scaled = FALSE
 )
 grid.newpage()
 grid.draw(venn.plot)
   #scaled needed to remove the weirdness of the venn diagram when the numbers are diff orders or magnitude
 
 #genes that are in both gro-seq and tt-seq but not rna-seq
 both_tt_gro <- ranking_groseq_untreated[ranking_groseq_untreated$Geneid %in% tt_expressed &!ranking_groseq_untreated$Geneid %in% rnaseq_expressed, ]
 write.csv(both_tt_gro, "genes_in_both_gro_and_tt.csv", row.names = FALSE)
 
 #genes that are in all three
 all_three <- ranking_OCT4_untreated[ranking_OCT4_untreated$Geneid %in% ranking_groseq_untreated$Geneid & ranking_OCT4_untreated$Geneid %in% ranking_tt_untreated$Geneid,]
 write.csv(all_three, "genes_in_all_three.csv", row.names = FALSE)

#take all of the IDs from all three conditions  
ids <- unique(c(ranking_OCT4_untreated$Geneid, ranking_groseq_untreated$Geneid, ranking_tt_untreated$Geneid))
all_three_table <- data.frame (Geneid = ids, gene = mapIds(org.Mm.eg.db, keys = gsub("\\.\\d+$", "", gsub("_up_1_.*", "", ids)), keytype = "REFSEQ", column ="SYMBOL"),
                              RNA = ranking_OCT4_untreated$rank[match(ids, ranking_OCT4_untreated$Geneid)],
                              GRO = ranking_groseq_untreated$rank[match(ids, ranking_groseq_untreated$Geneid)],
                              TT = ranking_tt_untreated$rank[match(ids, ranking_tt_untreated$Geneid)]
                              ) 
write.csv(all_three_table, "genes_in_three_conditions.csv", row.names = FALSE)
