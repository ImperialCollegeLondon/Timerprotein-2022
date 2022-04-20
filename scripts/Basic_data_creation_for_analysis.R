## Basic data import
library("DESeq2")
library("tximport")
library('dplyr')
library("ashr")
library("biomaRt")
library("magrittr")



## Clone names TBJ 3.60 and TBX4B; often abbreviated to TBJ and TBX


## Create a sample info datatable
metadata <- read.table(file="/home/hkiik/HTLV/metadata/WTCHG_617523_626991_sample_metadata.txt", header=TRUE, stringsAsFactors = FALSE)
read_name <- unique(substr(metadata$Oxford_Genomics_code, 14, 19))
sample_names <- paste(metadata$Sample_ID[1:16], c(unique(metadata$State)),sep = "_")
clones <- data.frame(clones = factor(c(metadata$Sample_ID[1:16])))
condition <- data.frame(condition=factor(c("DN", "Blue", "DP", "Red")))
replicate_nr <- data.frame(replicate=factor(c("1", "1", "1", "1", "2", "2", "2", "2")))
clone_condition <- as.factor(paste(clones[,1], condition[,1], sep="_"))
sample_info <- data.frame(condition, clones, replicate_nr, clone_condition, row.names = sample_names)
sample_info$condition <- relevel(sample_info$condition, ref = "DN")

## Import files
rsem_files <- paste("/home/hkiik/HTLV/data/RSEM_counts_2022/Quant_",read_name,"_Default.genes.results", sep = "")
names(rsem_files) <- sample_names
txi.rsem <- tximport(rsem_files, type="rsem", txIn=FALSE, txOut = FALSE, geneIdCol = "gene_id", importer = read.delim)

## Add a pseudo-length of 1.
txi.rsem$length[txi.rsem$length < 1] <- 1

## Make a SummarizedExperiment 
data_set_rsem <- DESeqDataSetFromTximport(txi.rsem, colData = sample_info, design = ~clone_condition)

## Filter some unnecessary entries for quicker biomart retrieval for full dataset.
keep <- rowSums(counts(data_set_rsem) >= 3) >= 2
data_set_rsem <- data_set_rsem[keep,]

# Separate clones before DESeq2 for clone-specific parameter estimation.
TBX_dds_rsem <- data_set_rsem[,1:8]
TBJ_dds_rsem <- data_set_rsem[,9:16]

TBX_dds_rsem$clones <- droplevels(TBX_dds_rsem$clones)
TBX_dds_rsem$clone_condition <- droplevels(TBX_dds_rsem$clone_condition)

TBJ_dds_rsem$clones <- droplevels(TBJ_dds_rsem$clones)
TBJ_dds_rsem$clone_condition <- droplevels(TBJ_dds_rsem$clone_condition)

## Filter out 0-counts for each clone
rsem_keep_tbx <- rowSums(counts(TBX_dds_rsem) >= 3) >= 2
rsem_keep_tbj <- rowSums(counts(TBJ_dds_rsem) >= 3) >= 2

TBX_dds_rsem <- TBX_dds_rsem[rsem_keep_tbx,]
TBJ_dds_rsem <- TBJ_dds_rsem[rsem_keep_tbj,]


## Annotation retrieval using BioMart
ensembl100 <- useEnsembl(biomart="ensembl", version=100)
ensembl100 <- useDataset("hsapiens_gene_ensembl", mart=ensembl100)

biomart_results_external_rsem <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol","description", "chromosome_name", "strand"), filters= "ensembl_gene_id",
                                       values = rownames(data_set_rsem),
                                       mart=ensembl100, 
                                       useCache = FALSE)

## Create dataframe for ensembl metadata from biomart
gene_id_rsem <- as.data.frame(rownames(data_set_rsem))
colnames(gene_id_rsem) <- "ensembl_gene_id"

## Obtain unique/non-duplicated ensembl codes
unique_symbols_rsem <- biomart_results_external_rsem[!duplicated(biomart_results_external_rsem$ensembl_gene_id),]

# Choose a representative of the duplicates instead of removing both
biomart_results_external_rsem[duplicated(biomart_results_external_rsem$ensembl_gene_id),]

biomart_results_external_rsem[biomart_results_external_rsem$ensembl_gene_id == "ENSG00000230417",]
unique_symbols_rsem[unique_symbols_rsem$ensembl_gene_id == "ENSG00000230417",] <- biomart_results_external_rsem[biomart_results_external_rsem$ensembl_gene_id == "ENSG00000230417",][1,]

biomart_results_external_rsem[biomart_results_external_rsem$ensembl_gene_id == "ENSG00000255374",]
unique_symbols_rsem[unique_symbols_rsem$ensembl_gene_id == "ENSG00000255374",] <- biomart_results_external_rsem[biomart_results_external_rsem$ensembl_gene_id == "ENSG00000255374",][2,]

biomart_results_external_rsem[biomart_results_external_rsem$ensembl_gene_id == "ENSG00000276085",]
unique_symbols_rsem[unique_symbols_rsem$ensembl_gene_id == "ENSG00000276085",] <- biomart_results_external_rsem[biomart_results_external_rsem$ensembl_gene_id == "ENSG00000276085",][1,]



## Merge into a look-up table
gene_names_rsem <- merge(gene_id_rsem, unique_symbols_rsem, by="ensembl_gene_id", all.x=TRUE, sort=TRUE)

## Replace NA values back to their original ensembl IDs for a more complete set.
gene_names_rsem$external_gene_name[is.na(gene_names_rsem$external_gene_name)] <- as.character(gene_names_rsem$ensembl_gene_id[is.na(gene_names_rsem$external_gene_name)])

# Fill in HTLV-1 information
gene_names_rsem$hgnc_symbol[is.na(gene_names_rsem$hgnc_symbol)] <- as.character(gene_names_rsem$external_gene_name[is.na(gene_names_rsem$hgnc_symbol)])
gene_names_rsem$hgnc_symbol[gene_names_rsem$hgnc_symbol == ""]<- as.character(gene_names_rsem$external_gene_name[gene_names_rsem$hgnc_symbol == ""])

gene_names_rsem$ensembl_gene_id <- as.character(gene_names_rsem$ensembl_gene_id)
gene_names_rsem$description[gene_names_rsem$ensembl_gene_id == "HTLVplus"] <- c("HTLV-1 plus strand")
gene_names_rsem$strand[gene_names_rsem$ensembl_gene_id == "HTLVplus"] <- c("1")
gene_names_rsem$chromosome_name[gene_names_rsem$ensembl_gene_id == "HTLVplus"] <- c("TBX:22, TJB:4")
gene_names_rsem$hgnc_symbol[gene_names_rsem$ensembl_gene_id == "HTLVplus"] <- c("HTLV-1 plus")


gene_names_rsem$description[gene_names_rsem$ensembl_gene_id == "HTLVminus"] <- c("HTLV-1 negative strand")
gene_names_rsem$chromosome_name[gene_names_rsem$ensembl_gene_id == "HTLVminus"] <- c("TBX:22, TJB:4")
gene_names_rsem$strand[gene_names_rsem$ensembl_gene_id == "HTLVminus"] <- c("-1")
gene_names_rsem$hgnc_symbol[gene_names_rsem$ensembl_gene_id == "HTLVminus"] <- c("HTLV-1 minus")

gene_names_rsem$description[gene_names_rsem$ensembl_gene_id == "Timer_protein"] <- c("Sam's timer protein")
gene_names_rsem$chromosome_name[gene_names_rsem$ensembl_gene_id == "Timer_protein"] <- c("Timerprotein")
gene_names_rsem$strand[gene_names_rsem$ensembl_gene_id == "Timer_protein"] <- c("Timerprotein")
gene_names_rsem$hgnc_symbol[gene_names_rsem$ensembl_gene_id == "Timer_protein"] <- c("Timer")



#######################################################################################
## Pairwise DE comparisons ############################################################
#######################################################################################
alpha_value <- 0.01
# Design
design(TBX_dds_rsem) <- formula(~0+condition)
TBX_dds_rsem_deseq_design <- DESeq(TBX_dds_rsem)

# Blue vs DN
rsem_condition_Blue_vs_DN_tbx <- results(TBX_dds_rsem_deseq_design,
                                         contrast=c("condition", "Blue", "DN"),
                                         alpha = alpha_value)

joined_Blue_vs_DN_table_tbx <-as.data.frame(rsem_condition_Blue_vs_DN_tbx) %>% na.omit
# write.csv(joined_Blue_vs_DN_table_tbx, file="./DE_lists/TBX4B_Blue_vs_DN_results.csv", quote = F, row.names = T, col.names = T )
joined_Blue_vs_DN_table_tbx <-as.data.frame(rsem_condition_Blue_vs_DN_tbx)

joined_Blue_vs_DN_table_tbx$SYMBOL <- gene_names_rsem[match(rownames(rsem_condition_Blue_vs_DN_tbx), gene_names_rsem$ensembl_gene_id),]$hgnc_symbol
joined_Blue_vs_DN_table_tbx$chr <- gene_names_rsem[match(rownames(rsem_condition_Blue_vs_DN_tbx), gene_names_rsem$ensembl_gene_id),]$chromosome_name

design(TBJ_dds_rsem) <- formula(~0+condition)
TBJ_dds_rsem_deseq_design <- DESeq(TBJ_dds_rsem)
rsem_condition_Blue_vs_DN_tbj <- results(TBJ_dds_rsem_deseq_design,
                                         contrast=c("condition", "Blue", "DN"),
                                         alpha = alpha_value)

joined_Blue_vs_DN_table_tbj <-as.data.frame(rsem_condition_Blue_vs_DN_tbj) %>% na.omit
# write.csv(joined_Blue_vs_DN_table_tbj, file="./DE_lists/360_Blue_vs_DN_results.csv", quote = F,  row.names = T, col.names = T )
joined_Blue_vs_DN_table_tbj <-as.data.frame(rsem_condition_Blue_vs_DN_tbj)

joined_Blue_vs_DN_table_tbj$SYMBOL <- gene_names_rsem[match(rownames(joined_Blue_vs_DN_table_tbj), gene_names_rsem$ensembl_gene_id),]$hgnc_symbol
joined_Blue_vs_DN_table_tbj$chr <- gene_names_rsem[match(rownames(joined_Blue_vs_DN_table_tbj), gene_names_rsem$ensembl_gene_id),]$chromosome_name

## DP vs DN
rsem_condition_DP_vs_DN_tbx <- results(TBX_dds_rsem_deseq_design,
                                       contrast=c("condition", "DP", "DN"),
                                       alpha = alpha_value)


joined_DP_vs_DN_table_tbx <-as.data.frame(rsem_condition_DP_vs_DN_tbx) %>% na.omit
# write.csv(joined_DP_vs_DN_table_tbx, file="./DE_lists/TBX4B_DP_vs_DN_results.csv", quote = F,  row.names = T, col.names = T )
joined_DP_vs_DN_table_tbx <-as.data.frame(rsem_condition_DP_vs_DN_tbx) 

joined_DP_vs_DN_table_tbx$SYMBOL <- gene_names_rsem[match(rownames(joined_DP_vs_DN_table_tbx), gene_names_rsem$ensembl_gene_id),]$hgnc_symbol
joined_DP_vs_DN_table_tbx$chr <- gene_names_rsem[match(rownames(joined_DP_vs_DN_table_tbx), gene_names_rsem$ensembl_gene_id),]$chromosome_name

rsem_condition_DP_vs_DN_tbj <- results(TBJ_dds_rsem_deseq_design,
                                       contrast=c("condition", "DP", "DN"),
                                       alpha = alpha_value)


joined_DP_vs_DN_table_tbj <-as.data.frame(rsem_condition_DP_vs_DN_tbj) %>% na.omit
# write.csv(joined_DP_vs_DN_table_tbj, file="./DE_lists/360_DP_vs_DN_results.csv", quote = F,  row.names = T, col.names = T )
joined_DP_vs_DN_table_tbj <-as.data.frame(rsem_condition_DP_vs_DN_tbj)

joined_DP_vs_DN_table_tbj$SYMBOL <- gene_names_rsem[match(rownames(joined_DP_vs_DN_table_tbj), gene_names_rsem$ensembl_gene_id),]$hgnc_symbol
joined_DP_vs_DN_table_tbj$chr <- gene_names_rsem[match(rownames(joined_DP_vs_DN_table_tbj), gene_names_rsem$ensembl_gene_id),]$chromosome_name

# Red vs DN
rsem_condition_Red_vs_DN_tbx <- results(TBX_dds_rsem_deseq_design,
                                        contrast=c("condition", "Red", "DN"),
                                        alpha = alpha_value)

joined_Red_vs_DN_table_tbx <-as.data.frame(rsem_condition_Red_vs_DN_tbx) %>% na.omit
# write.csv(joined_Red_vs_DN_table_tbx, file="./DE_lists/TBX4B_Red_vs_DN_results.csv", quote = F,  row.names = T, col.names = T )
joined_Red_vs_DN_table_tbx <-as.data.frame(rsem_condition_Red_vs_DN_tbx)

joined_Red_vs_DN_table_tbx$SYMBOL <- gene_names_rsem[match(rownames(joined_Red_vs_DN_table_tbx), gene_names_rsem$ensembl_gene_id),]$hgnc_symbol
joined_Red_vs_DN_table_tbx$chr <- gene_names_rsem[match(rownames(joined_Red_vs_DN_table_tbx), gene_names_rsem$ensembl_gene_id),]$chromosome_name

rsem_condition_Red_vs_DN_tbj <- results(TBJ_dds_rsem_deseq_design,
                                        contrast=c("condition", "Red", "DN"),
                                        alpha = alpha_value)

joined_Red_vs_DN_table_tbj <-as.data.frame(rsem_condition_Red_vs_DN_tbj) %>% na.omit
# write.csv(joined_Red_vs_DN_table_tbj, file="./DE_lists/360_Red_vs_DN_results.csv", quote = F,  row.names = T, col.names = T )
joined_Red_vs_DN_table_tbj <-as.data.frame(rsem_condition_Red_vs_DN_tbj)

joined_Red_vs_DN_table_tbj$SYMBOL <- gene_names_rsem[match(rownames(joined_Red_vs_DN_table_tbj), gene_names_rsem$ensembl_gene_id),]$hgnc_symbol
joined_Red_vs_DN_table_tbj$chr <- gene_names_rsem[match(rownames(joined_Red_vs_DN_table_tbj), gene_names_rsem$ensembl_gene_id),]$chromosome_name

#######################################################################################
## Sequential pairwise comparisons ####################################################
#######################################################################################

## Red vs Blue
rsem_condition_Red_vs_Blue_tbx <- results(TBX_dds_rsem_deseq_design,
                                          contrast=c("condition", "Red", "Blue"),
                                          alpha = alpha_value)

joined_Red_vs_Blue_table_tbx <-as.data.frame(rsem_condition_Red_vs_Blue_tbx) %>% na.omit
# write.csv(joined_Red_vs_Blue_table_tbx, file="./DE_lists/TBX4B_Red_vs_Blue_results.csv", quote = F,  row.names = T, col.names = T )
joined_Red_vs_Blue_table_tbx <-as.data.frame(rsem_condition_Red_vs_Blue_tbx)

joined_Red_vs_Blue_table_tbx$SYMBOL <- gene_names_rsem[match(rownames(rsem_condition_Red_vs_Blue_tbx), gene_names_rsem$ensembl_gene_id),]$hgnc_symbol
joined_Red_vs_Blue_table_tbx$chr <- gene_names_rsem[match(rownames(rsem_condition_Red_vs_Blue_tbx), gene_names_rsem$ensembl_gene_id),]$chromosome_name

rsem_condition_Red_vs_Blue_tbj <- results(TBJ_dds_rsem_deseq_design,
                                          contrast=c("condition", "Red", "Blue"),
                                          alpha = alpha_value)

joined_Red_vs_Blue_table_tbj <-as.data.frame(rsem_condition_Red_vs_Blue_tbj) %>% na.omit
# write.csv(joined_Red_vs_Blue_table_tbj, file="./DE_lists/360_Red_vs_Blue_results.csv", quote = F,  row.names = T, col.names = T )
joined_Red_vs_Blue_table_tbj <-as.data.frame(rsem_condition_Red_vs_Blue_tbj)

joined_Red_vs_Blue_table_tbj$SYMBOL <- gene_names_rsem[match(rownames(joined_Red_vs_Blue_table_tbj), gene_names_rsem$ensembl_gene_id),]$hgnc_symbol
joined_Red_vs_Blue_table_tbj$chr <- gene_names_rsem[match(rownames(joined_Red_vs_Blue_table_tbj), gene_names_rsem$ensembl_gene_id),]$chromosome_name

## DP vs Blue
rsem_condition_DP_vs_Blue_tbx <- results(TBX_dds_rsem_deseq_design,
                                         contrast=c("condition", "DP", "Blue"),
                                         alpha = alpha_value)

joined_DP_vs_Blue_table_tbx <-as.data.frame(rsem_condition_DP_vs_Blue_tbx) %>% na.omit
# write.csv(joined_DP_vs_Blue_table_tbx, file="./DE_lists/TBX4B_DP_vs_Blue_results.csv", quote = F,  row.names = T, col.names = T )
joined_DP_vs_Blue_table_tbx <-as.data.frame(rsem_condition_DP_vs_Blue_tbx)

joined_DP_vs_Blue_table_tbx$SYMBOL <- gene_names_rsem[match(rownames(joined_DP_vs_Blue_table_tbx), gene_names_rsem$ensembl_gene_id),]$hgnc_symbol
joined_DP_vs_Blue_table_tbx$chr <- gene_names_rsem[match(rownames(joined_DP_vs_Blue_table_tbx), gene_names_rsem$ensembl_gene_id),]$chromosome_name

rsem_condition_DP_vs_Blue_tbj <- results(TBJ_dds_rsem_deseq_design,
                                         contrast=c("condition", "DP", "Blue"),
                                         alpha = alpha_value)

joined_DP_vs_Blue_table_tbj <-as.data.frame(rsem_condition_DP_vs_Blue_tbj) %>% na.omit
# write.csv(joined_DP_vs_Blue_table_tbj, file="./DE_lists/360_DP_vs_Blue_results.csv", quote = F,  row.names = T, col.names = T )
joined_DP_vs_Blue_table_tbj <-as.data.frame(rsem_condition_DP_vs_Blue_tbj)

joined_DP_vs_Blue_table_tbj$SYMBOL <- gene_names_rsem[match(rownames(joined_DP_vs_Blue_table_tbj), gene_names_rsem$ensembl_gene_id),]$hgnc_symbol
joined_DP_vs_Blue_table_tbj$chr <- gene_names_rsem[match(rownames(joined_DP_vs_Blue_table_tbj), gene_names_rsem$ensembl_gene_id),]$chromosome_name

# Red vs DP
rsem_condition_Red_vs_DP_tbx <- results(TBX_dds_rsem_deseq_design,
                                        contrast=c("condition", "Red", "DP"),
                                        alpha = alpha_value)

joined_Red_vs_DP_table_tbx <-as.data.frame(rsem_condition_Red_vs_DP_tbx) %>% na.omit
# write.csv(joined_Red_vs_DP_table_tbx, file="./DE_lists/TBX4B_Red_vs_DP_results.csv", quote = F,  row.names = T, col.names = T )
joined_Red_vs_DP_table_tbx <-as.data.frame(rsem_condition_Red_vs_DP_tbx)

joined_Red_vs_DP_table_tbx$SYMBOL <- gene_names_rsem[match(rownames(joined_Red_vs_DP_table_tbx), gene_names_rsem$ensembl_gene_id),]$hgnc_symbol
joined_Red_vs_DP_table_tbx$chr <- gene_names_rsem[match(rownames(joined_Red_vs_DP_table_tbx), gene_names_rsem$ensembl_gene_id),]$chromosome_name

rsem_condition_Red_vs_DP_tbj <- results(TBJ_dds_rsem_deseq_design,
                                        contrast=c("condition", "Red", "DP"),
                                        alpha = alpha_value)

joined_Red_vs_DP_table_tbj <-as.data.frame(rsem_condition_Red_vs_DP_tbj) %>% na.omit
# write.csv(joined_Red_vs_DP_table_tbj, file="./DE_lists/360_Red_vs_DP_results.csv", quote = F,  row.names = T, col.names = T )
joined_Red_vs_DP_table_tbj <-as.data.frame(rsem_condition_Red_vs_DP_tbj)

joined_Red_vs_DP_table_tbj$SYMBOL <- gene_names_rsem[match(rownames(joined_Red_vs_DP_table_tbj), gene_names_rsem$ensembl_gene_id),]$hgnc_symbol
joined_Red_vs_DP_table_tbj$chr <- gene_names_rsem[match(rownames(joined_Red_vs_DP_table_tbj), gene_names_rsem$ensembl_gene_id),]$chromosome_name



#######################################################################################
## LRT ################################################################################
#######################################################################################
## LRT design
design(TBJ_dds_rsem) <- ~ condition
dds_TBJ <- DESeq(TBJ_dds_rsem, test="LRT", reduced = ~1)
design(TBX_dds_rsem) <- ~ condition
dds_TBX <- DESeq(TBX_dds_rsem, test="LRT", reduced = ~1)

############ Clone TBJ 3.60
dds_res_TBJ <- results(dds_TBJ, name="condition_Blue_vs_DN", alpha = 0.01)

## Filter significant results and add HGNC symbols.
dds_res_sig_data_TBJ <- as.data.frame(dds_res_TBJ)%>%
  subset(padj < 0.01)%>%
  arrange(desc(baseMean))%>%
  arrange(padj)

dds_res_sig_data_TBJ$symbol<-gene_names_rsem$hgnc_symbol[match(rownames(dds_res_sig_data_TBJ), gene_names_rsem$ensembl_gene_id)] 

# Deal with duplicates by taking the one with the highest basemean.
duplicated_tbj<-dds_res_sig_data_TBJ$symbol[duplicated(dds_res_sig_data_TBJ$symbol)]

for(i in 1:length(duplicated_tbj)){
  gene<-duplicated_tbj[i]
  duplicated_res<-dds_res_sig_data_TBJ[dds_res_sig_data_TBJ$symbol %in% gene,]
  gene_id<-row.names(duplicated_res)
  max_value <-which.max(dds_res_sig_data_TBJ[dds_res_sig_data_TBJ$symbol %in% gene,]$baseMean)
  dds_res_sig_data_TBJ[dds_res_sig_data_TBJ$symbol %in% gene,]<-duplicated_res[gene_id[max_value],]
}

# Unique entries, ordered based on baseMean and adjusted p-value.
dds_res_sig_TBJ <- dds_res_sig_data_TBJ%>%
  unique()%>%
  arrange(desc(baseMean))%>%
  arrange(padj)%>%
  row.names()

# Add gene names
dds_res_table_TBJ <- dds_res_TBJ
dds_res_table_TBJ$symbol <- gene_names_rsem$hgnc_symbol[match(row.names(dds_res_table_TBJ), gene_names_rsem$ensembl_gene_id)]

# Save data
data_to_save_TBJ <- dds_res_TBJ %>% na.omit
# write.csv(data_to_save_TBJ, file="./DE_lists/360_LRT_results.csv", quote = F, row.names = T, col.names = T )

############ Clone TBX4B

dds_res_TBX <- results(dds_TBX,name="condition_Blue_vs_DN", alpha = 0.01)
dds_res_table_TBX <- dds_res_TBX


## Filter significant results and add HGNC symbols
dds_res_sig_data_TBX <- as.data.frame(dds_res_TBX)%>%
  subset(padj < 0.01)%>%
  arrange(desc(baseMean))%>%
  arrange(padj)

dds_res_table_TBX$symbol <- gene_names_rsem$hgnc_symbol[match(row.names(dds_res_table_TBX), gene_names_rsem$ensembl_gene_id)]
duplicated_tbx<-dds_res_sig_data_TBX$symbol[duplicated(dds_res_sig_data_TBX$symbol)]
# No duplicates

# Order by baseMean and adjusted p-value.
dds_res_TBX_sig <- as.data.frame(dds_res_TBX)%>%
  subset(padj < 0.01)%>%
  arrange(desc(baseMean))%>%
  arrange(padj)%>%
  row.names()
data_to_save_TBX <- dds_res_TBX %>% na.omit
# write.csv(data_to_save_TBX, file="./DE_lists/TBX4B_LRT_results.csv", quote = F,  row.names = T, col.names = T )

save.image(file = "./Basic_data_for_analysis.RData")

# > sessionInfo()
# R version 4.1.2 (2021-11-01)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.6 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] remotes_2.4.2               extrafont_0.17              ggtext_0.1.1                PCAtools_2.6.0             
# [5] ggrepel_0.9.1               tibble_3.1.6                tidyr_1.2.0                 gprofiler2_0.2.1           
# [9] patchwork_1.1.1             ggplot2_3.3.5               ashr_2.2-54                 dplyr_1.0.8                
# [13] magrittr_2.0.3              tximport_1.22.0             BiocManager_1.30.16         DESeq2_1.34.0              
# [17] SummarizedExperiment_1.24.0 Biobase_2.54.0              MatrixGenerics_1.6.0        matrixStats_0.61.0         
# [21] GenomicRanges_1.46.1        GenomeInfoDb_1.30.1         IRanges_2.28.0              S4Vectors_0.32.4           
# [25] BiocGenerics_0.40.0         biomaRt_2.50.3             
# 
# loaded via a namespace (and not attached):
#   [1] colorspace_2.0-3          ellipsis_0.3.2            markdown_1.1              XVector_0.34.0           
# [5] gridtext_0.1.4            farver_2.1.0              bit64_4.0.5               AnnotationDbi_1.56.2     
# [9] fansi_1.0.3               xml2_1.3.3                splines_4.1.2             sparseMatrixStats_1.6.0  
# [13] cachem_1.0.6              geneplotter_1.72.0        knitr_1.38                jsonlite_1.8.0           
# [17] Rttf2pt1_1.3.10           annotate_1.72.0           dbplyr_2.1.1              png_0.1-7                
# [21] compiler_4.1.2            httr_1.4.2                dqrng_0.3.0               lazyeval_0.2.2           
# [25] assertthat_0.2.1          Matrix_1.4-1              fastmap_1.1.0             cli_3.2.0                
# [29] BiocSingular_1.10.0       htmltools_0.5.2           prettyunits_1.1.1         tools_4.1.2              
# [33] rsvd_1.0.5                gtable_0.3.0              glue_1.6.2                GenomeInfoDbData_1.2.7   
# [37] reshape2_1.4.4            rappdirs_0.3.3            Rcpp_1.0.8.3              vctrs_0.4.0              
# [41] Biostrings_2.62.0         extrafontdb_1.0           DelayedMatrixStats_1.16.0 xfun_0.30                
# [45] stringr_1.4.0             beachmat_2.10.0           lifecycle_1.0.1           irlba_2.3.5              
# [49] XML_3.99-0.9              zlibbioc_1.40.0           scales_1.1.1              hms_1.1.1                
# [53] parallel_4.1.2            RColorBrewer_1.1-3        yaml_2.3.5                curl_4.3.2               
# [57] memoise_2.0.1             yulab.utils_0.0.4         stringi_1.7.6             RSQLite_2.2.12           
# [61] SQUAREM_2021.1            genefilter_1.76.0         ScaledMatrix_1.2.0        filelock_1.0.2           
# [65] BiocParallel_1.28.3       truncnorm_1.0-8           rlang_1.0.2               pkgconfig_2.0.3          
# [69] bitops_1.0-7              evaluate_0.15             lattice_0.20-45           invgamma_1.1             
# [73] purrr_0.3.4               labeling_0.4.2            htmlwidgets_1.5.4         cowplot_1.1.1            
# [77] bit_4.0.4                 tidyselect_1.1.2          plyr_1.8.7                R6_2.5.1                 
# [81] generics_0.1.2            DelayedArray_0.20.0       DBI_1.1.2                 pillar_1.7.0             
# [85] withr_2.5.0               survival_3.3-1            KEGGREST_1.34.0           RCurl_1.98-1.6           
# [89] mixsqp_0.3-43             crayon_1.5.1              utf8_1.2.2                BiocFileCache_2.2.1      
# [93] plotly_4.10.0             rmarkdown_2.13            progress_1.2.2            locfit_1.5-9.5           
# [97] grid_4.1.2                data.table_1.14.2         blob_1.2.3                digest_0.6.29            
# [101] xtable_1.8-4              gridGraphics_0.5-1        munsell_0.5.0             viridisLite_0.4.0        
# [105] ggplotify_0.1.0   
