################################
# Gene-level differential expression analysis using DESeq2
# Created by Lisa Rodenburg
# November 2021
#
# Based on RNA seq analysis tutorial:
# https://github.com/hbctraining/DGE_workshop/tree/master/lessons
# More information: vignette("DESeq2")
################################

##### Install and load packages #####

library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)
library(ggplot2)
library(ggrepel)
library(here)

##### code below only necessary when need to install packages: #####

Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)

# To install/load package DEseq2
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("DESeq2")
  library("DEseq2")

# To install/load package apeglm
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("apeglm")
  library("apeglm")

# To install/load package DEGreport
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("DEGreport")
  library("DEGreport")

# To install/load package org.Hs.eg.db  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db")
  library("org.Hs.eg.db")
  
  # To install/load package org.Hs.eg.db  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("DOSE")
  library("DOSE")
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(c('pathview', 'clusterProfiler', 'AnnotationHub', 'ensembldb', 'annotables'))
  library('pathview')
  
  library('AnnotationHub')
  library('ensembldb')
 
  
  install.packages("devtools")
  devtools::install_github("stephenturner/annotables")
  library('annotables')
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("clusterProfiler")
  library('clusterProfiler')
  
# First specify the packages of interest
  packages = c('matrixStats', 'pheatmap', 'SummarizedExperiment', 'tidyverse', 'ggplot2', 'ggrepel', 'RColorBrewer')

# Now load or install&load all
  package.check <- lapply(
    packages,
    FUN = function(x) {
      if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE)
      }
    }
  )

  
##### Set parameters #####

# Set working direction, experiment name
  getwd()
  experiment = "_Nasal_Bronchial"
  
# Set thresholds
  padj.cutoff <- 0.01
  lfc.cutoff <- 0.58 #remember that we are working with log2 fold changes so this translates to an actual fold change of 1.5 which is pretty reasonable.
  
# Set design (variable on which groups need to be defined)
  # "The design formula should have all of the factors in your metadata that account for major sources of variation in your data. The last factor entered in the formula should be the condition of interest."
  design <- ~ Condition
  design_corrected <- ~ Donor + Condition #last factor needs to be the condition of interest

# Set contrast (first condition of interest, last baseline condition)  
  contrast <- c("Condition", "Nasal", "Bronchial") # betekent neus tov bronchiaal, dus bronchiaal is baseline
  
##### Output files #####
  
# Define names of output files
  filename_counts = paste0(("data/Counts"),experiment,(".csv"))
  filename_coldata = paste0(("meta/Coldata"),experiment,(".csv"))
  filename_normalized = paste0(("results/Normalized_counts"),experiment,(".csv"))
  filename_PCA = paste0(("results/PCA"),experiment,(".pdf"))
  filename_hierclust = paste0(("results/Hierclust"),experiment,(".pdf"))
  filename_results_table = paste0(("results/Results_table"),experiment)
  filename_DEG = paste0(("results/DEG"),experiment)
  filename_MA = paste0(("results/MA"),experiment,(".pdf"))
  filename_heatmap = paste0(("results/heatmap"),experiment,(".pdf"))
  filename_volcano = paste0(("results/volcano"),experiment,(".pdf"))
  filename_GEcluster = paste0(("results/GEcluster"),experiment,(".csv"))
  
##### Load input #####
# Note: row names of coldata should be in same order as column names of counts
  
# Load the counts
  counts = read.csv (file= here(filename_counts), header=T, sep =";", row.names=1,check.names=F)
  head(counts)
  dim(counts)
  class(counts) # should be data.frame
  
# Load the metadata
  coldata = read.csv (file= filename_coldata, header=T, sep =";", row.names=1, check.names=F)
  head(coldata)
  dim(coldata)
  class(coldata) # should be data.frame
  
# Check that sample names match in both files
  all(colnames(counts) %in% rownames(coldata))
  all(colnames(counts) == rownames(coldata))
  
##### DEseq #####
  
# Make a DESeq object, without and with donor/batch correction:
  dds <- DESeqDataSetFromMatrix(countData = counts , colData = coldata, design = design)
  dds_corrected <- DESeqDataSetFromMatrix(countData = counts , colData = coldata, design = design_corrected)
  
# Normalize counts (no difference for corrected dds)
  # Normalization is to account for differences in library sizes and RNA composition between samples
  # DEseq2 median ratio (= normalization outcome) can eb used for comparisons between samples, NOT for within sample comparisons!
  dds2 <- estimateSizeFactors(dds)
  sizeFactors(dds2) #to see normalization factors for each sample
  normalized_counts <- counts(dds2, normalized=TRUE)
  normalized_counts <- normalized_counts %>%
    as.data.frame() 
  normalized_counts_names <- normalized_counts %>%
    mutate(chromosome = gsub("__","", str_extract(rownames(normalized_counts), "__chr[a-zA-Z0-9]{1,2}"))) 
  rownames(normalized_counts_names)[grep("RPL21__chr10", rownames(normalized_counts_names), fixed = TRUE)] <- "RPL21_10" # this gene is from 2 chromosomes, so adapt name to avoid duplicate gene names
  rownames(normalized_counts_names)[grep("RPL21__chr13", rownames(normalized_counts_names), fixed = TRUE)] <- "RPL21_13"
  rownames(normalized_counts_names) <- gsub("__chr[a-zA-Z0-9]{1,2}", "", rownames(normalized_counts_names))
  write.csv2(normalized_counts_names, file= filename_normalized, quote=F)
  
##### QC plots (PCA plot and hierarchical clustering) #####  
  
  # log2-transformation of the normalized counts improves the distances/clustering for visualization. DESeq2 uses a regularized log transform (rlog) of the normalized counts for sample-level QC as it moderates the variance across the mean, improving the clustering
  # PCA plot can be used to find sources of variation for PC1 and PC2, which can be accounted for later in the model (except for variable of interest), such as donor or sex.
  
# PCA plot and save as pdf = quality control at sample level
    #by default looks at top 500 genes. ntop can be used to change number of genes
  rld <- rlog(dds2, blind = TRUE)
  pdf(file = filename_PCA, onefile = T)
  plotPCA(rld, intgroup = "Condition")
  plotPCA(rld, intgroup = "Donor")
  dev.off()
  
# Hierarchical clustering, to identify patterns and outliers
    # Samples below 0.80 may indicate an outlier or sample contamination
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat) 
  pdf(file = filename_hierclust)
  pheatmap(rld_cor)
  dev.off()
  
##### Results tables #####
  
# DEseq2 for statistics bewteeen groups
  # NB: note version number of DEseq2 package for publication
  dds2<-DESeq(dds)
  dds2_corrected<-DESeq(dds_corrected)
   
# Create results table for normal and corrected data, without and with lfcshrinking
  # To Generate more accurate log2 foldchange estimates: shrinkage of the LFC estimates toward zero
        # If the shrinkage estimator apeglm is used in published research, please cite:
        # Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics. 10.1093/bioinformatics/bty895
   
  Results_table <- results(dds2, contrast=contrast, alpha = padj.cutoff) #alpha = FDR, kan worden aangepast
   
  Results_table_shrinken <- lfcShrink(dds2, coef=2, res=Results_table)
   
  Results_table_corrected <- results(dds2_corrected, contrast=contrast, alpha = padj.cutoff) #alpha = FDR, kan worden aangepast
  
  Results_table_shrinken_corrected <- lfcShrink(dds2_corrected, coef=6, res=Results_table_corrected)
  
  # separate chromosome from gene name, better for visualization
  Results_table_shrinken_corrected_names <- Results_table_shrinken_corrected %>%
    as.data.frame() %>%
    mutate(chromosome = gsub("__","", str_extract(rownames(Results_table_shrinken_corrected), "__chr[a-zA-Z0-9]{1,2}"))) 
  rownames(Results_table_shrinken_corrected_names)[grep("RPL21__chr10", rownames(Results_table_shrinken_corrected_names), fixed = TRUE)] <- "RPL21_10" # this gene is from 2 chromosomes, so adapt name to avoid duplicate gene names
  rownames(Results_table_shrinken_corrected_names)[grep("RPL21__chr13", rownames(Results_table_shrinken_corrected_names), fixed = TRUE)] <- "RPL21_13"
  rownames(Results_table_shrinken_corrected_names) <- gsub("__chr[a-zA-Z0-9]{1,2}", "", rownames(Results_table_shrinken_corrected_names))
  
  # order for highest significance, then log2foldchange
  Results_table <- Results_table[order(Results_table$padj, Results_table$log2FoldChange), ] 
  Results_table_shrinken <- Results_table_shrinken[order(Results_table_shrinken$padj, Results_table_shrinken$log2FoldChange), ]
  Results_table_corrected <- Results_table_corrected[order(Results_table_corrected$padj, Results_table_corrected$log2FoldChange), ]
  Results_table_shrinken_corrected <- Results_table_shrinken_corrected[order(Results_table_shrinken_corrected$padj, Results_table_shrinken_corrected$log2FoldChange), ]
  
  # save
  write.csv2(Results_table, file = paste0(filename_results_table,'.csv'), quote = F)
  write.csv2(Results_table_shrinken, file = paste0(filename_results_table,'_shrinken.csv'), quote = F)
  write.csv2(Results_table_corrected, file = paste0(filename_results_table,'_corrected.csv'), quote = F)
  write.csv2(Results_table_shrinken_corrected, file = paste0(filename_results_table,'_shrinken_corrected.csv'), quote = F)
   
  # summary of significantly up- and downregulated Genes
  summary(Results_table)
  summary(Results_table_shrinken)
  summary(Results_table_corrected)
  summary(Results_table_shrinken_corrected)
   
# MA plot and save as pdf
   par(mfrow = c(2, 2))
   
   pdf(file = filename_MA, onefile = T)
   plotMA(Results_table, ylim=c(-10,10), main = "non-corrected")
   plotMA(Results_table_shrinken, ylim=c(-10,10), main = "non-corrected_shrunken")
   plotMA(Results_table_corrected, ylim=c(-10,10), main = "donor_corrected")
   plotMA(Results_table_shrinken_corrected, ylim=c(-10,10), main = "donor_corrected_shrunken")
   dev.off() 
   
##### Filter significant Genes #####
   
# Convert to tibble
   Res_shrinken_corrected_tb <- Results_table_shrinken_corrected_names %>%
     data.frame() %>%
     rownames_to_column(var="Gene") %>% 
     as_tibble()
   
# Filter significant Genes
   Sig_shrinken_corrected <- Res_shrinken_corrected_tb %>%
     filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) %>%
     arrange(padj, log2FoldChange)
   
   write.csv2(Sig_shrinken_corrected,paste0(filename_DEG,".csv"),quote = F)
   
# Separate positive versus negatively DEG's
   Sig_shrinken_corrected_positive <- Sig_shrinken_corrected %>%
     filter(log2FoldChange > 0) 
   
   Sig_shrinken_corrected_negative <- Sig_shrinken_corrected %>%
     filter(log2FoldChange < 0)
   
   write.csv2(Sig_shrinken_corrected_positive,paste0(filename_DEG,"_positive.csv"),quote = F)
   write.csv2(Sig_shrinken_corrected_negative,paste0(filename_DEG,"_negative.csv"),quote = F)
   
   
##### Visualizing results, data preparation #####
   
# Data preparation
   
# Create tibbles including row name
   coldata_tb <- coldata %>% 
     rownames_to_column(var="Sample") %>% 
     as_tibble()
   
   normalized_counts_tb <- normalized_counts_names %>% 
     data.frame(check.names=F) %>%
     rownames_to_column(var="gene") %>% 
     as_tibble()
   
# Order results by padj values and extract top 20 sign Genes, also separate for top 20 downregulated and top 20 upregulated genes
   top20_sig_genes_shr_corr <- Res_shrinken_corrected_tb %>% 
     arrange(padj) %>% 	#Arrange rows by padj values
     pull(Gene) %>% 		#Extract character vector of ordered Genes
     head(n=20) 		#Extract the first 20 Genes
   
   top20_sig_genes_shr_corr_positive <- Sig_shrinken_corrected_positive %>% 
     arrange(padj) %>% 	#Arrange rows by padj values
     pull(Gene) %>% 		#Extract character vector of ordered Genes
     head(n=20)
   
   top20_sig_genes_shr_corr_negative <- Sig_shrinken_corrected_negative %>% 
     arrange(padj) %>% 	#Arrange rows by padj values
     pull(Gene) %>% 		#Extract character vector of ordered Genes
     head(n=20)
   
   top20_sig_norm_shr_corr <- normalized_counts_tb %>%
     filter(gene %in% top20_sig_genes_shr_corr)
   
   top20_sig_norm_shr_corr_positive <- normalized_counts_tb %>%
     filter(gene %in% top20_sig_genes_shr_corr_positive)
   
   top20_sig_norm_shr_corr_negative <- normalized_counts_tb %>%
     filter(gene %in% top20_sig_genes_shr_corr_negative)
   
# Gathering the columns to have normalized counts to a single column
   gathered_top20_sig <- top20_sig_norm_shr_corr %>%
     gather(colnames(top20_sig_norm_shr_corr)[2:(nrow(coldata)+1)], key = "Sample", value = "normalized_counts") 
   
   gathered_top20_sig_positive <- top20_sig_norm_shr_corr_positive %>%
     gather(colnames(top20_sig_norm_shr_corr_positive)[2:(nrow(coldata)+1)], key = "Sample", value = "normalized_counts") 
  
    gathered_top20_sig_negative <- top20_sig_norm_shr_corr_negative %>%
     gather(colnames(top20_sig_norm_shr_corr_negative)[2:(nrow(coldata)+1)], key = "Sample", value = "normalized_counts") 
   
# Add sample names from column data
   gathered_top20_sig <- inner_join(coldata_tb, gathered_top20_sig)
   
   gathered_top20_sig_positive <- inner_join(coldata_tb, gathered_top20_sig_positive)
   
   gathered_top20_sig_negative <- inner_join(coldata_tb, gathered_top20_sig_negative)
   
##### Visualize result, making graphs #####

# Plot expression for single Genes
   plotCounts(dds2_corrected, gene="FOXJ1__chr17", intgroup="Condition") #Add returnData = T is you want to use it for ggplot with formula below
   plotCounts(dds2_corrected, gene="MUC5AC__chr11", intgroup="Condition")
  
# To create nicer gg plots with donor names
   individual_genes <- c("MUC5AC__chr11", "FOXJ1__chr17", "SCGB1A1__chr11", "TP63__chr3", "ITGA6__chr2", "KRT5__chr12", "GPX2__chr14", "SFTPB__chr2", "SPDEF__chr6", "CFTR__chr7", "PAX6__chr11", "SIX3__chr2", "SIX3-AS1__chr2", "FOXA2__chr20", "NKX2-1__chr14", "SOX2__chr3", "SOX2-OT__chr3", "FOXG1__chr14", "FOXG1-AS1__chr14", "OTX2__chr14", "OTX2-AS1__chr14")
   
   pdf(file = "results/individual_genes.pdf")

   for (i in 1:length(individual_genes)) {
     print(i)
     g <- plotCounts(dds2_corrected, gene=individual_genes[i], intgroup="Condition", returnData = T) 
     coldata_tb <- coldata %>% rownames_to_column(var="sample")
     g <- g %>% rownames_to_column(var="sample")
     g2 <- left_join (g, coldata_tb, by = "sample")
     
       g3 <- ggplot(g2, aes(x = Condition.x, y = count, color = Donor)) + 
       geom_point(position=position_jitter(w = 0.1,h = 0)) +
       #geom_text_repel(aes(label = Donor)) + 
       theme_bw() +
       ggtitle(individual_genes[i]) +
       theme(plot.title = element_text(hjust = 0.5)) +
       theme(axis.title.x = element_blank(), legend.title = element_blank())
     print(g3)
   }
   
   dev.off()
   
# Plot top20 significant Genes using ggplot2
   
   ggplot(gathered_top20_sig) +
     geom_point(aes(x = gene, y = normalized_counts, color = Condition)) +
     scale_y_log10() +
     xlab("Genes") +
     ylab("log10 Normalized Counts") +
     ggtitle("Top 20 Significant DE Genes") +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
     theme(plot.title = element_text(hjust = 0.5)) + 
     theme(legend.title = element_blank(), axis.title.x = element_blank())
   
   ggplot(gathered_top20_sig_positive) +
     geom_point(aes(x = gene, y = normalized_counts, color = Condition)) +
     scale_y_log10() +
     xlab("Genes") +
     ylab("log10 Normalized Counts") +
     ggtitle("Top 20 Significantly upregulated Genes") +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
     theme(plot.title = element_text(hjust = 0.5)) + 
     theme(legend.title = element_blank(), axis.title.x = element_blank())
   
   ggplot(gathered_top20_sig_negative) +
     geom_point(aes(x = gene, y = normalized_counts, color = Condition)) +
     scale_y_log10() +
     xlab("Genes") +
     ylab("log10 Normalized Counts") +
     ggtitle("Top 20 Significantly downregulated Genes") +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
     theme(plot.title = element_text(hjust = 0.5)) + 
     theme(legend.title = element_blank(), axis.title.x = element_blank())
   
# Heatmap of significantly different expressed genes
   
   # Filter all significant different normalized Genes
   norm_sig <- normalized_counts_names[,1:nrow(coldata)] %>% 
     rownames_to_column(var="Gene") %>%
     filter(Gene %in% Sig_shrinken_corrected$Gene) %>% 
     data.frame(check.names=F) %>%
     column_to_rownames(var = "Gene")
   
   annotation <- coldata_tb %>% 
     select(Sample, Condition) %>% 
     data.frame(row.names = "Sample")
   
   # Set a color palette
   heat_colors <- brewer.pal(6, "YlOrRd") 
   
   pdf(file = filename_heatmap)
   pheatmap(norm_sig, 
            #color = heat_colors, 
            cluster_rows = T, 
            show_rownames = F,
            annotation = annotation, 
            border_color = NA, 
            fontsize = 10, 
            scale = "row", 
            fontsize_row = 10, 
            height = 20,
            main = "Heatmap of DEG")
   dev.off()
   
# Volcano plot
   
   # Add column to define differentially expressed genes
   Res_shrinken_corrected_tb <- Res_shrinken_corrected_tb %>% 
   mutate(threshold = padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
   
   #First select which genes to label (top 10 with lowest asjusted p-values)
   Res_shrinken_corrected_tb <- Res_shrinken_corrected_tb %>% 
     arrange(padj) %>% 
     mutate(genelabels = "")
   Res_shrinken_corrected_tb$genelabels[1:20] <- Res_shrinken_corrected_tb$Gene[1:20]
  
   pdf(file = filename_volcano)
   
   ggplot(Res_shrinken_corrected_tb, aes(x = log2FoldChange, y = -log10(padj))) +
     geom_point(aes(colour = threshold)) +
     geom_text_repel(aes(label = genelabels)) +
     ggtitle(experiment) +
     xlab("log2 fold change") + 
     ylab("-log10 adjusted p-value") +
     #scale_y_continuous(limits = c(0,50)) +
     theme(legend.position = "none",
           plot.title = element_text(size = rel(1.5), hjust = 0.5),
           axis.title = element_text(size = rel(1.25)))  
   
   dev.off()
   
##### Session info #####   

   sessionInfo()
   


   

 
  
  
  #Heatmap of interesting Genes
  
  #Selection of DEG's
  is.double(resdds2_NB$padj)
  head(resdds2_NB)
  resdds2_NB_lowp <-   filter(resdds2_NB, padj != NA)
  head(resdds2_NB)
  
  x <- resdds2_NB$padj < 0.01
 resdds2_NB[,x]
 length(x)
  
  dds2_NB_plow0.01 <- resdds2_NB[,resdds2_NB(na.omit(resdds2_NB$padj < 0.01))]
  
  #volcano plot
  par(mfrow=c(1,1))
  with(resdds2_NB_donor_corrected, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-10,10), ylim=c(0,10)))
  
      # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
      with(subset(resdds2_NB_donor_corrected, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
      with(subset(resdds2_NB_donor_corrected, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
      text(resdds2_NB_donor_corrected$padj, resdds2_NB_donor_corrected$log2FoldChange, rownames(resdds2_NB_donor_corrected))
      

  

