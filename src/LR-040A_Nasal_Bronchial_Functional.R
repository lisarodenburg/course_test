
# Functional analysis of RNA seq data
# Created by Lisa Rodenburg
# November 2021
#
#Based on RNA seq analysis tutorial:
#https://github.com/hbctraining/DGE_workshop/tree/master/lessons
#

##### Install packages #####

library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(AnnotationHub)
library(ensembldb)
library(tidyverse)
library(annotables)
library(stats)
library(gprofiler2)
library(treemap)
library(enrichplot)
library(KEGGREST)

##### code below only necessary when need to install packages: #####

# First specify the packages of interest
  packages = c("org.Hs.eg.db", "DOSE", "pathview", "clusterProfiler", "AnnotationHub", "ensembldb", "tidyverse", "here")
  
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
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("GenomicFeatures")
  library("GenomicFeatures")
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("gprofiler2")
  library("gprofiler2")
  
  install.packages("treemap")
  
##### ClusterProfiler for over-representation analysis #####
    
## Explore the grch37 table loaded by the annotables library
grch37 <- grch37
  View(grch37)
  
## Set thresholds
  padj.cutoff <- 0.01
  lfc.cutoff <- 0.58

## Open results table
filename_results_table = paste0("results/Results_table_Nasal_Bronchial_shrinken_corrected.csv")

## Remove chromosome number from gene name  
Results = read.csv2 (file= filename_results_table, header=T, sep =";", row.names=1,check.names=F)
Results_names <- Results %>%
  as.data.frame() %>%
  mutate(chromosome = gsub("__","", str_extract(rownames(Results), "__chr[a-zA-Z0-9]{1,2}"))) 
rownames(Results_names)[grep("RPL21__chr10", rownames(Results_names), fixed = TRUE)] <- "RPL21_10" # this gene is from 2 chromosomes, so adapt name to avoid duplicate gene names
rownames(Results_names)[grep("RPL21__chr13", rownames(Results_names), fixed = TRUE)] <- "RPL21_13"
rownames(Results_names) <- gsub("__chr[a-zA-Z0-9]{1,2}", "", rownames(Results_names))

## Return the IDs for the gene symbols in the DE results
idx <- grch37$symbol %in% rownames(Results_names)
ids <- grch37[idx, ]

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates <- which(duplicated(ids$symbol) == FALSE)
ids <- ids[non_duplicates, ] 

## Merge the IDs with the results 
Results_names <- Results_names %>% rownames_to_column(var="gene")
res_ids <- inner_join(Results_names, ids, by=c("gene"="symbol"))  

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
all_genes <- as.character(res_ids$ensgene)
all_genes_entrez <- as.character(res_ids$entrez)

## Extract significant results
sig <- dplyr::filter(res_ids, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
sig_genes <- as.character(sig$ensgene)
sig_genes_entrez <- as.character(sig$entrez)

## Run GO over-representation analysis 
ego_bp <- enrichGO(gene = sig_genes, 
                universe = all_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

ego_all <- enrichGO(gene = sig_genes, 
                universe = all_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "ALL", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

## Run KEGG pathway over-representation analysis
eKEGG <- enrichKEGG(gene = sig_genes_entrez, 
                universe = all_genes_entrez,
                keyType = "kegg",
                organism = "hsa", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05)

## Output results from analysis to a table and save
cluster_summary_bp <- data.frame(ego_bp)

cluster_summary_all <- data.frame(ego_all)

cluster_summary_kegg <- data.frame(eKEGG) 

#t <- cluster_summary_kegg$geneID[1]

#t2 <- gsub("/","''", t)

#t2
#split(t2, ceiling(seq_along(t2)/1))

#getSYMBOL(c('3815', '3816', '2341'), data='org.Hs.eg')

write.csv2(cluster_summary_bp, "results/clusterProfiler_bp.csv")

write.csv2(cluster_summary_all, "results/clusterProfiler_all.csv")

write.csv2(cluster_summary_kegg, "results/clusterProfiler_kegg.csv")

## GO plots 
 # Dotplot: Top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not p-adjusted value.
pdf(file = "results/dotplot_over-representation_GO-bp.pdf", width = 8, height = 14)
dotplot(ego_bp, showCategory=50)
dev.off()

pdf(file = "results/dotplot_over-representation_GO-all.pdf", width = 8, height = 14)
dotplot(ego_all, showCategory=50)
dev.off()

pdf(file = "results/dotplot_over-representation_kegg.pdf", width = 8, height = 8)
dotplot(eKEGG, showCategory=50)
dev.off()

# Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
ego2_all <-pairwise_termsim(ego_all)
pdf(file = "results/emapplot_over-representation_GO-all.pdf", width = 24, height = 18)
emapplot(ego2_all, showCategory = 50)
dev.off()

ego2_bp <-pairwise_termsim(ego_bp)
pdf(file = "results/emapplot_over-representation_GO-bp.pdf", width = 24, height = 18)
emapplot(ego2_bp, showCategory = 50)
dev.off()

ekegg2 <-pairwise_termsim(eKEGG)
pdf(file = "results/emapplot_over-representation_kegg.pdf", width = 12, height = 8)
emapplot(ekegg2, showCategory = 50)
dev.off()

# Category netplot from top 5 GO terms

# To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
sig_foldchanges <- sig$log2FoldChange
names(sig_foldchanges) <- sig$gene

# Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego_all, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=sig_foldchanges, 
         vertex.label.font=6,
         cex_label_gene = 0.6)

# If some of the high fold changes are getting drowned out due to a large range, you could set a maximum fold change value
#sig_foldchanges <- ifelse(sig_foldchanges > 2, 2, sig_foldchanges)
#sig_foldchanges <- ifelse(sig_foldchanges < -2, -2, sig_foldchanges)
#cnetplot(ego_all, 
 #        categorySize="pvalue", 
  #       showCategory = 5, 
   #      foldChange=sig_foldchanges, 
    #     vertex.label.font=6)

# Category netplot from GO terms of interest (choose row numbers from cluster_summary Excel table)

# Subsetting the ego results without overwriting original `ego` variable
View(ego_all@result)
ego3 <- ego_all
ego3@result <- ego2@result[c(8,9,14,27),] #row number of GO terms of interest

# Plotting terms of interest
cnetplot(ego3, 
         categorySize="pvalue", 
         foldChange=sig_foldchanges, 
         showCategory = 5, 
         vertex.label.font=6,
         cex_label_gene = 0.6)

##### Functional class scoring with GSEA #####

## For GO terms ##

  # we want the log2 fold change 
     original_gene_list <- res_ids$log2FoldChange
  
  # name the vector
    names(original_gene_list) <- res_ids$ensgene
  
  # omit any NA values 
    gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
    gene_list = sort(gene_list, decreasing = TRUE)
  
  #gse <- gseGO(geneList=gene_list, 
                ont ="BP", 
                keyType = "ENSEMBL", 
              # nPerm = 10000, 
              # minGSSize = 20, 
               #maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE, 
                OrgDb = org.Hs.eg.db, 
                pAdjustMethod = "BH")
  
  gse_results <- gse@result
  
  write.csv2(gse_results, "results/gsea_GO.csv", quote=F)
  
  require(DOSE)
  #dotplot(gse, showCategory=50)
  dotplot(gse, showCategory=25, split=".sign") + facet_grid(.~.sign)
  
  

## For KEGG pathways ##

    ## Remove any NA values
    res_entrez <- filter(res_ids, entrez != "NA")
        
    ## Remove any Entrez duplicates
    res_entrez <- res_entrez[which(duplicated(res_entrez$entrez) == F), ]   
    
    ## Extract the foldchanges
    foldchanges <- res_entrez$log2FoldChange
    
    ## Name each fold change with the corresponding Entrez ID
    names(foldchanges) <- res_entrez$entrez
    
    ## Sort fold changes in decreasing order
    foldchanges <- sort(foldchanges, decreasing = TRUE)
    
    ## GSEA using gene sets from KEGG pathways
    # No significant genes???
    gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
                        organism = "hsa", # supported organisms listed below
                        #keyType = "kegg",
                        #nPerm = 1000, # default number permutations
                        minGSSize = 2, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                        pvalueCutoff = 0.5, # padj cutoff value
                        verbose = T)
    head(gseaKEGG)
    
    ## Extract the GSEA results
    gseaKEGG_results <- gseaKEGG@result
    
    ## Write GSEA results to file
    View(gseaKEGG_results)
    
    write.csv2(gseaKEGG_results, "results/gsea_kegg.csv", quote=F)
    
    ## Explore the GSEA plot of enrichment of the pathway-associated genes in the ranked list:
      
      ## Plot the GSEA plot for a single enriched pathway, `hsa03040`
    kegg = 'hsa04080'  
    gseaplot(gseaKEGG, geneSetID = kegg, title = kegg)
    kegg = 'hsa00980'
    kegg = 'hsa04750'
    kegg = 'hsa04020'
    kegg = 'hsa05146'
    kegg = 'hsa04550'
      
      ## Output images for any single significant KEGG pathway
      detach(package:tidyverse, unload=TRUE) # first unload dplyr to avoid conflicts
      pathview(gene.data = foldchanges,
               pathway.id = kegg,
               species = "hsa",
               limit = list(gene = 2, # value gives the max/min limit for foldchanges
                            cpd = 1))
     # NOTE: Printing out Pathview images for all significant pathways can be easily performed as follows:
        
        ## Output images for all significant KEGG pathways
        get_kegg_plots <- function(x) {
          pathview(gene.data = foldchanges, pathway.id = gseaKEGG_results$ID[x], species = "hsa", 
                   limit = list(gene = 2, cpd = 1))
        }
      
      purrr::map(1:length(gseaKEGG_results$ID), get_kegg_plots)
      
# GSEA using gene sets associated with BP Gene Ontology terms
      # This one does work! 
      # To do: save for all significant ones. Include name in title?
     # gseaGO <- gseGO(geneList = foldchanges, 
             #         OrgDb = org.Hs.eg.db, 
              #        ont = 'BP', 
               #       #nPerm = 1000, 
                #      minGSSize = 20, 
                 #     pvalueCutoff = 0.05,
                  #    verbose = FALSE) 
      
      #gseaGO_results <- gseaGO@result
      
      gseaplot(gseaGO, geneSetID = 'GO:0005216') # ion channel activity
      gseaplot(gseaGO, geneSetID = 'GO:0005254') # chloride channel activity
      
      #require(DOSE)
      #dotplot(gseaGO, showCategory=50)
      #dotplot(gseaGO, showCategory=25, split=".sign") + facet_grid(.~.sign)
      
      
##### Pathway topology ##### 
      
      # Set-up
      
      source("http://bioconductor.org/biocLite.R") 
      biocLite("SPIA")
      library(SPIA)
      
      if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      
      BiocManager::install("SPIA")
      library("SPIA")
      
      ## Significant genes is a vector of fold changes where the names are ENTREZ gene IDs. The background set is a vector of all the genes represented on the platform.
      
      background_entrez <- res_entrez$entrez
      sig_res_entrez <- res_entrez[which(res_entrez$padj < padj.cutoff), ]
      sig_entrez <- sig_res_entrez$log2FoldChange
      names(sig_entrez) <- sig_res_entrez$entrez
      head(sig_entrez)
      
      spia_result <- spia(de=sig_entrez, all=background_entrez, organism="hsa")
      
      head(spia_result, n=20)
      
      plotP(spia_result, threshold=0.05)
      
      ## Look at pathway 03013 and view kegglink
      subset(spia_result, ID == "04020")
      subset(spia_result, ID == "05146")
      subset(spia_result, ID == "04970")
      subset(spia_result, ID == "04960")
      
