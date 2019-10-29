library(DESeq2)
library(tidyverse)
library(pheatmap)
library(viridis)
library(clusterProfiler)

#load functions
source("bin/functions.R")

#protein coding genes data
                      counts <- read.csv("data/4veh-vs-3stz-mavsko-protein-coding-gene-counts.csv", row.names = 1)
                     
                      metadata <- data.frame(row.names=colnames(counts), condition=c(rep('Veh',4), rep('STZ',3)))
                      all(rownames(metadata) == colnames(counts))
                      
                      dds <- DESeqDataSetFromMatrix(countData = counts,
                                                    colData = metadata,
                                                    design = ~ condition)
                      dds <- estimateSizeFactors(dds)
                      sizeFactors(dds)
                      
                      #normalized counts_genes
                      normalized_counts_genes <- counts(dds, normalized=TRUE)
                      
                      #log transformation
                      vsd_wt <- vst(dds, blind=TRUE)
                      
                      # Extract the vst matrix from the object
                      vsd_cor_wt <- vsd_wt %>% 
                        assay() %>%
                        cor()
                      
                      # plot heatmap
                      pheatmap(vsd_cor_wt, annotation = select(metadata, condition))
                      
                      # Plot PCA
                      plotPCA(vsd_wt, intgroup="condition")
                      
                      
                      ## Run analysis
                      dds_wt <- DESeq(dds)
                      
                      ## Plot dispersion estimates
                      plotDispEsts(dds_wt)
                      
                      #
                      res <- results(dds_wt,
                                     contrast = c("condition", "STZ",
                                                  "Veh"),
                                     alpha = 0.05) #,lfcThreshold = 0.32
                      #MA plot
                      plotMA(res, ylim=c(-8,8))
                      
                      #shrink fc
                      res <- lfcShrink(dds_wt,
                                       contrast=c("condition", "STZ", "Veh"),
                                       res=res)
                      
                      res_all <- data.frame(res) %>%
                        rownames_to_column(var = "Gene") 
                      
                      #extract significant genes 
                      res_sig <- subset(res_all, padj < 0.05)
                      res_sig <- res_sig %>%
                        arrange(log2FoldChange)
                      
                      #plot sig genes
                      # Subset normalized counts_genes to significant genes
                      sig_norm_counts_genes <- normalized_counts_genes[res_sig$Gene, ]
                      
                      # Run pheatmap
                      pheat <- pheatmap(sig_norm_counts_genes,
                                        color = viridis(100, option = 4),
                                        cluster_rows = T,
                                        show_rownames = T,
                                        annotation = dplyr::select(metadata, condition),
                                        scale = "row")

#TE genes
                      counts <- read.csv("data/4veh-vs-3stz-mavsko-protein-coding-gene+transposable-elements-counts.csv", row.names = "Gene_name")[,-1]
                      
                      metadata <- data.frame(row.names=colnames(counts), condition=c(rep('Veh',4), rep('STZ',3)))
                      all(rownames(metadata) == colnames(counts))
                      
                      dds <- DESeqDataSetFromMatrix(countData = counts,
                                                    colData = metadata,
                                                    design = ~ condition)
                      dds <- estimateSizeFactors(dds)
                      sizeFactors(dds)
                      
                      #normalized counts
                      normalized_counts <- counts(dds, normalized=TRUE)
                      
                      #log transformation
                      vsd_wt <- vst(dds, blind=TRUE)
                      
                      # Extract the vst matrix from the object
                      vsd_cor_wt <- vsd_wt %>% 
                        assay() %>%
                        cor()
                      
                      # plot heatmap
                      pheatmap(vsd_cor_wt, annotation = select(metadata, condition))
                      
                      # Plot PCA
                      plotPCA(vsd_wt, intgroup="condition")
                      
                      ## Run analysis
                      dds_wt <- DESeq(dds)
                      
                      ## Plot dispersion estimates
                      plotDispEsts(dds_wt)
                      
                      #
                      res <- results(dds_wt,
                                     contrast = c("condition", "STZ",
                                                  "Veh"),
                                     alpha = 0.1) #,lfcThreshold = 0.32
                      #MA plot
                      plotMA(res, ylim=c(-8,8))
                      
                      #shrink fc
                      res <- lfcShrink(dds_wt,
                                       contrast=c("condition", "STZ", "Veh"),
                                       res=res)
                      
                      res_all <- data.frame(res) %>%
                        rownames_to_column(var = "Ensembl_Gene_ID") 
                      
                      #extract significant genes 
                      res_sig <- subset(res_all, pvalue < 0.05)
                      res_sig <- res_sig %>%
                        arrange(padj)
                      
                      #plot sig genes
                      # Subset normalized counts to significant genes
                      sig_norm_counts <- normalized_counts[res_sig$Ensembl_Gene_ID, ]
                      sig_norm_counts <- sig_norm_counts[grepl("\\:", rownames(sig_norm_counts)),]
                      
                      # Run pheatmap
                      pheat <- pheatmap(sig_norm_counts,
                                        color = viridis(100),
                                        cluster_rows = T,
                                        show_rownames = T,
                                        annotation = select(metadata, condition),
                                        scale = "row")
                      
                      #donut plots
                      genes <- res_sig %>%
                        filter(grepl("\\:",Ensembl_Gene_ID)) %>%
                        tidyr::separate(Ensembl_Gene_ID, c("geneid", "family", "class"), sep=':') %>%
                        mutate(condition = factor(ifelse(log2FoldChange > 0, "WT", "muMT"), levels = c("WT", "muMT"))) 
                      
                      genes$family <- factor(genes$family, levels = c("ERV1", "ERVK", "ERVL", "ERVL-MaLR", "Gypsy", "L1", "Alu", "B2", "hAT-Charlie",  "hAT-Tip100", "DNA","PiggyBac", "TcMar-Tigger", "UCON26"))
                      
                      genes %>%
                        group_by(condition, family) %>%
                        summarise(freq = length(family)) %>%
                        ggplot(aes(x=2, y=freq,fill=family)) +
                        geom_bar(position = 'fill', stat = 'identity', color='black', lwd=0.1) +
                        coord_polar(theta='y', start=0) +
                        theme_void() +
                        scale_fill_manual(values=colors_many) +
                        facet_wrap(~condition) +
                        xlim(0.5, 2.5)
