###############################################################
# Script: 02_differential_expression.R
# Purpose: Identify differentially expressed genes between experimental groups.
# Inputs:
# 1. Working directory with Count values generated from HTseqcount.
# 2. Name specifying the the generation to be compared (e.g F0_5dpf_)
# Usage:
#   Rscript 02_differential_expression.R \
#       <working_directory> \
#       <name>
## Output:
# Differential expression tables, summary statistics, Figures...
# GSEA and rrvgo enrichment results, summary statistics, Figures...

###################################################
##To ensure that R script accepts command-line arguments using commandArgs().
args <- commandArgs(trailingOnly = TRUE)
count_dir <- args[1] #Directory containing the HTSeq count files
nam <- args[2]   #Name specifying the the generation

##################################################
#count_dir <- "/path_to/.../F0_HTSeq_count"
#nam <- "F0_5dpf_"

#####################
gen <- strsplit(nam, "_")[[1]][1]

################################################################################
## List of required CRAN packages
cran_packages <- c("edgeR","readxl","writexl","sva","ggplot2",
                   "openxlsx","limma","tidyr","pheatmap","DESeq2",
                   "matrixStats","tidyverse","variancePartition",
                  "dplyr","ggvenn")

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")}

## List of required Bioconductor packages
bioc_packages <- c("annotatr","GenomicRanges","genomation","rrvgo",
                   "org.Rn.eg.db","BiocParallel","EnhancedVolcano",
                   "clusterProfiler","enrichplot","AnnotationHub")

## Install missing CRAN packages
missing_cran <- cran_packages[!sapply(cran_packages, requireNamespace, quietly = TRUE)]
if (length(missing_cran) > 0) {
  install.packages(missing_cran)}

## Install missing Bioconductor packages
missing_bioc <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc)}

##################################################
## Load library
library(edgeR)
library(readxl)
library(writexl)
library(sva)
library(ggplot2)
library(openxlsx)
library(limma)
library(tidyr)
library(pheatmap)
library(annotatr)
library(GenomicRanges)
library(genomation)
library("DESeq2")
library("tidyverse")
library('variancePartition')
library('BiocParallel')
library("dplyr")
library("EnhancedVolcano")
library("org.Rn.eg.db", character.only = TRUE)
keytypes(org.Rn.eg.db)

###############################################################
## List all count files
files <- list.files(count_dir, pattern = "\\.counts\\.txt$", full.names = TRUE)

## Read each file
count_list <- lapply(files, function(f) {
  dat <- read.table(f, header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE)
  
  ## Rename count column using filename
  colnames(dat) <- c("gene",
                     sub("\\.counts\\.txt$", "", basename(f)))
  dat
})

# Merge all files by gene ID
count_matrix <- Reduce(function(x, y) merge(x, y, by = "gene"), count_list)

## Remove HTSeq summary rows if present
count_matrix <- count_matrix[!grepl("^__", count_matrix$gene), ]

## Save
write.csv(count_matrix,
          file = paste0(nam,"HTSeq_counts_matrix.csv"),
          row.names = FALSE)

######### SECTION:2 #######################
## IMPLEMENTING SVA
set.seed(1234567)
#filter data by removing all the genes with 0 count across all samples
count_matrix <- count_matrix %>%
  filter(rowSums(.) != 0)
##########
counts <- count_matrix
head(counts, 2)
nrow(counts)
ncol(counts)

#filter data by removing all the genes with less than 10 count across all samples
#counts <- counts %>% filter(rowSums(.) !=0)
counts <- counts %>% filter(rowSums(.) > 10)
nrow(counts)
ncol(counts)
head(counts)

#Make a sample dataframe containing samples info.
################################################################
if (gen == "F0") {
  metadata <- data.frame(
    sample = c("F0_Control1","F0_Control2","F0_Control3",
               "F0_PFBSHC1","F0_PFBSHC2","F0_PFBSHC3",
               "F0_PFBSLC1","F0_PFBSLC2","F0_PFBSLC3",
               "F0_PFOSHC1","F0_PFOSHC2","F0_PFOSHC3",
               "F0_PFOSLC1","F0_PFOSLC2","F0_PFOSLC3"),
    condition = c("Control","Control","Control",
                  "PFBSHC","PFBSHC","PFBSHC",
                  "PFBSLC","PFBSLC","PFBSLC",
                  "PFOSHC","PFOSHC","PFOSHC",
                  "PFOSLC","PFOSLC","PFOSLC"),
    row.names = "sample")
} else if (gen == "F1") {metadata <- data.frame(
  sample = c("F1_Control1","F1_Control2","F1_Control3",
             "F1_PFBSHC1","F1_PFBSHC2","F1_PFBSHC3",
             "F1_PFBSLC1","F1_PFBSLC2","F1_PFBSLC3",
             "F1_PFOSHC1","F1_PFOSHC2","F1_PFOSHC3",
             "F1_PFOSLC1","F1_PFOSLC2","F1_PFOSLC3"),
  condition = c("Control","Control","Control",
                "PFBSHC","PFBSHC","PFBSHC",
                "PFBSLC","PFBSLC","PFBSLC",
                "PFOSHC","PFOSHC","PFOSHC",
                "PFOSLC","PFOSLC","PFOSLC"),
  row.names = "sample")
} else{metadata <- data.frame(
  sample = c("F2_Control1","F2_Control2","F2_Control2",
             "F2_PFBSHC1","F2_PFBSHC2","F2_PFBSHC3",
             "F2_PFBSLC1","F2_PFBSLC2","F2_PFBSLC3",
             "F2_PFOSHC1","F2_PFOSHC2","F2_PFOSHC3",
             "F2_PFOSLC1","F2_PFOSLC2","F2_PFOSLC3"),
  condition = c("Control","Control","Control",
                "PFBSHC","PFBSHC","PFBSHC",
                "PFBSLC","PFBSLC","PFBSLC",
                "PFOSHC","PFOSHC","PFOSHC",
                "PFOSLC","PFOSLC","PFOSLC"),
  row.names = "sample")
}

################################################################
metadata$condition <- as.factor(metadata$condition)
metadata

##Check if all samples are included in both count data and sample dataframe
all(rownames(metadata) == colnames(counts))

# Convert condition to factor (DO NOT include batch because it's unknown)
metadata$condition <- as.factor(metadata$condition) 
metadata

##Check if all samples are included in both count data and sample dataframe
#all(rownames(coldata) %in% colnames(count_matrix)) #or
all(rownames(metadata) == colnames(counts))

counts <- round(counts, digits = 0) #convert data to integer
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, 
                              design = ~ condition)
dds

# Remove genes with low counts
keep <- rowSums(counts(dds)) > 10 #keep only genes with more than 10 reads
dds <- dds[keep,]
dds

##Estimate size factor and plot ##not important.
dds <- estimateSizeFactors(dds)
###############################################################
###############================ SVA PART =================#####
#create a the transformed data matrix with the variables in rows and samples in columns
dat <- counts(dds, normalized=TRUE)
nrow(dat)

#Filter the data further in required
idx <- rowMeans(dat) > 0
dat <- dat[idx,]
nrow(dat)

#Initiate the model matrix to be used to fit the data
mod <- model.matrix(~ condition, colData(dds))
#Initiate the null model to be compared when fitting the data
mod0 <- model.matrix(~ 1, colData(dds))

# Estimate the number of surrogate variables
n.sv <- num.sv(dat, mod, method = "be", B=100,seed = 100000)
#n.sv <- 5
cat("Optimal number of surrogate variables: ", n.sv)

# Run SVA to estimate surrogate variables
svseq <- svaseq(dat, mod, mod0, n.sv = n.sv,B=100)
# View surrogate variables
svseq$sv

# Add SVs to metadata (colData)
sv_columns <- paste0("SV", 1:n.sv)  # Create column names: SV1, SV2, ...

ddssva <- dds
# Add surrogate variables to ddssva
for (i in seq_len(n.sv)) {
  ddssva[[sv_columns[i]]] <- svseq$sv[, i]
}

# Construct the design formula dynamically
design_formula <- as.formula(paste("~", paste(sv_columns, collapse = " + "), " + condition"))

# Assign the design formula to the DESeq2 object
design(ddssva) <- design_formula

# Run DESeq2
ddssva <- DESeq(ddssva)
#view result
(ressva <- results(ddssva))
res_af <- results(ddssva,c("condition","PFBSHC","Control"))

mcols(ressva, use.names=TRUE) ##Get meaning of the res columns
ressva.05 <- results(ddssva, alpha=.05)
table(ressva.05$padj < .05)

### Each group VS Control
sample <- list("PFBSHC","PFBSLC","PFOSHC","PFOSLC")
for (group in sample) {
  ## Using Contrast to obtain sample pairs
  ressva <- results(ddssva, contrast = c('condition', group, 'Control'))
  #ressva <- ressva[order(ressva$padj),] #sort result by padj
  summary(ressva)
  
  #Examine read count of a single gene
  plotCounts(ddssva, gene=which.min(ressva$padj), intgroup="condition")
  
  ## Map Ensembl to gene symbol and make dataframe with all features
  library(org.Dr.eg.db)
  ressva.df <- as.data.frame(ressva)
  ressva.df$gene_symbol <- mapIds(org.Dr.eg.db, key = row.names(ressva.df), keytype = "ENSEMBL", column = "SYMBOL")
  nrow(ressva.df)
  head(ressva.df)
  write.csv(as.data.frame(ressva.df), file=paste0(nam,"DEG_",group,"_vs_Control_sva_all.csv"))
  
  #save all the result from the analysis
  res_ordered <- ressva.df[order(ressva.df$padj),] #sorting data using the pvalue
  write.csv(as.data.frame(res_ordered), file=paste0(nam,"DEG_",group,"_vs_Control_sva_all_order.csv"))
  
  ## View and save only the significant DEGs
  resFilt <- ressva.df[which(ressva.df$padj < 0.05 & abs(ressva.df$log2FoldChange) > 0.0), ]
  write.csv(resFilt, file=paste0(nam,"DEG_",group,"_vs_Control_sva_sig.csv")) #OR
  
  resOrdered <- resFilt[order(resFilt$pvalue),] #sorting data using the pvalue
  #write.csv(resOrdered, file=paste0(nam,"DEG_",group,"_vs_Control_sva_sig_order.csv")) #OR
  
  ## Count the number of DEG using padj
  sum(ressva.df$padj < 0.05, na.rm=T) #How many adjusted p-values were less than 0.05
  
  ## Count the No. of upreg & No. of downreg
  resSigind = ressva.df[ which(ressva.df$padj < 0.05 & ressva.df$log2FoldChange > 0), ] #upregulated or induced
  resSigrep = ressva.df[ which(ressva.df$padj < 0.05 & ressva.df$log2FoldChange < 0), ] #downregulated or repressed
  resSig = rbind(resSigind, resSigrep)
  dim(resSigind)
  dim(resSigrep)
  head(resSigind, 3)
  head(resSigrep, 3)
  
  ## Adjusted p value histogram
  hist(ressva.df$padj, 
       breaks = seq(0, 1, 0.05), 
       col = "purple1",
       main = paste("Histogram of padj for", group,"vs Control"),
       xlab = "padj")
  
  ## MAKE VOLCANO PLOTS
  #Determine point colors based on significance and sign of the logFC
  #Use adj.P.Value instead of P.Value
  ressva.df <- ressva.df %>% 
    mutate(point_color = case_when(
      padj < 0.05 & log2FoldChange < 0 ~ "down", # significantly down
      padj < 0.05 & log2FoldChange > 0 ~ "up", # significantly up
      TRUE ~ "NS") # not significant
    )
  
  ## Create the Volcano plot without lable
  ##https://biostatsquid.com/volcano-plots-r-tutorial/
  volcano_plot <- ggplot(data = ressva.df, aes(x = log2FoldChange, y = -log10(padj), col = point_color)) +
    #geom_vline(xintercept = c(-1, 1), col = "green4", linetype = 'dashed')+
    geom_hline(yintercept = -log10(0.05), col = "green4", linetype = 'dashed')+
    #geom_point(size = 2)+
    geom_point(color = ifelse(ressva.df$log2FoldChange < -0.001 & ressva.df$padj<0.05, "slateblue", ifelse(ressva.df$log2FoldChange > 0.001 & ressva.df$padj<0.05, "lightpink3", "grey")), size = 1.5) +
    geom_text(aes(label = ifelse(ressva.df$log2FoldChange < -2 & -log10(padj) >20 | ressva.df$log2FoldChange > 2 & -log10(padj) >20, ressva.df$gene_symbol, ""), size = 3),
              vjust = -0.5, hjust = 0, size = 3)+
    scale_color_manual(values = c("slateblue", "grey", "lightpink"),
                       labels = c("Down", "NS", "Up"))+
    #coord_cartesian(ylim = c(0, 5), xlim = c(-10, 10))+
    labs(color = 'Annotation', #legend_title, 
         x = expression("log"[2]*"FC"), y = expression("-log"[15]*"Adjp-value"), size = 5)+
    #scale_x_continuous(breaks = seq(-10, 10, 2)) +
    ggtitle(paste('F2:',group,'vs Control')) +
    guides(shape = guide_legend(override.aes = list(color = NULL, size = 4)))+
    scale_shape_manual(values = c(16, 17, 18))+
    #guides(color = guide_legend(override.aes = list(shape = NA))) +
    theme(axis.line = element_line(),panel.background = element_blank())
  #theme(panel.background = element_blank(),legend.position = "top", legend.justification = "right") #change legend location
  print(volcano_plot)
  ## Export the plot to a file (e.g., PNG)
  ggsave(paste0(nam,"DEG_",group,"_vs_control_volcano_plot_sva.png"), plot = volcano_plot, width = 8, height = 6, dpi = 300)
  
  # VOLCANO
  # Add a column to indicate significant genes with specific labels
  ressva.df$Significant <- with(ressva.df, ifelse(padj < 0.05 & log2FoldChange > 0, "Upregulated",
                                                  ifelse(padj < 0.05 & log2FoldChange < -0, "Downregulated", "Not Significant")))
  ##########################
  ##https://biostatsquid.com/volcano-plots-r-tutorial/
  library(ggplot2)
  library(ggrepel)
  
  ressva.df$Significance <- "Not Significant"
  ressva.df$Significance[ressva.df$padj < 0.05 & ressva.df$log2FoldChange > 0] <- "Upregulated"
  ressva.df$Significance[ressva.df$padj < 0.05 & ressva.df$log2FoldChange < 0] <- "Downregulated"
  ################
  ressva.df$negLogP <- -log10(ressva.df$padj)
  ##############
  top_up <- ressva.df %>% filter(log2FoldChange > 0) %>% arrange(padj) %>% slice(1:5)
  top_down <- ressva.df %>% filter(log2FoldChange < 0) %>% arrange(padj) %>% slice(1:5)
  top_genes <- bind_rows(top_up, top_down)
  ##################
  volcano_plot <- ggplot(ressva.df, aes(x = log2FoldChange, y = negLogP, color = Significance)) +
    geom_point(alpha = 0.7) +
    geom_text_repel(data = top_genes,aes(label = gene_symbol),size = 4, max.overlaps = Inf) +
    scale_color_manual(values = c("Upregulated" = "#e66101","Downregulated" = "#5e3c99",
                                  "Not Significant" = "grey"
    )) +
    coord_cartesian(ylim = c(0.0, 10), xlim = c(-10,10))+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_classic() +
    theme_minimal() +
    theme(
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      axis.title.x = element_text(size = 18),  # x-axis title
      axis.title.y = element_text(size = 18),
      axis.text.x  = element_text(size = 14),  # x-axis ticks
      axis.text.y  = element_text(size = 14),   # y-axis ticks
      legend.position = "none"
    ) +
    labs(title = paste0(gen,' :',group,' vs Control'), x = "log2 Fold Change",y = "-log10(Adjusted P-value)")
  print(volcano_plot)
  # Export the plot to a file (e.g., PNG)
  ggsave(paste0(nam,"DEG_",group,"_vs_control_volcano_plot2_sva.png"), plot = volcano_plot, width = 8, height = 6, dpi = 300)
  
  ##############################################################
  ###=============== GSEA Analysis ==========================###
  library("clusterProfiler")
  library("enrichplot")
  library(AnnotationHub)
  
  ## Recall all differential expression dataframe (res_ordered)
  diff = res_ordered
  head(diff, 2)
  
  colnames(diff)[1] <- "Gene_stable_ID" #change column name
  head(diff)
  
  # Create a named vector for GSEA
  ranked_genes <- setNames(diff$stat, diff$Gene_stable_ID)
  head(ranked_genes)
  
  # Ensure the vector is sorted in decreasing order
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  head(ranked_genes)
  
  # Do enrichment abalysis
  ##GENENAME for gene symbol or ENSEMBL for ensembl
  gse <- gseGO(geneList=ranked_genes, 
               ont ="BP", 
               keyType = "ENSEMBL",  
               nPerm = 10000, 
               minGSSize = 5, 
               maxGSSize = 500, 
               verbose = TRUE,
               pvalueCutoff=0.1,
               eps = 0,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH")
  require(DOSE)
  dot_plot <- dotplot(gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
  dot_plot
  ggsave(paste0(nam,"DEG_",group,"_GSEA_stat.png"), plot = dot_plot, 
         width = 12, height = 6, dpi = 400)
  
  ##### save to dataframe
  enrich_pathways<-as.data.frame(gse) #convert object to dataframe
  nrow(enrich_pathways)
  write_xlsx(enrich_pathways, paste0(nam,"DEG_",group,"_GSEA_stat.xlsx"))
         
  ##sig terms
  sig_bcg = enrich_pathways[which(enrich_pathways$p.adjust < 0.05),]
  nrow(sig_bcg)
  write_xlsx(sig_bcg, paste0(nam,"DEG_",group,"sig_GSEA_stat.xlsx"))
  
  ########################################################
  ####================== RRVGO =========================##
  sig_bp <- sig_bcg
  library(rrvgo)
  simMatrix <- calculateSimMatrix(sig_bp$ID,
                                  orgdb="org.Dr.eg.db",
                                  ont="BP",
                                  method="Rel")
  head(simMatrix)
  
  ######################################
  sig_bp$p.adjust <- as.numeric(sig_bp$p.adjust) #convert to numeric
  scores <- setNames(-log10(sig_bp$p.adjust), sig_bp$ID)
  head(scores)
  
  ######################################
  reducedTerms <- reduceSimMatrix(
    simMatrix,
    scores,
    threshold = 0.7,
    orgdb="org.Dr.eg.db"
  )
  head(reducedTerms)
  
  ## check for overlaps between the IDs in the enriched GO and RRvGO 
  library(ggvenn)
  set1 <- unique(sig_bp$ID)
  set2 <- unique(reducedTerms$go)
  
  ## Combine into a list for ggvenn
  mylist <- list(
    EnrichedGO = set1,
    RRVGO = set2
  )
  # Draw Venn diagram
  ggvenn(mylist, fill_color = c("skyblue", "pink"))
  
  ### Merge the two data enrichedGO and RRVGO data
  merged_GO <- merge(sig_bp, reducedTerms, by.x = "ID", by.y = "go", all.x = TRUE)
  write_xlsx(merged_GO, paste0(nam,"DEG_",group,"_RRVGO.xlsx"))
  
  ### The representative term is the one with:termDispensability = 0
  repTerms <- merged_GO %>% filter(termDispensability == 0)
  dim(repTerms)
  length(unique(repTerms$ID))
  write_xlsx(repTerms, paste0(nam,"DEG_",group,"_RRVGO_repterms.xlsx"))
}

#######################################################
## ================ BEFORE SVA ====================####
# Perform variance stabilizing transformation (VST)
vsd_before <- vst(ddssva, blind = TRUE)
# Compute PCA before batch correction
pca_before <- prcomp(t(assay(vsd_before)))
# Prepare data for plotting
df_before <- data.frame(
  PC1 = pca_before$x[, 1],  # First principal component
  PC2 = pca_before$x[, 2],  # Second principal component
  Condition = metadata$condition  # Assign colors based on condition groups
)
# Define custom colors for the 5 conditions
condition_colors <- c("PFOSHC" = "red", "PFOSLC" = "blue", "PFBSLC" = "green", 
                      "PFBSHC" = "purple", "Control" = "orange")
# Define distinct shapes for each group
group_shapes <- c("Control" = 16, # filled circle
  "PFOSLC" = 17, # triangle
  "PFOSHC" = 15, # square
  "PFBSLC" = 18, # diamond
  "PFBSHC" = 8 # star
)
# Calculate explained variance
var_explained <- round(100 * summary(pca_before)$importance[2, 1:2], 1)

# Scatter plot before batch correction
sva_bf <- ggplot(df_before, aes(x = PC1, y = PC2, color = Condition, shape = Condition)) +
  geom_point(size = 5, alpha = 0.9) +
  scale_shape_manual(values = group_shapes) +
  theme_minimal(base_size = 20) +
  labs(title = paste0(gen,": Before Batch Correction"),
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  #scale_color_manual(values = condition_colors) +  # Assign colors
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5))
print(sva_bf)
ggsave(paste0(nam,"DEG_",group,"_vs_control_uncorrected.png"), plot = sva_bf,
       width = 8, height = 6.5, dpi = 300)

###############################################################
##================= AFTER SVA ===========================######
# Remove batch effects using SVA surrogate variables
corrected_counts <- removeBatchEffect(assay(vsd_before), covariates = svseq$sv,
                                      group = metadata$condition)
# Compute PCA after batch correction
pca_after <- prcomp(t(corrected_counts),center = TRUE)
# Prepare data for plotting
df_after <- data.frame(
  PC1 = pca_after$x[, 1],
  PC2 = pca_after$x[, 2],
  Condition = metadata$condition  # Assign colors based on condition groups
)

# Calculate explained variance
var_explained <- round(100 * summary(pca_after)$importance[2, 1:2], 1)

# Scatter plot after batch correction
sva_af <- ggplot(df_after, aes(x = PC1, y = PC2, color = Condition, shape = Condition)) +
  geom_point(size = 5, alpha = 0.9) +
  scale_shape_manual(values = group_shapes) +
  theme_minimal(base_size = 20) +
  labs(title = paste0(gen,": SVA Corrected (SV = 1)"),
       
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  #scale_color_manual(values = condition_colors) +  # Assign colors
  theme(
    #legend.position = "bottom",
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
    )
print(sva_af)
ggsave(paste0(nam,"DEG_",group,"_vs_control_sva_corrected.png"), plot = sva_af,
       dpi = 400, width = 8, height = 6.5)