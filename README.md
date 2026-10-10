# Nanopore-Methylation-RNASeq-pipeline
Reproducible workflow to perform an integrated multi-omics analysis of RNA sequencing (RNA-seq) and Oxford Nanopore Technologies (ONT) sequencing data. The workflow includes differential gene expression analysis, DNA methylation analysis, functional enrichment, methylation-expression correlation analysis, association of methylation with phenotypic outcomes, and downstream visualization.

The analyses were performed primarily in **R**, with ONT basecalling, alignment, and methylation extraction performed using **Dorado** and **Modkit**. Further downstream visualization was performed in **Python**

---

## Workflow

### Part 1: RNA-seq Differential Expression Analysis
#### Section 1: Indexing, Alignment, and Counting Pipeline
*The repository includes a shell script (.sh) to automate the RNA-seq processing workflow, including:
- Indexing: Build the reference genome index
- Alignment: Align sequencing reads to the reference genome
- Counting: Quantify gene-level read counts for downstream analysis.
*Ensure that the required tools, reference genome, and annotation files are available before running the script.
RNA-seq count data were analyzed using **DESeq2** to identify differentially expressed genes between experimental groups.

#### Section 2:  Differential Expression Analysis Steps
The 02_differential_expression.R script performs the following:
- Import raw gene count matrix
- Perform normalization and quality control
- Differential expression analysis using DESeq2
- Identify significantly differentially expressed genes
- Perform Gene Set Enrichment Analysis (GSEA)
- Summarize redundant Gene Ontology (GO) terms using RRvGo

---
### Part 2: Oxford Nanopore DNA Methylation Analysis
Oxford Nanopore sequencing data were processed to identify differentially methylated genomic positions.

#### Section 1: Dorado Basecalling and 5mC Detection
The 03A_dorado.slum script performs basecalling, alignment, and 5-methylcytosine (5mC) detection in the CpG context using Dorado.
##### Requirements
Run the script from the directory containing the following files and folders:
- dorado (Dorado binary)
- dorado_models/ (Dorado basecalling and modification models)
- pod5/ (input POD5 files)
- reference.fa (reference genome in FASTA format)
#### Output
- Basecalled reads with 5mC modification calls in the CpG context.
- Reads aligned to the reference genome.
- Sorted BAM file for downstream analysis.

#### Section 2: Differential methylation analysis
The 03B_modkit.slum script performs pileup generation and differential methylation analysis using Modkit.
##### Requirements
Run the script with the directory containing the Sorted, aligned BAM files (.bam).
##### Output
- Pileup file for each exposure in each generation.
- Differential methylation results for each exposure in each generation.

#### Section 3: Annotation of the methylated positions
- Annotate methylated positions to genomic features (e.g., promoter, exon, intron, gene body)
#### Section 4: Functional enrichment analysis
- Perform functional enrichment analysis of annotated genes
- Summarize enriched GO terms using RRvGo

#### Software
- Dorado
- Modkit
- clusterProfiler
- RRvGo
- R

---
### Part 3: Correlation Between DNA Methylation and Gene Expression
The relationship between DNA methylation and transcriptional activity was evaluated using Spearman's rank correlation.

#### Analysis Steps
- Match methylation measurements with corresponding gene expression values
- Compute Spearman's rank correlation coefficient
- Assess statistical significance of correlations

#### Statistical Method
- Spearman's rank correlation

---
### Part 4: Association Between DNA Methylation and Phenotypic Outcomes
DNA methylation levels were tested for association with phenotypic variables using beta regression.

#### Analysis Steps
- Fit beta regression models using methylation proportions as the response variable
- Identify significantly associated methylation sites
- Annotate significant methylation sites to genomic features
- Perform Over-Representation Analysis (ORA)
- Summarize enriched Gene Ontology terms using RRvGo

---
## Repository Structure

```text
.
├── RNAseq/
│   ├── Differential_expression/
│   ├── GSEA/
│   └── RRvGo/
│
├── ONT_Methylation/
│   ├── Dorado/
│   ├── Modkit/
│   ├── Annotation/
│   ├── Enrichment/
│   └── RRvGo/
│
├── Correlation/
│   └── Spearman/
│
├── Phenotype_Association/
│   ├── BetaRegression/
│   ├── Annotation/
│   ├── ORA/
│   └── RRvGo/
│
├── scripts/
├── results/
└── README.md
```
### R Packages
- DESeq2, clusterProfiler, RRvGo, betareg, dplyr, tidyr, data.table, ggplot2, GenomicRanges, AnnotationDbi...
---

## Analysis Summary

| Analysis | Method | Output |
|----------|--------|--------|
| Differential expression | DESeq2 | Differentially expressed genes |
| Functional enrichment | GSEA | Enriched biological pathways |
| GO term summarization | RRvGo | Non-redundant GO terms |
| DNA methylation analysis | Dorado + Modkit | Differentially methylated positions |
| Methylation annotation | R | Gene feature annotations |
| Functional enrichment | clusterProfiler | Enriched biological processes |
| Methylation-expression integration | Spearman correlation | Correlation statistics |
| Phenotype association | Beta regression | Associated methylation sites |
| Functional interpretation | ORA + RRvGo | Summarized enriched GO terms |

---

## Reproducibility

The analyses are organized into modular workflows corresponding to each stage of the study. Intermediate outputs generated by one module serve as inputs for downstream analyses. Each script can be executed independently, provided that the required input files and software dependencies are available.

---

## Citation

If you use this repository or adapt parts of this workflow, please cite the original software packages used in the analysis:

- Love MI, Huber W, Anders S. **DESeq2: Moderated estimation of fold change and dispersion for RNA-seq data.**
- Yu G, Wang LG, Han Y, He QY. **clusterProfiler: An R package for comparing biological themes among gene clusters.**
- Sayols S. **RRvGo: A Bioconductor package to reduce and visualize Gene Ontology terms.**
- Oxford Nanopore Technologies. **Dorado.**
- Oxford Nanopore Technologies. **Modkit.**
- Cribari-Neto F, Zeileis A. **Beta Regression in R.**

---

## Author

**Ogunleye Adeolu**

**Uppsala University, Sweden**

---
