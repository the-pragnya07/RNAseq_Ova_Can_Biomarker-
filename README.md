# RNAseq_Ova_Can_Biomarker

## Overview
This project focuses on transcriptomic analysis of ovarian cancer using RNA-seq data to identify differentially expressed genes, co-expression networks, enriched biological pathways, and potential biomarker candidates associated with ovarian cancer progression.

The study integrates differential expression analysis, batch correction, visualization, WGCNA, and functional enrichment analysis using R and Bioconductor packages.

---

# Objectives

- Identify differentially expressed genes (DEGs) in ovarian cancer
- Compare gene expression between:
  - Early Stage vs Normal
  - Late Stage vs Normal
  - Early Stage vs Late Stage
- Perform batch correction and normalization
- Construct co-expression gene networks using WGCNA
- Identify hub genes associated with ovarian cancer progression
- Perform GO and KEGG enrichment analysis
- Discover potential ovarian cancer biomarkers

---

# Workflow

```text
TCGA Ovarian RNA-seq Dataset
              ↓
Data Cleaning & Preprocessing
              ↓
Differential Expression Analysis (DESeq2)
              ↓
Visualization & Quality Assessment
              ↓
WGCNA Network Construction
              ↓
GO & KEGG Pathway Enrichment
              ↓
Hub Gene Identification
              ↓
Biomarker Discovery
```

---

# Dataset

- Source: TCGA Ovarian Cancer RNA-seq Dataset
- Data Type: Transcriptomic RNA-seq Data
- Disease: Ovarian Cancer

---

# Methodology

## 1. Data Cleaning & Preprocessing
- Removal of low-expression genes
- Data normalization
- Expression matrix preparation
- Sample quality assessment

---

## 2. Differential Expression Analysis
Differential expression analysis was performed using the DESeq2 package in R.

### Comparisons Performed
- Early Stage vs Normal
- Late Stage vs Normal
- Early Stage vs Late Stage

### Criteria
- Adjusted p-value < 0.05
- |log2FoldChange| > 1

---

## 3. Visualization
Several plots were generated to visualize differential expression patterns and sample clustering.

### Visualizations Included
- Volcano plots
- Heatmaps
- Venn diagrams
- WGCNA plots
- KEGG enrichment plots

---

## 4. WGCNA Analysis
Weighted Gene Co-expression Network Analysis (WGCNA) was used to identify gene modules and hub genes associated with ovarian cancer progression.

### WGCNA Steps
- Soft threshold selection
- Adjacency matrix construction
- TOM calculation
- Module identification
- Hub gene extraction

---

## 5. Functional Enrichment Analysis

### GO Analysis
- Biological Process (BP)
- Molecular Function (MF)
- Cellular Component (CC)

### KEGG Analysis
Pathway enrichment analysis was performed to identify significant biological pathways involved in ovarian cancer.

---

# Tools & Packages

## Programming Language
- R

## R Packages Used
- DESeq2
- WGCNA
- clusterProfiler
- ggplot2
- pheatmap
- EnhancedVolcano
- edgeR
- org.Hs.eg.db

---

# Repository Structure

```text
RNAseq_Ova_Can_Biomarker/
│
├── scripts/
│   ├── Deseq.R
│   ├── Kegg.R
│   ├── WGCNA.R
│   └── very_clean.R
│
├── figures/
│   ├── E_VS_L_PLOT.png
│   ├── E_vs_N_plot.png
│   ├── EARLY_STAGE_HUB.png
│   ├── Early_vs_Normal_kegg.png
│   ├── Early_vs_Normal_modheatmap.png
│   ├── L_VS_N_PLOT.png
│   ├── LATE_STgE_HUB_GENES.png
│   ├── Late_vs_Normal_kegg.png
│   ├── Late_vs_Normal_modheatmap.png
│   ├── softpower.png
│   └── venn27.png
│
├── Results/
│   ├── early_only_DEGs.xls
│   ├── early_vs_late_only_DEGs.xls
│   ├── Early_vs_Normal_all_genes.xls
│   ├── late_only_DEGs.xls
│   ├── Late_vs_Normal_all_genes.xls
│   └── tcga_genes_common.xls
│
├── docs/
│   └── dissertation_MSc.pdf
│
├── README.md
└── .gitignore
```

---

# Figures

## Volcano Plots
Volcano plots were generated to visualize significantly upregulated and downregulated genes.

### Included Plots
- Early vs Late
- Early vs Normal
- Late vs Normal

---

## Heatmaps
Heatmaps were used to visualize expression patterns of significant genes and module relationships.

---

## KEGG Enrichment Plots
KEGG pathway enrichment plots highlight pathways associated with ovarian cancer progression and dysregulated biological processes.

---

## WGCNA Plots
WGCNA visualizations include:
- Soft-threshold power selection
- Module heatmaps
- Hub gene identification

---

## Venn Diagram
The Venn diagram represents common and unique DEGs among different comparisons.

---

# Results

## Differentially Expressed Genes
Significant DEGs were identified across all comparisons using DESeq2 analysis.

### Result Files
- early_only_DEGs.xls
- early_vs_late_only_DEGs.xls
- Early_vs_Normal_all_genes.xls
- late_only_DEGs.xls
- Late_vs_Normal_all_genes.xls

---

## Common Genes
Shared genes across comparisons were identified and stored in:

```text
tcga_genes_common.xls
```

---

## Hub Genes
Potential hub genes associated with ovarian cancer progression were identified using WGCNA.

---

# Key Findings

- Identification of stage-specific dysregulated genes
- Detection of co-expression modules associated with ovarian cancer
- Discovery of significant pathways related to cancer progression
- Identification of hub genes as potential biomarkers
- Detection of common DEGs across multiple comparisons

---

# Future Directions

- Clinical validation of hub genes
- Survival analysis using clinical metadata
- Machine learning-based biomarker prediction
- Multi-omics integration
- Drug-target interaction studies

---

# Applications

This project demonstrates applications of:
- Bioinformatics
- Transcriptomics
- Cancer Genomics
- Computational Biology
- Biomarker Discovery
- Network Biology

---

# Author

Kollapur Dhaksha Pragnya Sree  
M.Sc Bioinformatics  
Pondicherry University

---

# Acknowledgements

I sincerely thank my project guide Dr.Ayaluru Murali, faculty members, research scholars Ahana Roy Choudhary and Lukkani Laxman Kumar and laboratory members for their continuous guidance and support throughout this research work.

---

# License

This project is licensed under the MIT License.
