library(TCGAbiolinks)
library(SummarizedExperiment)
library(recount3)

library(limma)
library(pheatmap)
library(DESeq2)
library(dplyr)
if (!require("variancePartition"))
  BiocManager::install("variancePartition")
library(variancePartition)
#Download TCGA-OV Tumor count
query <- GDCquery(
project = "TCGA-OV",
data.category = "Transcriptome Profiling",
data.type = "Gene Expression Quantification",
workflow.type = "STAR - Counts",
sample.type = "Primary Tumor",
access = "open"
)

GDCdownload(query)
ov_data  <- GDCprepare(query)
 #Extract Counts + Metadata
tcga_counts <- assay(ov_data , "unstranded")
tcga_meta   <- as.data.frame(colData(ov_data))
 tcga_genes  <- as.data.frame(rowData(ov_data))
 
#Remove Ensembl versions
rownames(tcga_counts) <- sub("\\..*", "", rownames(tcga_counts))
tcga_genes$gene_id <- sub("\\..*", "", tcga_genes$gene_id)
#4.Keep Protein-Coding Genes Only
protein_coding <- tcga_genes$gene_type == "protein_coding"

tcga_counts <- tcga_counts[protein_coding, ] #19k genes
tcga_genes  <- tcga_genes[protein_coding, ]
tcga_counts
dim(tcga_genes)
View(tcga_genes)
saveRDS(tcga_counts, "tcga_counts.rds")
saveRDS(ov_data, "ov_data.rds")

cat("Genes after protein coding filter:", nrow(tcga_counts), "\n")
# 5.Remove samples with tissue not reported
valid_tissue <- tcga_meta$tissue_or_organ_of_origin != "Not Reported"
tcga_meta   <- tcga_meta[valid_tissue, ]
tcga_counts <- tcga_counts[, valid_tissue]

cat("Samples after tissue filter:", ncol(tcga_counts), "\n")
View(tcga_counts)

# Clean stage text
tcga_meta$stage_clean <- toupper(gsub(" ", "", tcga_meta$figo_stage))

# Create Early vs Late stage groups
tcga_meta$stage_group <- ifelse(
  tcga_meta$stage_clean %in% c("STAGEIA","STAGEIB","STAGEIC",
                               "STAGEIIA","STAGEIIB","STAGEIIC"),
  "Early",
  ifelse(
    tcga_meta$stage_clean %in% c("STAGEIIIA","STAGEIIIB","STAGEIIIC","STAGEIV"),
    "Late",
    NA
  )
)

table(tcga_meta$figo_stage, useNA="always")
# Remove samples without stage
valid <- !is.na(tcga_meta$stage_group)

tcga_meta <- tcga_meta[valid, ]
tcga_counts <- tcga_counts[, valid]
table(tcga_meta$stage_clean)
#Early vs Late count

stage_counts <- table(tcga_meta$stage_group)

print(stage_counts)
#print total numbers

cat("Total Early Stage samples (IIA,IIB,IIC):", stage_counts["Early"], "\n")
cat("Total Late Stage samples (IIIA,IIIB,IIIC,IV):", stage_counts["Late"], "\n")

#See full table with stages
View(tcga_meta[, c("barcode","figo_stage","stage_clean","stage_group")])



View(tcga_meta[, c("barcode", "figo_stage", "stage_group")])
#Keep only important metadata columns
tcga_meta_clean <- tcga_meta[, c(
  "barcode",
  "stage_group"
)]
View(tcga_meta_clean)
#Rename barcode to sample
colnames(tcga_meta_clean)[colnames(tcga_meta_clean) == "barcode"] <- "sample"
View(tcga_meta_clean)

#Add condition and dataset labels
tcga_meta_clean$condition <- "Tumor"
tcga_meta_clean$dataset <- "TCGA"   #is it necessary to keep dataset? should i remove

#Ensure sample order matches count matrix
tcga_meta_clean <- tcga_meta_clean[match(colnames(tcga_counts), tcga_meta_clean$sample), ]

#Check stage distribution
table(tcga_meta_clean$stage_group)

#Check total samples
cat("Total TCGA samples:", nrow(tcga_meta_clean), "\n")

# 9.Count samples in each stage group
stage_counts <- table(tcga_meta$stage_group)
print(stage_counts)
View(tcga_meta_clean)
View(tcga_counts)

cat("Early Stage samples:", stage_counts["Early"], "\n")
cat("Late Stage samples:", stage_counts["Late"], "\n")

tcga_counts <- tcga_counts[rowSums(tcga_counts) > 0, ]
View(tcga_counts)
sum(rowSums(tcga_counts) == 0)
dim(tcga_genes)
# 10.Filter low count genes
keep_genes <- rowSums(tcga_counts >= 10) >= (0.1 * ncol(tcga_counts))
final_counts <- tcga_counts[keep_genes, ]
cat("Genes after low count filtering:", nrow(final_counts), "\n") #17707

sum(rowSums(tcga_counts) == 0)
#Add condition + batch:
tcga_meta$condition <- "Tumor"
tcga_meta$dataset   <- "TCGA"
View(tcga_meta$condition)


tcga_ids_clean <- sub("\\..*", "", rownames(tcga_genes))
tcga_genes <- tcga_genes[!duplicated(tcga_ids_clean), ]
tcga_ids_clean <- tcga_ids_clean[!duplicated(tcga_ids_clean)]
rownames(tcga_genes) <- tcga_ids_clean
any(duplicated(rownames(tcga_genes)))
count_ids_clean <- sub("\\..*", "", rownames(final_counts))
final_counts <- final_counts[!duplicated(count_ids_clean), ]
count_ids_clean <- count_ids_clean[!duplicated(count_ids_clean)]
rownames(final_counts) <- count_ids_clean

length(intersect(rownames(tcga_genes), rownames(final_counts)))

tcga_genes <- tcga_genes[!duplicated(rownames(tcga_genes)), ]
final_counts <- final_counts[!duplicated(rownames(final_counts)), ]
common_genes <- intersect(rownames(tcga_genes), rownames(final_counts))
tcga_genes_common   <- tcga_genes[common_genes, ]
final_counts_common <- final_counts[common_genes, ]
# 12.Save final files
dim(final_counts_common)
dim(tcga_genes_common)

View(final_counts)
View(tcga_genes_common)
# 12.Save final files
write.csv(tcga_meta_clean, "tcga_metadata_clean.csv")
write.csv(final_counts, "tcga_counts_clean.csv")
write.csv(tcga_genes_common, "tcga_genes_common.csv")
View(tcga_genes_common)
View(final_counts)


library(recount3)

# -------------------------------------------
# Download GTEx Gene-level data
# -------------------------------------------
#not important start


# 1. Download/Create the RSE object
gtex_rse <- create_rse_manual(
  project      = "OVARY",
  project_home = "data_sources/gtex",
  organism     = "human",
  annotation   = "gencode_v26",
  type         = "gene"
)




# 1. Reset the full object with counts (ensure it has 195 samples)
assay(gtex_rse, "counts") <- compute_read_counts(gtex_rse)

# 2. Define your protein coding filter based on the ORIGINAL rowData
# rowData(gtex_rse) contains the gene types for all 63k genes
pcg_mask <- rowData(gtex_rse)$gene_type == "protein_coding"

# 3. Define your tissue filter (Ovary)
tissue_mask <- colData(gtex_rse)$gtex.smtsd == "Ovary"

# 4. SUBSET BOTH AT ONCE: [rows, columns]
# This keeps the 22k protein coding genes AND the 195 Ovary samples
gtex_final_rse <- gtex_rse[pcg_mask, tissue_mask]

# 5. EXTRACT the synced components
gtex_counts <- assay(gtex_final_rse, "counts")
gtex_genes  <- as.data.frame(rowData(gtex_final_rse))

# 6. VERIFY: Should be [22321, 195] and [22321, 10]
dim(gtex_counts)
dim(gtex_genes)


# 6. FINAL CHECK: This should no longer give a "subscript out of bounds" error
print(dim(gtex_counts))
head(gtex_counts[, 1:3])


# Remove Ensembl versions
rownames(gtex_counts) <- sub("\\..*", "", rownames(gtex_counts))
gtex_genes$gene_id    <- sub("\\..*", "", gtex_genes$gene_id)
# 1. Create a logical mask for genes with at least 10 reads in at least 20 samples
# Adjust '10' and '20' based on how strict you want to be
sum(rowSums(gtex_counts) == 0)
gtex_counts <- gtex_counts[rowSums(gtex_counts) > 0, ]
#2️⃣ Filter low-count genes
keep_genes <- rowSums(gtex_counts >= 10) >= (0.1 * ncol(gtex_counts))
gtex_counts_filtered <- gtex_counts[keep_genes, ]
View(gtex_counts_filtered)
#3️⃣ Keep only required columns from annotation
gtex_genes_small <- gtex_genes[, c("gene_id", "gene_type", "gene_name")]
#Remove version numbers from count matrix
rownames(gtex_counts_filtered) <- sub("\\..*", "", rownames(gtex_counts_filtered))
#Match annotation to filtered counts
gtex_genes_filtered <- gtex_genes_small[
  gtex_genes_small$gene_id %in% rownames(gtex_counts_filtered),
]
#Reorder annotation to match counts
gtex_genes_filtered <- gtex_genes_filtered[
  match(rownames(gtex_counts_filtered), gtex_genes_filtered$gene_id),
]
all(rownames(gtex_counts_filtered) == gtex_genes_filtered$gene_id)
# 4. Check how many genes are left
dim(gtex_counts_filtered)
dim(gtex_genes_filtered)
View(gtex_genes_filtered)

# Remove version numbers
gtex_ids_clean  <- sub("\\..*", "", rownames(gtex_genes_filtered))
count_ids_clean <- sub("\\..*", "", rownames(gtex_counts_filtered))

# Assign back
rownames(gtex_genes_filtered)  <- gtex_ids_clean
rownames(gtex_counts_filtered) <- count_ids_clean

head(rownames(gtex_genes_filtered))
head(rownames(gtex_counts_filtered))




gtex_genes <- as.data.frame(rowData(gtex_final_rse))
write.csv(
  gtex_counts_filtered,
  file = "gtex_counts_filtered.csv")
write.csv(
  gtex_genes_filtered,
  file = "gtex_genes_filtered.csv")


# 6.FINAL VERIFICATION
print(paste("Samples found:", ncol(gtex_counts_filtered)))
print(paste("Max expression value:", max(gtex_counts_filtered)))
print(paste("Non-zero values:", sum(gtex_counts_filtered > 0)))

length(intersect(rownames(tcga_genes), rownames(final_counts)))
length(intersect(rownames(gtex_counts_filtered),rownames(gtex_genes_filtered)))
 
#---block ends---- 
# 1. Extract the sample metadata (colData) from your filtered object
gtex_meta <- as.data.frame(colData(gtex_final_rse))

# 2. Now you can create the 'stage_group' column
gtex_meta$stage_group <- "Normal"

# 3. Double-check it worked
head(gtex_meta)

#Assign GTEx Metadata
gtex_meta$stage_group <- "Normal"
gtex_meta$dataset     <- "GTEx"
rownames(gtex_meta) <- colnames(gtex_counts_filtered)

View(gtex_counts_filtered)
# 1. Find common genes between the FILTERED datasets
common_genes <- intersect(rownames(final_counts_common), 
                          rownames(gtex_counts_filtered))
View(common_genes)
View(tcga_genes)
length(common_genes)
# 2. Subset both to common genes
tcga_subset <- final_counts[common_genes, ]
gtex_subset <- gtex_counts_filtered[common_genes, ]

# 3. Final dimensions check# 16145
print(dim(tcga_subset))
print(dim(gtex_subset))


# Merge count matrices
combined_counts <- cbind(tcga_subset, gtex_subset)
print(dim(combined_counts))  


# Ensure metadata rownames match sample IDs
rownames(tcga_meta_clean) <- colnames(tcga_subset)
rownames(gtex_meta) <- colnames(gtex_subset)


# Use the subsetted versions that both have 16,420 rows
#merge counts
combined_counts <- cbind(tcga_subset, gtex_subset)

# Check dimensions to confirm
print(dim(combined_counts)) 
# Should be [16420, 471] (276 TCGA + 195 GTEx samples)

combined_meta <- rbind(
  tcga_meta[, c("stage_group","dataset")],
  gtex_meta[, c("stage_group","dataset")]
)
View(combined_meta)
#Ensure metadata order matches counts
combined_meta <- combined_meta[match(colnames(combined_counts),
                                     rownames(combined_meta)), ]

#check
all(rownames(combined_meta) == colnames(combined_counts))

#Convert factors

combined_meta$stage_group <- factor(
  combined_meta$stage_group,
  levels = c("Normal","Early","Late")
)
combined_meta$condition <- ifelse(
  combined_meta$stage_group == "Normal",
  "Normal",
  "Tumor"
)

View(combined_counts)
sum(rowSums(combined_counts) == 0)
combined_counts <- combined_counts[rowSums(combined_counts) > 0, ]

write.csv(
  combined_counts,
  file = "combined_counts.csv")

early_samples <- rownames(combined_meta[combined_meta$stage_group == "Early", ])
early_counts <- combined_counts[, early_samples]
early_meta <- combined_meta[early_samples, ]
View(early_counts)
write.csv(
 early_counts,
  file = "early_counts.csv")

late_samples <- rownames(combined_meta[combined_meta$stage_group == "Late", ])
late_counts <- combined_counts[, late_samples]
late_meta <- combined_meta[late_samples, ]
View(late_samples)
write.csv(
late_counts,
  file = "late_counts.csv")

#9.Filter Low-Expressed Genes
keep <- rowSums(combined_counts >= 10) >= 15
combined_counts <- combined_counts[keep, ]

vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "stage_group")
#10.Run DESeq2 (Batch Corrected)
View(combined_counts)
dds <- DESeqDataSetFromMatrix(
  countData = combined_counts,
  colData   = combined_meta,
  design    = ~ stage_group
)
dim(dds)
# Filter low count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

library(sva)
keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep,]
# ... (Previous code: combined_meta setup and filtering) ...

# --- NEW: START SVA CODE HERE ---
# A. Get normalized counts for SVA
dds <- estimateSizeFactors(dds)
dat <- counts(dds, normalized = TRUE)
dat <- dat[rowMeans(dat) > 1, ] # Extra noise filter for SVA

# B. Run SVA
mod  <- model.matrix(~ stage_group, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
svseq <- svaseq(dat, mod, mod0)

# Keep only the first 15 Surrogate Variables
n_sv_to_keep <- 10
sv_reduced <- svseq$sv[, 1:n_sv_to_keep]

# Rename them for clarity
sv_names <- paste0("SV", 1:ncol(sv_reduced))

# Add only these 15 to your colData
colData(dds) <- cbind(colData(dds), sv_reduced)

# Assign valid names to the new columns in colData
# (This ensures DESeq2 can "see" them in the design formula)
colnames(colData(dds))[ (ncol(colData(dds)) - ncol(sv_reduced) + 1) : ncol(colData(dds)) ] <- sv_names

# Update design formula
design(dds) <- as.formula(
  paste("~", paste(c(sv_names, "stage_group"), collapse=" + "))
)
View(design(dds))
# Run DESeq (This will now be much faster)
dds <- DESeq(dds, parallel = TRUE) # Remember to use parallel=TRUE if you set up the cluster!

table(colData(dds)$stage_group)

View(dds)
View(sv_reduced)
View(colData(dds))
# --- NEW: END SVA CODE ---

# 10. Run DESeq2 (Now Batch Corrected)


# ... (Continue with results extraction) ...

dds <- DESeq(dds)
svseq$n.sv
vsd <- vst(dds)

plotPCA(vsd, intgroup = "stage_group")

dim(dds)
####refer Deseq2.R code 
