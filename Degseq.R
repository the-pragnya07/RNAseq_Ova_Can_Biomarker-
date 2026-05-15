############################################################
# 1. LOAD LIBRARIES
############################################################

library(DESeq2)
library(sva)
library(BiocParallel)
library(EnhancedVolcano)
library(pheatmap)
library(dplyr)


############################################################
# 2. CREATE DESEQ DATASET
############################################################

dds <- DESeqDataSetFromMatrix(
  countData = combined_counts,
  colData   = combined_meta,
  design    = ~ stage_group
)

############################################################
# 3. FILTER LOW COUNT GENES
############################################################

keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep, ]



############################################################
# 7. RUN DESEQ2
############################################################

dds <- DESeq(dds)

dim(dds)

table(colData(dds)$stage_group)


############################################################
# 8. VARIANCE STABILIZATION + PCA
############################################################

vsd <- vst(dds)

plotPCA(vsd, intgroup = "stage_group")


############################################################
# 9. EARLY vs NORMAL
############################################################

res_early_vs_normal <- results(
  dds,
  contrast = c("stage_group","Early","Normal")
)
dim(res_early_vs_normal)
View(res_early_vs_normal)

res_early_vs_normal <- as.data.frame(res_early_vs_normal)
res_early_vs_normal$gene <- rownames(res_early_vs_normal)

write.csv(res_early_vs_normal,
          "Early_vs_Normal_all_genes.csv")


############################################################
# 10. LATE vs NORMAL
############################################################

res_late_vs_normal <- results(
  dds,
  contrast = c("stage_group","Late","Normal")
)

res_late_vs_normal <- as.data.frame(res_late_vs_normal)
res_late_vs_normal$gene <- rownames(res_late_vs_normal)

write.csv(res_late_vs_normal,
          "Late_vs_Normal_all_genes.csv")
nrow(res_late_vs_normal)

############################################################
# 12. LATE vs EARLY  ⭐ MOST IMPORTANT
############################################################

res_late_vs_early <- results(
  dds,
  contrast = c("stage_group","Late","Early")
)

res_late_vs_early <- as.data.frame(res_late_vs_early)
res_late_vs_early$gene <- rownames(res_late_vs_early)

write.csv(res_late_vs_early,
          "Late_vs_Early_all_genes.csv")


############################################################
# 12. EARLY vs LATE  ⭐ (Late = base/control)
############################################################

res_early_vs_late <- results(
  dds,
  contrast = c("stage_group","Early","Late")
)

res_early_vs_late <- as.data.frame(res_early_vs_late)
res_early_vs_late$gene <- rownames(res_early_vs_late)

write.csv(res_early_vs_late,
          "Early_vs_Late_all_genes.csv")






library(ggplot2)

res_early_vs_normal$Significance <- "Not Significant"

res_early_vs_normal$Significance[
  res_early_vs_normal$padj < 0.05 &
    res_early_vs_normal$log2FoldChange >= 2
] <- "Up"

res_early_vs_normal$Significance[
  res_early_vs_normal$padj < 0.05 &
    res_early_vs_normal$log2FoldChange <= -2
] <- "Down"

# Handle zero p-values (important)

# Handle zero p-values
res_early_vs_normal$padj[res_early_vs_normal$padj == 0] <- 1e-300

ggplot(res_early_vs_normal,
       aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  
  geom_point(alpha = 0.6, size = 1.5) +
  
  # ✅ Correct threshold (match your filtering)
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  
  scale_color_manual(values = c(
    "Up" = "red",
    "Down" = "blue",
    "Not Significant" = "grey"
  )) +
  
  # keep zoom, don’t cut data
  coord_cartesian(ylim = c(0, 150)) +
  
  scale_y_continuous(
    breaks = seq(0, 150, by = 50)
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Volcano Plot: Early vs Normal",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )
#late_vs_normal

res_late_vs_normal$Significance <- "Not Significant"

res_late_vs_normal$Significance[
  res_late_vs_normal$padj < 0.05 &
    res_late_vs_normal$log2FoldChange >= 2
] <- "Up"

res_late_vs_normal$Significance[
  res_late_vs_normal$padj < 0.05 &
    res_late_vs_normal$log2FoldChange <= -2
] <- "Down"

# Replace infinite values (important!)
res_late_vs_normal$padj[res_late_vs_normal$padj == 0] <- 1e-300

ggplot(res_late_vs_normal,
       aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  
  geom_point(alpha = 0.6, size = 1.5) +
  
  # Use correct thresholds (match your filtering)
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  scale_color_manual(values = c("Up" = "red", 
                                "Down" = "blue", 
                                "Not Significant" = "grey")) +
  # keep zoom, don’t cut data
  coord_cartesian(ylim = c(0, 150)) +
  
  scale_y_continuous(
    breaks = seq(0, 150, by = 50)
  ) +
  
  theme_minimal(base_size = 14) +
  
  labs(title = "Volcano Plot: Late vs Normal",
       x = "Log2 Fold Change",
       y = "-log10 Adjusted P-value")
#late_vs_early

res_late_vs_early$Significance <- "Not Significant"

res_late_vs_early$Significance[
  res_late_vs_early$padj < 0.05 &
    res_late_vs_early$log2FoldChange >= 2
] <- "Up in Late"

res_late_vs_early$Significance[
  res_late_vs_early$padj < 0.05 &
    res_late_vs_early$log2FoldChange <= -2
] <- "Up in Early"

ggplot(res_late_vs_early,
       aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("Up in Late" = "red",
                                "Up in Early" = "blue",
                                "Not Significant" = "grey")) +
  coord_cartesian(ylim = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 20, by = 50)) +
  theme_minimal() +
  labs(title = "Volcano Plot: Late vs Early",
       x = "Log2 Fold Change",
       y = "-log10 Adjusted P-value")


library(ggplot2)

# Add significance column
res_early_vs_late$Significance <- "Not Significant"

res_early_vs_late$Significance[
  res_early_vs_late$padj < 0.05 &
    res_early_vs_late$log2FoldChange >= 2
] <- "Up"

res_early_vs_late$Significance[
  res_early_vs_late$padj < 0.05 &
    res_early_vs_late$log2FoldChange <= -2
] <- "Down"
res_early_vs_late$padj[res_early_vs_late$padj == 0] <- 1e-300
# Replace 0 or extremely small values with minimum threshold

res_early_vs_normal$padj <- pmax(res_early_vs_normal$padj, 1e-300)
res_late_vs_normal$padj  <- pmax(res_late_vs_normal$padj, 1e-300)
# Volcano plot
# Handle zero p-values


# Handle zero p-values
res_early_vs_late$padj[res_early_vs_late$padj == 0] <- 1e-300

ggplot(res_early_vs_late,
       aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  
  geom_point(alpha = 0.6, size = 1.5) +
  
  # ✅ Correct threshold
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  scale_color_manual(values = c(
    "Up" = "red",
    "Down" = "blue",
    "Not Significant" = "grey"
  )) +
  
  # keep zoom, don't cut data
  coord_cartesian(ylim = c(0, 150)) +
  
  scale_y_continuous(
    breaks = seq(0, 150, by = 50)
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Volcano Plot: Early vs Late",
    x = "Log2 Fold Change (Early / Late)",
    y = "-log10 Adjusted P-value"
  )


#-------------------------------
#deg filtering for all results genes
res_early_vs_normal <- res_early_vs_normal[!is.na(res_early_vs_normal$padj), ]
res_late_vs_normal  <- res_late_vs_normal[!is.na(res_late_vs_normal$padj), ]
res_late_vs_early   <- res_late_vs_early[!is.na(res_late_vs_early$padj), ]
res_early_vs_late <- res_early_vs_late[!is.na(res_early_vs_late$padj), ]

############################################################
# 11. DEG FILTERING
############################################################

deg_early <- res_early_vs_normal %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 2)
nrow(deg_early)
deg_late <- res_late_vs_normal %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 2)
nrow(deg_late)
deg_late_vs_early <- res_late_vs_early %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 2)

deg_early_vs_late <- res_early_vs_late %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 2)

nrow(deg_late_vs_early)
nrow(deg_early_vs_late)
nrow(deg_early)
nrow(deg_late)


############################################################
# 12. UP & DOWN REGULATED GENES
############################################################

up_early <- deg_early %>%
  filter(log2FoldChange >= 2) %>%
  arrange(desc(log2FoldChange))

down_early <- deg_early %>%
  filter(log2FoldChange <= -2) %>%
  arrange(log2FoldChange)

up_reg <- deg_late_vs_early %>%
  filter(log2FoldChange >= 2)

down_reg <- deg_late_vs_early %>%
  filter(log2FoldChange < -1)

up_late <- deg_late %>%
  filter(log2FoldChange >= 2) %>%
  arrange(desc(log2FoldChange))

down_late <- deg_late %>%
  filter(log2FoldChange <= -2) %>%
  arrange(log2FoldChange)

up_reg_early <- deg_early_vs_late %>%
  filter(log2FoldChange >= 2)

down_reg_late <- deg_early_vs_late %>%
  filter(log2FoldChange <= -2)

nrow(up_reg_early)
nrow(down_reg_late)
nrow(up_early)
nrow(down_early)
nrow(up_late)
nrow(down_late)
nrow(up_reg)

# Early vs Normal
write.csv(up_early, "Upregulated_Early_vs_Normal.csv", row.names = FALSE)
write.csv(down_early, "Downregulated_Early_vs_Normal.csv", row.names = FALSE)

# Late vs Normal
write.csv(up_late, "Upregulated_Late_vs_Normal.csv", row.names = FALSE)
write.csv(down_late, "Downregulated_Late_vs_Normal.csv", row.names = FALSE)

# Late vs Early (old direction)
write.csv(up_reg, "Upregulated_Late_vs_Early.csv", row.names = FALSE)
write.csv(down_reg, "Downregulated_Late_vs_Early.csv", row.names = FALSE)

# Early vs Late (correct direction for your goal)
write.csv(up_reg_early, "Upregulated_Early_vs_Late.csv", row.names = FALSE)
write.csv(down_reg_late, "Downregulated_Early_vs_Late.csv", row.names = FALSE)

View(up_reg_early)


common_all <- Reduce(intersect, list(
  genes_up_early,
  genes_up_late,
  genes_up_early_vs_late
))

length(common_all)

early_only <- setdiff(
  genes_up_early,
  union(genes_up_late, genes_up_early_vs_late)
)

length(early_only)

genes_up_early <- rownames(up_early)
genes_up_late  <- rownames(up_late)
genes_up_early_vs_late <- up_reg_early$gene

common_genes <- intersect(genes_up_early, genes_up_late)

diff_genes <- setdiff(
  setdiff(genes_up_early, genes_up_late),
  early_only
)

length(diff_genes)  # should be 4
diff_genes


early_only <- setdiff(genes_up_early, genes_up_late)

late_only <- setdiff(genes_up_late, genes_up_early)


length(common_genes) #early and late
length(early_only)
length(late_only)
length(reg_early_only)
nrow(up_early)

#Genes uniquely upregulated in Early vs Late
reg_early_only <- setdiff(
  genes_up_early_vs_late,
  union(genes_up_early, genes_up_late)
)

nrow(reg_early_only)
dim(reg_early_only)

#only in early late 
early_late <- setdiff(
  intersect(genes_up_early, genes_up_late),#it has genes which are presengt in both late and early
  genes_up_reg
)
length(early_late)
length(reg_only)

#comparision within
# Late ∩ Progression ONLY
late_prog <- setdiff(
  intersect(genes_up_late, genes_up_reg),
  genes_up_early
)
length(late_prog )
# Early ∩ Progression ONLY 
early_prog <- (
  intersect(genes_up_early, genes_up_reg))
length(early_prog) 
length(early_late) 

length(common_genes)
length(early_only) 
length(late_only) 
dim(up_early)


#venn diagram for detailes within comparision
genes_up_early <- rownames(up_early)
genes_up_late  <- rownames(up_late)

genes_up_early_vs_late <-rownames(up_reg_early)
View(genes_up_early)
library(VennDiagram)
library(grid)

venn.plot <- venn.diagram(
  x = list(
    UpEarly_vs_Normal  = genes_up_early,
    UpLate_vs_Normal = genes_up_late,
    UpEarly_vs_late= genes_up_early_vs_late
  ),
  filename = "Upregulated_Venn.png",
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  margin = 0.1
)

grid.draw(venn.plot)


# All three
common_all_three <- Reduce(intersect, list(
  genes_up_early, genes_up_late, genes_up_early_vs_late
))
length(common_all_three)
View(common_all_three)
early_deg <- up_early[rownames(up_early) %in% early_only, ]
late_deg <- up_late[rownames(up_late) %in% late_only, ]
prog_deg <- up_reg[rownames(up_reg) %in% reg_only, ]
early_vs_late_deg <-up_reg_early [rownames(up_reg_early) %in% reg_early_only,  ]
View(early_deg)
#save
write.csv(common_all_three, "common_all_three_up.csv", row.names = FALSE)
write.csv(common_early_late, "common_early_late_up.csv", row.names = FALSE)

write.csv(early_only, "early_only_up.csv", row.names = FALSE)
write.csv(late_only, "late_only_up.csv", row.names = FALSE)
write.csv(early_late, "early_late.csv", row.names = FALSE)
write.csv(prog_deg, "progression_only_DEGs.csv")
write.csv(early_deg, "early_only_DEGs.csv")
write.csv(late_deg, "late_only_DEGs.csv")
write.csv(early_vs_late_deg , "early_vs_late_only_DEGs.csv")

common_combined <- data.frame(
  gene = common_all,
  logFC_early = up_early[common_all, "log2FoldChange"],
  pval_early  = up_early[common_all, "pvalue"],
  
  logFC_late  = up_late[common_all, "log2FoldChange"],
  pval_late   = up_late[common_all, "pvalue"],
  
  logFC_reg   =up_reg_early[common_all, "log2FoldChange"],
  pval_reg    = up_reg_early[common_all, "pvalue"]
)
length(common_combined)
nrow(common_combined)
write.csv(common_combined, "common_50_DEGs.csv")


library(biomart)


#if not wor
# Clean IDs
up_early$ENSEMBL <- sub("\\..*", "", up_early$gene)

# Map to symbols
up_early$SYMBOL <- mapIds(
  org.Hs.eg.db,
  keys = up_early$ENSEMBL,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
up_early <- up_early[!is.na(up_early$SYMBOL), ]

write.csv(up_early, "up_early_with_symbols.csv")

#early_vs_late
# Clean Ensembl IDs
up_reg_early$ENSEMBL <- sub("\\..*", "", up_reg_early$gene)

# Map to gene symbols
up_reg_early$SYMBOL <- mapIds(
  org.Hs.eg.db,
  keys = up_reg_early$ENSEMBL,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Optional: remove NA
up_reg_early <- up_reg_early[!is.na(up_reg_early$SYMBOL), ]

write.csv(up_reg_early, "up_reg_early_with_symbols.csv")


#late vs normal
# Clean Ensembl IDs
up_late$ENSEMBL <- sub("\\..*", "", up_late$gene)

# Map to gene symbols
up_late$SYMBOL <- mapIds(
  org.Hs.eg.db,
  keys = up_late$ENSEMBL,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Optional: remove NA
up_late <- up_late[!is.na(up_late$SYMBOL), ]

write.csv(up_late, "up_late_with_symbols.csv")


#----
late_vs_early_map_df<- data.frame(
  ENSEMBL = names(prog_map),
  SYMBOL  = prog_map
)
write.csv(late_vs_early_map_df, "late_vs_early_map_symbols.csv", row.names = FALSE)

gene_map_early <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = early_only,
  mart = mart
)
gene_map_early <- gene_map_early[gene_map_early$hgnc_symbol != "", ]
write.csv(gene_map_early,
          "gene_map_early_HGNC.csv",
          row.names = FALSE)
#if its not working then use this
early_map <- mapIds(org.Hs.eg.db,
                    keys = early_only,
                    column = "SYMBOL",
                    keytype = "ENSEMBL",
                    multiVals = "first")

early_df <- data.frame(ENSEMBL = names(early_map), SYMBOL = early_map)
write.csv(early_df, "early_df_names.csv")

gene_map_late <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = late_only,
  mart = mart
)
gene_map_late <- gene_map_late[gene_map_late$hgnc_symbol != "", ]
write.csv(gene_map_late,
          "gene_map_late_HGNC.csv",
          row.names = FALSE)
#if not working above run thsi

late_map <- mapIds(org.Hs.eg.db,
                   keys = late_only,
                   column = "SYMBOL",
                   keytype = "ENSEMBL",
                   multiVals = "first")

late_df <- data.frame(ENSEMBL = names(late_map), SYMBOL = late_map)
write.csv(late_df, "late_df_names.csv")
#top 50
gene_map_common <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = common_all,
  mart = mart
)
#if not working then
common_map <- mapIds(org.Hs.eg.db,
                     keys = common_all,
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")

common_df <- data.frame(ENSEMBL = names(common_map), SYMBOL = common_map)
View(common_df)
write.csv(common_df, "common_df_names.csv")
# 13. TOP GENES
############################################################

top_up_early   <- head(up_early,50)   #look here 
top_down_early <- head(down_early,50)

top_up_late   <- head(up_late,50)
top_down_late <- head(down_late,50)

top_up_reg <-head(up_reg,100)
top_down_reg <-head(down_reg,100)
############################################################
# 14. SAVE RESULTS
############################################################

write.csv(top_up_early,
          "Top50_Up_Early_vs_Normal.csv") #extarct top hundered

write.csv(top_down_early,
          "Top50_Down_Early_vs_Normal.csv")

write.csv(top_up_late,
          "Top50_Up_Late_vs_Normal.csv")

write.csv(top_down_late,
          "Top50_Down_Late_vs_Normal.csv")

#tumor
############################################################
# ADD: TUMOR vs NORMAL VARIABLE
############################################################

#if no sva added
combined_meta$tumor_status <- ifelse(
  combined_meta$stage_group == "Normal",
  "Normal",
  "Tumor"
)

combined_meta$tumor_status <- factor(combined_meta$tumor_status)

#dds_run
dds_tumor <- DESeqDataSetFromMatrix(
  countData = combined_counts,
  colData   = combined_meta,
  design    = ~ tumor_status
)

dds_tumor <- DESeq(dds_tumor)

res_tumor <- results(
  dds_tumor,
  contrast = c("tumor_status","Tumor","Normal")
)

res_tumor <- as.data.frame(res_tumor)
res_tumor$gene <- rownames(res_tumor)

deg_tumor <- res_tumor %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.05 & abs(log2FoldChange) >=2 )

nrow(deg_tumor)
nrow(res_tumor)
up_tumor <- deg_tumor %>%
  filter(log2FoldChange >=2 ) %>%
  arrange(desc(log2FoldChange))

down_tumor <- deg_tumor %>%
  filter(log2FoldChange <= -2) %>%
  arrange(log2FoldChange)

dim(up_tumor)
dim(down_tumor)
############################################################
#add gene names

top_up_tumor <- head(up_tumor, 100)
top_down_tumor <- head(down_tumor, 100)


top_up_tumor$ensembl_gene_id <- gsub("\\..*", "", rownames(top_up_tumor))
top_down_tumor$ensembl_gene_id <- gsub("\\..*", "", rownames(top_down_tumor))


library(biomaRt)

mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "asia"   # try "useast" if this fails
)
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = unique(c(top_up_tumor$ensembl_gene_id,
                    top_down_tumor$ensembl_gene_id)),
  mart = mart
)
library(dplyr)

# Sort by adjusted p-value (best practice)
up_tumor_sorted <- up_tumor %>% arrange(padj)
down_tumor_sorted <- down_tumor %>% arrange(padj)

# Take top 100
top_up_tumor <- head(up_tumor_sorted, 100)
top_down_tumor <- head(down_tumor_sorted, 100)


top_up_tumor$ensembl_gene_id <- gsub("\\..*", "", rownames(top_up_tumor))
top_down_tumor$ensembl_gene_id <- gsub("\\..*", "", rownames(top_down_tumor))


library(biomaRt)

mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "useast"   # more stable than asia
)
#merge
library(dplyr)

top_up_tumor_named <- top_up_tumor %>%
  left_join(gene_map, by = "ensembl_gene_id")

top_down_tumor_named <- top_down_tumor %>%
  left_join(gene_map, by = "ensembl_gene_id")

sum(top_up_tumor_named$hgnc_symbol == "")
write.csv(top_up_tumor_named, "Top100_up_tumor_names.csv", row.names = FALSE)
write.csv(top_down_tumor_named, "Top100_down_tumor_names.csv", row.names = FALSE)

############################################################
# 18. EXPORT RESULTS
############################################################

write.csv(res_tumor,
          "Tumor_vs_Normal_DESeq2_results.csv")

write.csv(sig,
          "Tumor_vs_Normal_significant_DEGs.csv")


############################################################
# END PIPELINE
############################################################

common_genes <- intersect(up_early, up_late)
length(common_genes)
early_only <- setdiff(up_early, up_late)
late_only  <- setdiff(up_late,up_early)

length(common_genes)
length(early_only)
length(late_only)
############################################################
# 3-WAY VENN DIAGRAM
############################################################

library(VennDiagram)
library(grid)

# Apply significance filtering

early_genes <- gsub("\\..*", "", rownames(deg_early))
late_genes  <- gsub("\\..*", "", rownames(deg_late))
tumor_genes <- gsub("\\..*", "", deg_tumor$gene_id)

venn.plot <- venn.diagram(
  x = list(
    Early = early_genes,
    Late  = late_genes,
    Tumor = tumor_genes
  ),
  filename = NULL,
  fill = c("red","blue","green"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5
)

grid.draw(venn.plot)
#extract the u ique gene names
# 1. Genes UNIQUE to 'Early' (the 556)
unique_early <- setdiff(early_genes, union(late_genes, tumor_genes))

# 2. Genes UNIQUE to 'Late' (the 140)
unique_late <- setdiff(late_genes, union(early_genes, tumor_genes))

# 3. Genes UNIQUE to 'Tumor' (the 150)
unique_tumor <- setdiff(tumor_genes, union(early_genes, late_genes))

# To see the names of the unique early genes:
print(unique_early)
length(unique_early)

#Ma plot
# Early vs Normal
plotMA(deg_early, main="MA Plot: Early vs Normal", ylim=c(-5,5))

# Late vs Normal
plotMA(deg_late, main="MA Plot: Late vs Normal", ylim=c(-5,5))

# Tumor vs Normal
plotMA(deg_tumor, main="MA Plot: Tumor vs Normal", ylim=c(-5,5))

library(EnhancedVolcano)
#early valcano

deg_early_df <- as.data.frame(deg_early)

# remove version numbers (VERY IMPORTANT)
deg_early_df$ensembl_gene_id <- gsub("\\..*", "", rownames(deg_early_df))
#connect to ensembel id"s
library(biomaRt)
View(deg_early)
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "asia"
)
#get gene names 
gene_early <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = deg_early_df$ensembl_gene_id,   # ✅ CORRECT
  mart = mart
)
#merge properly
deg_early_named <- merge(
  deg_early_df,
  gene_early,
  by = "ensembl_gene_id",
  all.x = TRUE
)
#late valcano
EnhancedVolcano(
  deg_late,
  lab = rownames(deg_late),
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Volcano Plot: Late vs Normal',
  pCutoff = 0.05,
  FCcutoff = 1
)
#





############################################################
# common_genes <- intersect(rownames(deg_early), rownames(deg_late))

length(common_genes)


dds <- DESeq(dds)
res_all_tumor_vs_normal <- results(
  dds,
  contrast = c("tumor_status","Tumor","Normal")
)

res_all_tumor_vs_normal <- as.data.frame(res_all_tumor_vs_normal)

res_all_tumor_vs_normal$gene <- rownames(res_all_tumor_vs_normal)



############################################################
# 16. SIGNIFICANT DEGs
############################################################

sig <- res_all_tumor_vs_normal[
  res_all_tumor_vs_normal$padj < 0.05 &
    abs(res_all_tumor_vs_normal$log2FoldChange) > 1,
]

nrow(sig)


############################################################
# 17. HIGH CONFIDENCE DEGs
############################################################

res_highconf <- res_all_tumor_vs_normal %>%
  filter(padj < 0.001 & abs(log2FoldChange) > 2.5)

highconf_genes <- res_highconf$gene

length(highconf_genes)
head(highconf_genes)


############################################################
# 18. EXPORT RESULTS
############################################################

write.csv(res_all_tumor_vs_normal,
          "Tumor_vs_Normal_DESeq2_results.csv")

write.csv(sig,
          "Tumor_vs_Normal_significant_DEGs.csv")

write.csv(res_highconf,
          "Tumor_vs_Normal_highconf_DEGs.csv")

############################################################
# END PIPELINE
############################################################