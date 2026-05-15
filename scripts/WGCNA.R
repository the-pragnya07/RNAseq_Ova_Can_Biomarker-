# 1. Extract and Transpose
datExpr <- t(assay(vsd))

# Check if sample names match exactly
all(rownames(datExpr) == rownames(combined_meta))


# 2. Final check for "Bad" Genes/Samples
library(WGCNA)
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}
# Create a comparison table of the first 10 samples
check_alignment <- data.frame(
  Expression_Samples = rownames(datExpr)[1:10],
  Metadata_Samples = rownames(combined_meta)[1:10]
)
View(check_alignment)

# View it
print(check_alignment)
#don"t run now wait and go to next block
# 1. Create a copy of your metadata for traits
datTraits <- combined_meta

# 2. Convert 'stage_group' to a numeric factor
# Normal = 1, Early = 2, Late = 3 (or however you want to order them)
datTraits$stage_numeric <- as.numeric(factor(datTraits$stage_group, 
                                             levels = c("Normal", "Early", "Late")))

# 3. (Optional) Convert 'dataset' to numeric to check for batch effects later
# GTEx = 1, TCGA = 2
datTraits$dataset_numeric <- as.numeric(factor(datTraits$dataset))

# 4. Keep only the numeric columns for WGCNA trait correlation
datTraits <- datTraits[, c("stage_numeric", "dataset_numeric")]


#Step B: Pick the Soft Threshold (The Power \(\beta \))
# 3. Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# 4. Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, 
                         powerVector = powers, 
                         verbose = 5)

# 5. Plot the results to pick your power
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2",
     type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red")
abline(h=0.80, col="red") # Aim for a power that crosses this line

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", 
     type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

# Define your power based on the results above
softPower = 2
####run the other block
# Run the network construction
net = blockwiseModules(
  datExpr, 
  power = softPower,
  TOMType = "unsigned",      # "unsigned" is standard for combined datasets
  minModuleSize = 30,        # Minimum number of genes in a module
  remergeCutHeight = 0.25,   # Merges modules that are very similar
  numericLabels = TRUE,      # Uses numbers for modules (we'll change to colors later)
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,           # Save the Topological Overlap Matrix
  saveTOMFileBase = "OvarianCancerTOM",
  verbose = 3
)
#run this not up one
net = blockwiseModules(
  datExpr,
  softpower = 2,
  TOMType = "unsigned",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  verbose = 3
)

moduleColors <- labels2colors(net$colors)

MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs)
# See how many modules were found
table(net$colors)
signedKME(datExpr, MEs)

plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)

table(moduleColors)
# 1. Convert numeric labels to colors for visualization
moduleColors = labels2colors(net$colors)

# 2. See how many modules you have and their sizes
table(moduleColors)

# 3. Plot the Cluster Dendrogram (The famous WGCNA "tree" plot)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guide
                  Hang = 0.05)


traitData <- data.frame(
  Early = ifelse(combined_meta$stage_group == "Early", 1, 0),
  Late  = ifelse(combined_meta$stage_group == "Late", 1, 0)
)

moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))


textMatrix <- paste(
  signif(moduleTraitCor, 2), "\n(",
  signif(moduleTraitPvalue, 1), ")",
  sep = ""
)

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(traitData),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.7,
  zlim = c(-1,1),
  main = "Module-trait relationships"
)


#Convert ENSEMBL → SYMBOL

library(org.Hs.eg.db)

module_symbol_list = list()

for (m in names(module_gene_list)) {
  
  genes = module_gene_list[[m]]
  
  symbols = mapIds(org.Hs.eg.db,
                   keys = genes,
                   column = "SYMBOL",
                   keytype = "ENSEMBL",
                   multiVals = "first")
  
  module_symbol_list[[m]] = na.omit(symbols)
}

for (m in names(module_symbol_list)) {
  
  df = data.frame(Gene_Symbol = module_symbol_list[[m]])
  
  write.csv(df,
            paste0(m, "_genes_symbols.csv"),
            row.names = FALSE)
}

# 2. Extract genes for ME3, ME2, and ME16 (Top positive correlations)
module3_genes = genes[net$colors == 3]
module2_genes = genes[net$colors == 2]
module16_genes = genes[net$colors == 16]

# 3. Check the size (just to be sure)
length(module3_genes) # Should match your table (1140)

#to find hub genes inmodule 3 
kME = signedKME(datExpr, MEs)
#main point only significant modules are selected 
#Extract kME for ME3
kME_ME3 = kME[, "kME3"]

#Create dataframe
hub_df = data.frame(
  gene = colnames(datExpr),
  module = net$colors,
  kME = kME_ME3
)
#Filter only ME3 genes
hub_ME3 = hub_df[hub_df$module == 3, ]
#Sort by importance
hub_ME3_sorted = hub_ME3[order(-hub_ME3$kME), ]
#Get top hub genes
top_hubs = head(hub_ME3_sorted, 10)
print(top_hubs)
res_sorted
# Correlation of each gene with stage
GS = cor(datExpr, stage_numeric, use = "p")

GS = cor(datExpr, stage_numeric, use = "pairwise.complete.obs")

GS = as.vector(GS)
names(GS) = colnames(datExpr)
# Combine with hub data
hub_ME3$GS = GS[hub_ME3$gene]
hub_ME3$GS = GS[hub_ME3$gene]
sum(!is.na(hub_ME3$GS))
# Sort by both kME and GS

hub_ME3_clean = hub_ME3_clean[order(-abs(hub_ME3_clean$kME), -abs(hub_ME3_clean$GS)), ]

top_hubs_final = head(hub_ME3_clean, 10)

head(top_hubs_final)
common_ME3 = intersect(module3_genes, deg_genes)
length(common_ME3)

#plot
png("ME3_GS_vs_kME.png", width = 800, height = 800)

plot(abs(hub_ME3_clean$GS), abs(hub_ME3_clean$kME),
     xlab = "Gene Significance (GS)",
     ylab = "Module Membership (kME)",
     main = "ME3: GS vs kME",
     pch = 19)

dev.off()
library(clusterProfiler)
library(org.Hs.eg.db)

# Run GO enrichment for Biological Process (BP)
ego_BP <- enrichGO(gene          = module3_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENSEMBL",
                   ont           = "BP", 
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   readable      = TRUE)

# View the top results
head(ego_BP)


#early vs normal
early_vs_normal_trait = ifelse(
  combined_meta$stage_group %in% c("Early", "Normal"),
  ifelse(combined_meta$stage_group == "Early", 1, 0),
  NA
)

# Subset data
keep = !is.na(early_vs_normal_trait)

datExpr_sub = datExpr[keep, ]
trait_sub = early_vs_normal_trait[keep]
MEs_sub = MEs[keep, ]
#corelate
cor_EN = cor(MEs_sub, trait_sub, use = "p")
p_EN = corPvalueStudent(cor_EN, nrow(datExpr_sub))

res_EN = data.frame(
  Module = names(MEs),
  Correlation = cor_EN,
  Pvalue = p_EN
)

res_EN = res_EN[order(-res_EN$Correlation), ]
print(res_EN)


early_modules = res_EN$Module[res_EN$Correlation > 0.7]

early_modules

#extract
extract_hubs = function(module_name) {
  
  cat("\n========================\n")
  cat("Module:", module_name, "\n")
  
  module_number = as.numeric(sub("ME", "", module_name))
  
  genes = colnames(datExpr)
  module_genes = genes[net$colors == module_number]
  
  cat("Total genes:", length(module_genes), "\n")
  
  # kME
  kME = signedKME(datExpr, MEs)
  kME_module = kME[, paste0("kME", module_number)]
  
  hub_df = data.frame(
    gene = colnames(datExpr),
    module = net$colors,
    kME = kME_module
  )
  
  hub_module = hub_df[hub_df$module == module_number, ]
  
  # GS
  GS = cor(datExpr_sub, trait_sub, use = "pairwise.complete.obs")
  GS = as.vector(GS)
  names(GS) = colnames(datExpr_sub)
  
  hub_module$GS = GS[hub_module$gene]
  hub_module = hub_module[!is.na(hub_module$GS), ]
  
  # Sort
  hub_module = hub_module[
    order(-abs(hub_module$kME), -abs(hub_module$GS)), 
  ]
  
  # Top hubs
  top_hubs = head(hub_module, 20)
  
  # Key genes
  key_genes = hub_module[
    abs(hub_module$kME) > 0.9 &
      abs(hub_module$GS) > 0.6,
  ]
  
  
  # 🔥 Hidden regulators (FIXED)
  hidden_regulators = setdiff(key_genes$gene, deg_genes)
  
  cat("Top hubs:", nrow(top_hubs), "\n")
  cat("Key genes:", nrow(key_genes), "\n")
  cat("Hidden regulators:", length(hidden_regulators), "\n")
  
  return(list(
    module = module_name,
    module_genes = module_genes,
    top_hubs = top_hubs,
    key_genes = key_genes,
    hidden_regulators = hidden_regulators
  ))
}
deg_genes = rownames(res_early_vs_normal[
   !is.na(res_early_vs_normal$padj) &
         res_early_vs_normal$padj < 0.05 &
         abs(res_early_vs_normal$log2FoldChange) > 2,
])

moduleColors   # vector of module assignments for each gene
datExpr        # expression matrix (rownames = samples, colnames = genes)
moduleTraitCor
#run early modules
early_results = lapply(early_modules, extract_hubs)

#early biomarkers
early_hub_genes = unlist(lapply(early_results, function(x) x$key_genes$gene))

early_biomarkers = intersect(early_hub_genes, up_early)

length(early_biomarkers)
head(early_biomarkers)

library(clusterProfiler)
library(org.Hs.eg.db)

get_module_kegg = function(module_name) {
  
  cat("\n========================\n")
  cat("Module:", module_name, "\n")
  
  module_number = as.numeric(sub("ME", "", module_name))
  
  # Extract genes
  module_genes = colnames(datExpr)[net$colors == module_number]
  
  # Remove version numbers
  module_genes = sub("\\..*", "", module_genes)
  
  # Convert IDs
  gene_df = bitr(
    module_genes,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  # Remove duplicates
  gene_df = gene_df[!duplicated(gene_df$ENTREZID), ]
  
  # KEGG
  kegg = enrichKEGG(
    gene = gene_df$ENTREZID,
    organism = "hsa",
    pvalueCutoff = 0.05
  )
  
  # If no pathways
  if (is.null(kegg) || nrow(kegg@result) == 0) {
    cat("No significant KEGG pathways\n")
    return(NULL)
  }
  
  # Print top pathways
  cat("Top pathways:\n")
  print(head(kegg@result[, c("Description", "p.adjust")], 5))
  
  return(kegg)
}
#for early
early_modules = c("ME4", "ME14", "ME9", "ME11", "ME2")

kegg_results = lapply(early_modules, get_module_kegg)

names(kegg_results) = early_modules

#create clean module pathway
module_summary = data.frame()

for (m in names(kegg_results)) {
  
  kegg = kegg_results[[m]]
  
  if (!is.null(kegg)) {
    
    top_pathways = head(kegg@result$Description, 3)
    
    temp = data.frame(
      Module = m,
      Pathway = top_pathways
    )
    
    module_summary = rbind(module_summary, temp)
  }
}

module_summary

library(org.Hs.eg.db)

get_module_symbols = function(module_name) {
  
  module_number = as.numeric(sub("ME", "", module_name))
  
  module_genes = colnames(datExpr)[net$colors == module_number]
  
  # Remove version numbers
  module_genes = sub("\\..*", "", module_genes)
  
  # Convert to SYMBOL
  symbols = mapIds(
    org.Hs.eg.db,
    keys = module_genes,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  symbols = na.omit(symbols)
  
  return(symbols)
}#get all modules
module_symbol_list = lapply(early_modules, get_module_symbols)
names(module_symbol_list) = early_modules

#save
for (m in names(module_symbol_list)) {
  
  df = data.frame(Gene = module_symbol_list[[m]])
  
  write.csv(df, paste0(m, "_symbols.csv"), row.names = FALSE)
}



#early vs late
EL_trait = ifelse(
  combined_meta$stage_group %in% c("Early", "Late"),
  ifelse(combined_meta$stage_group == "Late", 1, 0),
  NA
)

keep = !is.na(EL_trait)

datExpr_EL = datExpr[keep, ]
trait_EL = EL_trait[keep]
MEs_EL = MEs[keep, ]
#COREALTE
cor_EL = cor(MEs_EL, trait_EL, use = "p")
p_EL = corPvalueStudent(cor_EL, nrow(datExpr_EL))

res_EL = data.frame(
  Module = names(MEs),
  Correlation = cor_EL,
  Pvalue = p_EL
)

res_EL = res_EL[order(-abs(res_EL$Correlation)), ]
print(res_EL)
#SELECT
late_modules  = res_EL$Module[res_EL$Correlation > 0.3]
early_modules = res_EL$Module[res_EL$Correlation < -0.3]

late_modules
early_modules

#LATE VS NORMAL
trait_LN = ifelse(combined_meta$stage_group == "Late", 1,
                  ifelse(combined_meta$stage_group == "Normal", 0, NA))

keep = !is.na(trait_LN)

datExpr_LN = datExpr[keep, ]
MEs_LN = MEs[keep, ]
trait_LN = trait_LN[keep]
#COREALTE
#corelate
cor_LN = cor(MEs_LN, trait_LN, use = "p")
p_LN = corPvalueStudent(cor_LN, nrow(datExpr_LN))

res_LN = data.frame(
  Module = names(MEs),
  Correlation = cor_LN,
  Pvalue = p_LN
)

res_LN = res_LN[order(-res_LN$Correlation), ]
print(res_LN)

late_modules = res_LN$Module[res_LN$Correlation > 0.7]
late_modules

#extract
extract_hubs_late = function(module_name) {
  
  module_number = as.numeric(sub("ME", "", module_name))
  
  genes = colnames(datExpr)
  
  module_genes = genes[net$colors == module_number]
  
  kME = signedKME(datExpr, MEs)
  kME_module = kME[, paste0("kME", module_number)]
  
  hub_df = data.frame(
    gene = colnames(datExpr),
    module = net$colors,
    kME = kME_module
  )
  
  hub_module = hub_df[hub_df$module == module_number, ]
  
  # GS for late
  GS = cor(datExpr_sub_L, trait_sub_L, use = "pairwise.complete.obs")
  GS = as.vector(GS)
  names(GS) = colnames(datExpr_sub_L)
  
  hub_module$GS = GS[hub_module$gene]
  
  hub_module = hub_module[!is.na(hub_module$GS), ]
  
  hub_module = hub_module[
    order(-abs(hub_module$kME), -abs(hub_module$GS)), 
  ]
  
  key_genes = hub_module[
    abs(hub_module$kME) > 0.9 &
      abs(hub_module$GS) > 0.6,
  ]
  
  return(key_genes)
}
late_results = lapply(late_modules, extract_hubs_late)

#late biomarket
late_hub_genes = unlist(lapply(late_results, function(x) x$gene))

late_biomarkers = intersect(late_hub_genes, up_late)

length(late_biomarkers)
head(late_biomarkers)

#to extarct the genes from modules and compare with DEG 
# Get numeric labels (IMPORTANT)
moduleLabels = net$colors   # numeric labels (0,1,2,...)

# Convert to colors
moduleColors_correct = labels2colors(moduleLabels)

# Check mapping
table(moduleLabels, moduleColors_correct)

all_genes = colnames(datExpr)

# Use numeric labels (MOST RELIABLE)
ME4_genes = all_genes[moduleLabels == 4]
ME2_genes = all_genes[moduleLabels == 2]

length(ME4_genes)
length(ME2_genes)
library(clusterProfiler)
library(org.Hs.eg.db)

 ME4_df = bitr(ME4_genes,
                     fromType = "ENSEMBL",
                     toType = c("SYMBOL"),
                  OrgDb = org.Hs.eg.db)

 ME2_df = bitr(ME2_genes,
                       fromType = "ENSEMBL",
                      toType = c("SYMBOL"),
                      OrgDb = org.Hs.eg.db) 
 
 ME4_deg = intersect(ME4_genes, deg_genes)
  ME2_deg = intersect(ME2_genes, deg_genes)
 length(ME4_deg)

length(ME2_deg)
