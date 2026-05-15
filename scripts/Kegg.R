library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)


top300_up_early <- head(up_early, 300)
top300_up_late <- head(up_late, 300)
top23_up_reg_early <- head(up_reg_early, 23)

#Early vs Normal
top300_up_early$ENSEMBL <- sub("\\..*", "", top300_up_early$gene)

early_entrez <- mapIds(
  org.Hs.eg.db,
  keys = top300_up_early$ENSEMBL,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

early_entrez <- na.omit(early_entrez)

#🔹 Late vs Normal
top300_up_late$ENSEMBL <- sub("\\..*", "", top300_up_late$gene)

late_entrez <- mapIds(
  org.Hs.eg.db,
  keys = top300_up_late$ENSEMBL,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

late_entrez <- na.omit(late_entrez)

#🔹 Early vs Late
top23_up_reg_early$ENSEMBL <- sub("\\..*", "", top23_up_reg_early$gene)

early_vs_late_entrez <- mapIds(
  org.Hs.eg.db,
  keys = top23_up_reg_early$ENSEMBL,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

early_vs_late_entrez <- na.omit(early_vs_late_entrez)

#✅ STEP 3: KEGG analysis
kegg_early <- enrichKEGG(gene = early_entrez, organism = "hsa", pvalueCutoff = 0.05)
kegg_late <- enrichKEGG(gene = late_entrez, organism = "hsa", pvalueCutoff = 0.05)
kegg_early_vs_late <- enrichKEGG(gene = early_vs_late_entrez, organism = "hsa", pvalueCutoff = 0.05)
kegg_early_vs_late <- enrichKEGG(
  gene = early_vs_late_entrez,
  organism = "hsa",
  pvalueCutoff = 0.2
)
#STEP 4: Visualization
dotplot(kegg_early, showCategory = 10, title = "Early vs Normal")
dotplot(kegg_late, showCategory = 10, title = "Late vs Normal")
dotplot(kegg_early_vs_late, showCategory = 10, title = "Early vs Late")


#to extract genes involved in specific pathway 
kegg_df <-kegg_early@result
kegg_df <- kegg_late@result

# Extract that pathway row
pathway_genes <- kegg_df[kegg_df$Description == "Cornified envelope formation", "geneID"]

#
genes_split <- unlist(strsplit(pathway_genes, "/"))
#convert to gene symbols
library(org.Hs.eg.db)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = genes_split,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)
#view genes
gene_symbols
#saving
write.csv(data.frame(Symbol = gene_symbols),
          "Cornified_envelope_late.csv",
          row.names = FALSE)

#GO anlaysis
library(clusterProfiler)
library(org.Hs.eg.db)

ego_early <- enrichGO(
  gene = early_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "MF",              # Biological Process
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
ego_late <- enrichGO(
  gene = late_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
ego_early_vs_late <- enrichGO(
  gene = early_vs_late_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(ego_early, showCategory = 20, title = "GO: Early vs Normal")
dotplot(ego_late, showCategory = 10, title = "GO: Late vs Normal")
dotplot(ego_early_vs_late, showCategory = 20, title = "GO: Early vs Late")
ego_early_vs_late

#to extract specific pathways 
go_df <- ego_late@result

genes_go <- go_df[go_df$Description == "serine-type peptidase activity", "geneID"]

#split into individual genes
genes_list <- unlist(strsplit(genes_go, "/"))
genes_list
write.csv(data.frame(Gene = genes_list),
          "MFserinetypelate_genes.csv",
          row.names = FALSE)
