library(DESeq2)
library(tidyverse)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(enrichplot)
library(ReactomePA)
library(clusterProfiler)
library(ReactomeGSA)
library(GSVA)
library(ReactomePA)
library(AnnotationDbi)
library(dplyr)
library(limma)
library(msigdbr)
library(pheatmap)
library(RColorBrewer)
library(ggVennDiagram)
library(ggvenn)


# ---------------- Paths ----------------
counts_file <- "C:/Users/User/Documents/Damola/counts/gene_counts.txt"  ## counts file
meta_file   <- "C:/Users/User/Documents/Damola/counts/samples.csv"  ## metadata samples file
out_dir     <- "C:/Users/User/Documents/Damola/deseq2_doxy_0.56-v2"    ## output directory
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

# ---------------- Load counts ----------------
fc <- read.delim(counts_file, comment.char = "#", check.names = FALSE)   ## read the featureCounts file

count_mat <- fc[, 7:ncol(fc)] ##extract only the count columns (column 7 onwards)
rownames(count_mat) <- fc$Geneid   ## set gene IDs as rownames

# Clean sample names to match the meta data
colnames(count_mat) <- sub(".*/", "", colnames(count_mat))
colnames(count_mat) <- sub("\\.bam$", "", colnames(count_mat))

count_mat <- as.matrix(count_mat)
storage.mode(count_mat) <- "integer"
write.csv(count_mat,"gene_counts_formatted.csv")
# ---------------- Load metadata ----------------
meta <- read.csv(meta_file, stringsAsFactors = FALSE)
meta <- meta[match(colnames(count_mat), meta$sample), ]
rownames(meta) <- meta$sample
## check to ensure alignment betwen samples and count matrix
stopifnot(all(colnames(count_mat) == meta$sample))

### subset to doxy treated only
doxy_samples <- meta %>% filter(treatment == "doxy") %>% select(sample,genotype,replicate)   ### select only doxy treated samples
count_doxy <- count_mat[, doxy_samples$sample]   ### extract the counts for the doxy treated samples
stopifnot(all(colnames(count_doxy) == doxy_samples$sample))
write.csv(doxy_samples,file="samples.csv",row.names = F)

# Factors and reference levels
doxy_samples$genotype  <- factor(doxy_samples$genotype,  levels = c("EV", "G0", "G1", "G2"))
doxy_samples$replicate <- factor(doxy_samples$replicate)

# DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_doxy, colData = doxy_samples, design = ~ replicate + genotype)
dds <- dds[rowSums(counts(dds)) >= 10, ] ### select only genes with more than 10 counts
dds <- DESeq(dds)
saveRDS(dds, "dds_full.rds")

# ==============================================================================
# HELPER FUNCTION: ANNOTATION OF ENSEMBL ID TO GENE SYMBOLS
# ==============================================================================
annotate_results <- function(res) {
  res_df <- as.data.frame(res) %>% rownames_to_column("gene_id")   ### move the gene id from row names to a new column
  gene_ids <- gsub("\\..*", "", res_df$gene_id) ##remove version number from the ensemble ID
  ### query org.Hs.eg.db for the gene symbol
  annot <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gene_ids,
    columns = c("SYMBOL"),
    keytype = "ENSEMBL"
  )
  ### jon symbol to the results
  res_annot <- res_df %>%
    mutate(ENSEMBL =  gene_ids) %>%
    left_join(annot, by = "ENSEMBL")
  
  res_annot
}

# ==============================================================================
# HELPER FUNCTION: VOLCANO
# ==============================================================================
plot_volcano <- function(res_annot, title, file) {
volc <- res_annot %>%
    mutate(
      neglog10padj = -log10(padj),
      sig = case_when(
        padj < 0.05 & abs(log2FoldChange) >= 0.56 ~ "Significant",
        TRUE ~ "Not significant"
      )
    )
df <- volc %>% filter(sig == "Significant")
write.csv(df,paste0(title,"_filtered.csv"),row.names = F)
y_cap <- 50
volc$neglog10padj_capped <- pmin(volc$neglog10padj, y_cap)

label_genes <- volc %>%
  filter(sig == "Significant", !is.na(SYMBOL)) %>%
  arrange(padj) %>%
  slice(1:20)

p <- ggplot(volc, aes(x = log2FoldChange, y = neglog10padj_capped)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 1.8) +
  scale_color_manual(values = c("Significant" = "red", "Not significant" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_vline(xintercept = c(-0.56, 0.56), linetype = "dashed", color = "grey") +
  geom_text_repel(
    data = label_genes,
    aes(label = SYMBOL),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    x = "log2 Fold Change",
    y = expression(-log[10]("adjusted p-value")),
    title = title,
    color = "Significance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(file, p, width = 6, height = 5, dpi = 300)
}

# ==============================================================================
# HELPER FUNCTION: PCA
# ==============================================================================
plot_pca <- function(dds, title, file) {
vsd <- vst(dds, blind = FALSE)
pca <- prcomp(t(assay(vsd)))

pc_df <- as.data.frame(pca$x) %>%
  rownames_to_column("sample")

meta_df <- as.data.frame(colData(dds))
meta_df$sample <- rownames(meta_df)
pc_df <- left_join(pc_df, meta_df, by = "sample")
percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
p <- ggplot(pc_df, aes(
  x = PC1, y = PC2,
  color = replicate,      
  shape =genotype
)) +
  geom_point(size = 3, alpha = 0.9) +
  
  labs(
    x = paste0("PC1 (", percentVar[1], "%)"),
    y = paste0("PC2 (", percentVar[2], "%)"),
    color = "Replicate",
    shape = "Genotype",
    title = title
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
ggsave(file, p, width = 6, height = 5, dpi = 300)
}

# ==============================================================================
# HELPER FUNCTION: HEATMAP
# ==============================================================================
# This function plots a heatmap of the top DE genes using VST‑normalized expression values, 
## scaled by row (Z‑score), to highlight relative up‑ and down‑regulation across samples
plot_heatmap <- function(dds, res_annot, title, file, n_genes = 50) {
  # 1. Select top genes
  top_sig <- res_annot %>%
    filter(!is.na(SYMBOL)) %>%
    #arrange(padj) %>%
    filter(padj < 0.05, abs(log2FoldChange) >= 0.56) %>% 
    head(n_genes)
  
  # Ensure genes exist in the VST matrix
  vsd <- vst(dds, blind = FALSE)
  vsd_mat <- assay(vsd)
  top_sig <- top_sig[top_sig$gene_id %in% rownames(vsd_mat), ]
  
  # 2. Extract expression matrix
  plot_mat <- vsd_mat[top_sig$gene_id, ]
  rownames(plot_mat) <- top_sig$SYMBOL
  
  # 3. Annotation
  anno_col <- as.data.frame(colData(dds)[, "genotype", drop = FALSE])

  # Drop unused levels in the annotation column
  anno_col$genotype <- droplevels(anno_col$genotype)

  # Full palette
  full_colors <- c(
    "EV" = "#D4AF37",
    "G0" = "#00CED1",
    "G1" = "brown",
    "G2" = "orange"
  )

  # Restrict palette to genotypes present in this subset
  present_genotypes <- levels(anno_col$genotype)
  anno_colors <- list(genotype = full_colors[present_genotypes])
  
  # 4. Plot
  png(file, width = 800, height = 1000, res = 120)
  pheatmap(plot_mat,
           main = title,
           annotation_col = anno_col,
           annotation_colors = anno_colors,
           scale = "row",
           clustering_method = "ward.D2",
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           show_colnames = FALSE,
           fontsize_row = 8)
  dev.off()
}


## ============================================================================================
### Differential analysis
## ==========================================================================================
#### Visualize clustering
plot_pca(dds,"Doxycline treated samples","pca_doxy_only_0.56.png")  ### plot PCA

res_G0_vs_EV <- results(dds, contrast = c("genotype", "G0", "EV"))
res_G1_vs_EV <- results(dds, contrast = c("genotype", "G1", "EV"))
res_G2_vs_EV <- results(dds, contrast = c("genotype", "G2", "EV"))

# G0 vs EV
res_G0_EV <- lfcShrink(dds,coef = "genotype_G0_vs_EV",res = res_G0_vs_EV,type = "apeglm")   ## shrink the log fold change
res_G0_EV_annot <- annotate_results(res_G0_EV)    ## annotate the result with the Gene symbol
write.csv(res_G0_EV_annot, "doxy_G0_vs_EV_results_0.56.csv", row.names = FALSE)  ### save to a csv file
plot_volcano(res_G0_EV_annot,"Doxycycline_treated_G0_vs_EV","volcano_doxy_G0_vs_EV_0.56.png")   ### plot volcano

# G1 vs EV
res_G1_EV <- lfcShrink(dds,coef = "genotype_G1_vs_EV",res = res_G1_vs_EV,type = "apeglm")   ## shrink the log fold change
res_G1_EV_annot <- annotate_results(res_G1_EV)    ## annotate the result with the Gene symbol
write.csv(res_G1_EV_annot, "doxy_G1_vs_EV_results_0.56.csv", row.names = FALSE)  ### save to a csv file
plot_volcano(res_G1_EV_annot,"Doxycycline_treated_G1_vs_EV","volcano_doxy_G1_vs_EV_0.56.png")   ### plot volcano


# G2 vs EV
res_G2_EV <- lfcShrink(dds,coef = "genotype_G2_vs_EV",res = res_G2_vs_EV,type = "apeglm")   ## shrink the log fold change
res_G2_EV_annot <- annotate_results(res_G2_EV)    ## annotate the result with the Gene symbol
write.csv(res_G2_EV_annot, "doxy_G2_vs_EV_results_0.56.csv", row.names = FALSE)  ### save to a csv file
plot_volcano(res_G2_EV_annot,"Doxycycline_treated_G2_vs_EV","volcano_doxy_G2_vs_EV_0.56.png")   ### plot volcano

###################################################
### Venn diagram
### Significant genes within doxy treated 
sig_G0 <- res_G0_EV_annot %>% filter(padj < 0.05, abs(log2FoldChange) >= 0.56) %>% pull(SYMBOL) %>% na.omit()
sig_G1 <- res_G1_EV_annot %>% filter(padj < 0.05, abs(log2FoldChange) >= 0.56) %>% pull(SYMBOL) %>% na.omit()
sig_G2 <- res_G2_EV_annot %>% filter(padj < 0.05, abs(log2FoldChange) >= 0.56) %>% pull(SYMBOL) %>% na.omit()
write.csv(sig_G0,"list_of_genes_in_Venn_plot_for_G0_vs_EV.csv")
write.csv(sig_G1,"list_of_genes_in_Venn_plot_for_G1_vs_EV.csv")
write.csv(sig_G2,"list_of_genes_in_Venn_plot_for_G2_vs_EV.csv")

### variations of venn diagram
### option 1
venn.plot <- venn.diagram( x = list(G0 = unique(sig_G0),G1 = unique(sig_G1),G2 = unique(sig_G2)),
  filename = "venn_doxy_0.56.png",
  fill = c("pink", "grey", "purple"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  main = "Overlap of DEGs in the doxycycline treated samples"
)

# Option 2. Create the plot
gene_list <- list(G0 = unique(sig_G0),G1 = unique(sig_G1),G2 = unique(sig_G2))
p_venn <- ggVennDiagram(gene_list, label_alpha = 0) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + # Professional Blue Gradient
  theme(legend.position = "none") +
  labs(title = "Overlap of DEGs in the doxycycline treated samples",
       subtitle = "Significant genes")
#Display and Save
ggsave("venn_doxy_treated.png", p_venn, width = 7, height = 6, dpi = 300)

## option 3
genotype_colors <- c("G0" = "#00CED1", "G1" = "pink", "G2" = "grey")
p_ggvenn <- ggvenn(
  gene_list, 
  fill_color = genotype_colors,
  stroke_size = 0.5, 
  set_name_size = 5,
  text_size = 4
) +
  ggtitle("Overlap of DEGs in the doxycycline treated samples")

##### https://www.sc-best-practices.org/introduction/prior_art.html
# ==============================================================================
# APPLYING HEATMAPS TO CREATE PLOTS
# ==============================================================================
# G2 vs EV
dds_G2_EV <- dds[, dds$genotype %in% c("G2", "EV")] ### create a subset
dds_G2_EV$genotype <- droplevels(dds_G2_EV$genotype)  ## droplevels removes factors that no longer appear in the subset
dds_G2_EV$replicate <- droplevels(dds_G2_EV$replicate)
plot_heatmap(dds = dds_G2_EV, res_annot  = res_G2_EV_annot, title = "DEGs G2_vs_EV", file = "heatmap_G2_vs_EV.png")

## G1 vs EV
dds_G1_EV <- dds[, dds$genotype %in% c("G1", "EV")] ### create a subset
dds_G1_EV$genotype <- droplevels(dds_G1_EV$genotype)  ## droplevels removes factors that no longer appear in the subset
dds_G1_EV$replicate <- droplevels(dds_G1_EV$replicate)
plot_heatmap(dds = dds_G1_EV,res_annot = res_G1_EV_annot, title = "DEGs G1 vs EV", file = "heatmap_G1_vs_EV.png")

## G0 vs EV
dds_G0_EV <- dds[, dds$genotype %in% c("G0", "EV")] ### create a subset
dds_G0_EV$genotype <- droplevels(dds_G0_EV$genotype)  ## droplevels removes factors that no longer appear in the subset
dds_G0_EV$replicate <- droplevels(dds_G0_EV$replicate)
plot_heatmap(dds = dds_G0_EV,res_annot  = res_G0_EV_annot, title = "DEGs G0 vs EV", file = "heatmap_G0vsEV.png")

## ====================
####### Genotype progression
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
rownames(mat) <- sub("\\..*", "", rownames(mat))  ### remove decimal values 
### convert to data frames
res_G0_df <- as.data.frame(res_G0_EV_annot)
res_G1_df <- as.data.frame(res_G1_EV_annot)
res_G2_df <- as.data.frame(res_G2_EV_annot)
### replace the gene_id column with the formatted values in the ENSEMBL column
res_G0_df$gene_id <- res_G0_df$ENSEMBL
res_G1_df$gene_id <- res_G1_df$ENSEMBL
res_G2_df$gene_id <- res_G2_df$ENSEMBL
### filter using p_adj and lfc
sig_G0 <- res_G0_df[ abs(res_G0_df$log2FoldChange) >= 0.56 &res_G0_df$padj < 0.05, ]
sig_G1 <- res_G1_df[abs(res_G1_df$log2FoldChange) >= 0.56 & res_G1_df$padj < 0.05, ]
sig_G2 <- res_G2_df[abs(res_G2_df$log2FoldChange) >= 0.56 & res_G2_df$padj < 0.05, ]
### combine into a single df
res_all <- dplyr::bind_rows(sig_G0, sig_G1, sig_G2)
## rank by statistical significance
top100 <- res_all %>% arrange(padj) %>% head(100)  ## sort by pdj in ascending order and take the first 100
#top100 <- res_all %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(100)
top100_unique <- top100[!duplicated(top100$SYMBOL), ]
heatmap_mat <- mat[top100_unique$gene_id, ]
rownames(heatmap_mat) <- top100_unique$SYMBOL
### reorder the columns to show progresion of DEGs
desired_order <- c("G2", "G1", "G0", "EV")
sample_info <- as.data.frame(colData(dds))
ordered_samples <- rownames(sample_info)[order(factor(sample_info$genotype, levels = desired_order))]
heatmap_mat <- heatmap_mat[, ordered_samples]

# Annotation for heatmap
anno_col <- as.data.frame(colData(dds)[, "genotype", drop = FALSE])
anno_col <- anno_col[ordered_samples,,drop = FALSE]
# Full palette
full_colors <- c("EV" = "#D4AF37","G0" = "#00CED1", "G1" = "brown", "G2" = "orange")
# Restrict palette to genotypes present in this subset
present_genotypes <- levels(anno_col$genotype)
anno_colors <- list(genotype = full_colors[present_genotypes])

## plot
file <- "genotype_progression_heatmap.png"
png(file, width = 800, height = 1000, res = 120)
pheatmap(
  heatmap_mat,
  scale = "row",  ## convert expression values to Zscore
  cluster_cols = FALSE,   # keeps the order of the genotypes
  cluster_rows = TRUE,
  annotation_col = anno_col,
  annotation_colors = anno_colors,
  clustering_method = "ward.D2",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  show_colnames = F,
  fontsize_row = 5,
  main = "Progression of DEGs in Doxycycline treated samples"
)
dev.off()


## ================================================================================
#  REACTOME PATHWAY ANALYSIS
# =================================================================================
# 1. Load Reactome pathways (gene symbols)
reactome_df <- msigdbr(species = "Homo sapiens",category = "C2",subcategory = "CP:REACTOME")
reactome_gs <- split(reactome_df$gene_symbol, reactome_df$gs_name)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================
# Convert ENSEMBL → SYMBOL
convert_to_symbol <- function(mat) {
  ens <- sub("\\..*", "", rownames(mat))

  mapping <- AnnotationDbi::select(org.Hs.eg.db,keys = ens,columns = "SYMBOL", keytype = "ENSEMBL")

  mapping <- mapping[!is.na(mapping$SYMBOL), ]
  mapping <- mapping[!duplicated(mapping$ENSEMBL), ]
  mapping <- mapping[mapping$ENSEMBL %in% ens, ]

  mat2 <- mat[match(mapping$ENSEMBL, ens), , drop = FALSE]
  rownames(mat2) <- mapping$SYMBOL
  mat2
}

# Collapse duplicate SYMBOL rows
collapse_duplicates <- function(mat) {
  df <- as.data.frame(mat)
  df$symbol <- rownames(df)

  df <- df %>%
    group_by(symbol) %>%
    summarise(across(everything(), sum)) %>%
    as.data.frame()

  rownames(df) <- df$symbol
  df$symbol <- NULL
  as.matrix(df)
}

# ==============================================================================
# HEATMAP FUNCTION FOR PATHWAY ACTIVITY
# ==============================================================================

plot_gsva_heatmap <- function(gsva_results, title, file) {
  png(file, width = 1000, height = 800, res = 120)
  # Select top 30 pathways by variance
  pathway_vars <- apply(gsva_results, 1, var)
  top_pathways <- gsva_results[order(pathway_vars, decreasing = TRUE)[1:30], ]
  pheatmap(top_pathways,main = title,scale = "row",color = colorRampPalette(c("blue", "white", "red"))(100),fontsize_row = 7)
  dev.off()
}
# ==============================================================================
# PREPARE EXPRESSION MATRIX FOR GSVA
# ==============================================================================
# Use VST-normalized counts (best practice for GSVA)
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
# Convert ENSEMBL → SYMBOL and collapse duplicates
mat_sym <- convert_to_symbol(mat)
mat_sym <- collapse_duplicates(mat_sym)
# ==============================================================================
# RUN GSVA (ONE TIME ONLY)
# ==============================================================================
param <- ssgseaParam(mat_sym, reactome_gs)
gsva_all <- GSVA::gsva(param)   # pathways × samples matrix
# ==============================================================================
# SUBSET GSVA SCORES BY GENOTYPE
# ==============================================================================
gsva_G0 <- gsva_all[, dds$genotype == "G0", drop = FALSE]
gsva_G1 <- gsva_all[, dds$genotype == "G1", drop = FALSE]
gsva_G2 <- gsva_all[, dds$genotype == "G2", drop = FALSE]
gsva_EV <- gsva_all[, dds$genotype == "EV", drop = FALSE]


# ==============================================================================
# PLOT GENOTYPE-SPECIFIC PATHWAY ACTIVITY
# ==============================================================================
plot_gsva_heatmap(gsva_EV, "Top Reactome Pathways: EV Activity", "gsva_heatmap_GL.png")
plot_gsva_heatmap(gsva_G2, "Top Reactome Pathways: G2 Activity", "gsva_heatmap_G2.png")
plot_gsva_heatmap(gsva_G1, "Top Reactome Pathways: G1 Activity", "gsva_heatmap_G1.png")
plot_gsva_heatmap(gsva_G0, "Top Reactome Pathways: G0 Activity", "gsva_heatmap_G0.png")

library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(enrichplot)
library(GSVA)
library(GSEABase)
library(tidyverse)
convert_ids <- function(df) {
  bitr(df$gene_id,
       fromType = "ENSEMBL",
       toType = "ENTREZID",
       OrgDb = org.Hs.eg.db) %>%
    pull(ENTREZID)
}

entrez_G0 <- convert_ids(sig_G0)
entrez_G1 <- convert_ids(sig_G1)
entrez_G2 <- convert_ids(sig_G2)

reactome_G0 <- enrichPathway(entrez_G0, organism = "human", pvalueCutoff = 0.05, readable = TRUE)
reactome_G1 <- enrichPathway(entrez_G1, organism = "human", pvalueCutoff = 0.05, readable = TRUE)
reactome_G2 <- enrichPathway(entrez_G2, organism = "human", pvalueCutoff = 0.05, readable = TRUE)

dotplot(reactome_G0, showCategory = 15) + ggtitle("Reactome – G0 vs EV")
dotplot(reactome_G1, showCategory = 15) + ggtitle("Reactome – G1 vs EV")
dotplot(reactome_G2, showCategory = 15) + ggtitle("Reactome – G2 vs EV")

## Enrich GO
ego_G0 <- enrichGO(entrez_G0, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP", pvalueCutoff = 0.05, readable = TRUE)
ego_G1 <- enrichGO(entrez_G1, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP", pvalueCutoff = 0.05, readable = TRUE)
ego_G2 <- enrichGO(entrez_G2, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP", pvalueCutoff = 0.05, readable = TRUE)
dotplot(ego_G0, showCategory = 15)
dotplot(ego_G1, showCategory = 15)
dotplot(ego_G2, showCategory = 15)

### KEGG pathway
kegg_G0 <- enrichKEGG(entrez_G0, organism = "hsa", pvalueCutoff = 0.05)
kegg_G1 <- enrichKEGG(entrez_G1, organism = "hsa", pvalueCutoff = 0.05)
kegg_G2 <- enrichKEGG(entrez_G2, organism = "hsa", pvalueCutoff = 0.05)
dotplot(kegg_G0, showCategory = 15)
dotplot(kegg_G1, showCategory = 15)
dotplot(kegg_G2, showCategory = 15)

## GSEA
rank_genes <- function(df) {
  # Convert ENSEMBL → ENTREZ
  conv <- bitr(df$gene_id,
               fromType = "ENSEMBL",
               toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
  
  # Merge ENTREZ IDs back into the results
  df2 <- df %>%
    dplyr::inner_join(conv, by = c("gene_id" = "ENSEMBL"))
  
  # Remove duplicated ENTREZ IDs
  df2 <- df2[!duplicated(df2$ENTREZID), ]
  
  # Create ranked list
  ranks <- df2$log2FoldChange
  names(ranks) <- df2$ENTREZID
  
  sort(ranks, decreasing = TRUE)
}

rank_G0 <- rank_genes(res_G0_df)
rank_G1 <- rank_genes(res_G1_df)
rank_G2 <- rank_genes(res_G2_df)

gsea_G0 <- gsePathway(rank_G0, organism = "human", pvalueCutoff = 0.05)
gsea_G1 <- gsePathway(rank_G1, organism = "human", pvalueCutoff = 0.05)
gsea_G2 <- gsePathway(rank_G2, organism = "human", pvalueCutoff = 0.05)

ridgeplot(gsea_G0)
ridgeplot(gsea_G1)
ridgeplot(gsea_G2)

# emapplot(gsea_G0)
# emapplot(gsea_G1)
# emapplot(gsea_G2)

####GSVA (Gene Set Variation Analysis)
reactome_gmt <- getGmt("ReactomePathways.gmt")
gs_list <- lapply(reactome_gmt, function(gs) {
  geneIds(gs)
})
names(gs_list) <- sapply(reactome_gmt, function(gs) gs@setName)
gsva_scores <- gsva(
  expr = as.matrix(mat),
  gset.idx.list = gs_list,
  method = "gsva",
  kcdf = "Gaussian"   # correct for VST data
)


pheatmap(gsva_scores,cluster_cols = FALSE,show_rownames = FALSE,main = "Reactome Pathway Activity (GSVA)")

































## ================================================================================
## 1. Load Reactome pathways (msigdbr, SYMBOL-based)
## ================================================================================
reactome_df <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcategory = "CP:REACTOME"
)

reactome_gs <- split(reactome_df$gene_symbol, reactome_df$gs_name)


## ================================================================================
## 2. Helper functions
## ================================================================================

### ENSEMBL → SYMBOL
convert_to_symbol <- function(mat) {
  ens <- sub("\\..*", "", rownames(mat))
  
  mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = ens,
    columns = "SYMBOL",
    keytype = "ENSEMBL"
  )
  
  # Remove NA symbols
  mapping <- mapping[!is.na(mapping$SYMBOL), ]
  
  # Keep only ENSEMBL IDs that exist in the matrix
  mapping <- mapping[mapping$ENSEMBL %in% ens, ]
  
  # Remove duplicated ENSEMBL IDs
  mapping <- mapping[!duplicated(mapping$ENSEMBL), ]
  
  # Reorder mapping to match matrix order
  mapping <- mapping[match(ens, mapping$ENSEMBL), ]
  mapping <- mapping[!is.na(mapping$ENSEMBL), ]
  
  # Subset matrix safely
  mat2 <- mat[mapping$ENSEMBL, , drop = FALSE]
  
  # Assign SYMBOL rownames
  rownames(mat2) <- mapping$SYMBOL
  
  mat2
}


### Collapse duplicate SYMBOL rows
collapse_duplicates <- function(mat) {
  df <- as.data.frame(mat)
  df$symbol <- rownames(df)
  
  df <- df %>%
    group_by(symbol) %>%
    summarise(across(everything(), sum)) %>%
    as.data.frame()
  
  rownames(df) <- df$symbol
  df$symbol <- NULL
  as.matrix(df)
}

### Convert ENSEMBL → ENTREZ (vector input)
convert_ids <- function(ensembl_ids) {
  conv <- bitr(
    ensembl_ids,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  unique(conv$ENTREZID)
}

### Create ranked vector for GSEA
rank_genes <- function(df) {
  conv <- bitr(
    df$gene_id,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  df2 <- df %>%
    inner_join(conv, by = c("gene_id" = "ENSEMBL")) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  
  ranks <- df2$stat
  names(ranks) <- df2$ENTREZID
  sort(ranks, decreasing = TRUE)
}


## ================================================================================
## 3. Prepare VST matrix for GSVA
## ================================================================================
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)

mat_sym <- convert_to_symbol(mat)
mat_sym <- collapse_duplicates(mat_sym)


## ================================================================================
## 4. Run ssGSEA (GSVA)
## ================================================================================
param <- ssgseaParam(mat_sym, reactome_gs)
gsva_all <- GSVA::gsva(param)


## ================================================================================
## 5. Subset GSVA scores by genotype
## ================================================================================
gsva_G0 <- gsva_all[, dds$genotype == "G0", drop = FALSE]
gsva_G1 <- gsva_all[, dds$genotype == "G1", drop = FALSE]
gsva_G2 <- gsva_all[, dds$genotype == "G2", drop = FALSE]
gsva_EV <- gsva_all[, dds$genotype == "EV", drop = FALSE]


## ================================================================================
## 6. Plot GSVA heatmaps
## ================================================================================
plot_gsva_heatmap <- function(gsva_results, title, file) {
  png(file, width = 1000, height = 800, res = 120)
  pathway_vars <- apply(gsva_results, 1, var)
  top_pathways <- gsva_results[order(pathway_vars, decreasing = TRUE)[1:30], ]
  pheatmap(
    top_pathways,
    main = title,
    scale = "row",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    fontsize_row = 7
  )
  dev.off()
}

plot_gsva_heatmap(gsva_EV, "Reactome Pathways – EV", "gsva_EV.png")
plot_gsva_heatmap(gsva_G0, "Reactome Pathways – G0", "gsva_G0.png")
plot_gsva_heatmap(gsva_G1, "Reactome Pathways – G1", "gsva_G1.png")
plot_gsva_heatmap(gsva_G2, "Reactome Pathways – G2", "gsva_G2.png")


## ================================================================================
## 7. Over-representation analysis (Reactome, GO, KEGG)
## ================================================================================

### DEG lists must contain ENSEMBL IDs
entrez_G0 <- convert_ids(sub("\\..*", "", sig_G0$gene_id))
entrez_G1 <- convert_ids(sub("\\..*", "", sig_G1$gene_id))
entrez_G2 <- convert_ids(sub("\\..*", "", sig_G2$gene_id))


### Universe = all genes tested in DESeq2

universe_ens <- sub("\\..*", "", rownames(dds))
valid_ens <- universe_ens[universe_ens %in% keys(org.Hs.eg.db, keytype = "ENSEMBL")]
universe_entrez <- convert_ids(valid_ens)



### ReactomePA
reactome_G0 <- enrichPathway(entrez_G0, organism = "human",
                             universe = universe_entrez,
                             pvalueCutoff = 0.05, readable = TRUE)

reactome_G1 <- enrichPathway(entrez_G1, organism = "human",
                             universe = universe_entrez,
                             pvalueCutoff = 0.05, readable = TRUE)

reactome_G2 <- enrichPathway(entrez_G2, organism = "human",
                             universe = universe_entrez,
                             pvalueCutoff = 0.05, readable = TRUE)


### GO Biological Process
ego_G0 <- enrichGO(entrez_G0, OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID", ont = "BP",
                   universe = universe_entrez,
                   pvalueCutoff = 0.05, readable = TRUE)

ego_G1 <- enrichGO(entrez_G1, OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID", ont = "BP",
                   universe = universe_entrez,
                   pvalueCutoff = 0.05, readable = TRUE)

ego_G2 <- enrichGO(entrez_G2, OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID", ont = "BP",
                   universe = universe_entrez,
                   pvalueCutoff = 0.05, readable = TRUE)


### KEGG
kegg_G0 <- enrichKEGG(entrez_G0, organism = "hsa",
                      universe = universe_entrez,
                      pvalueCutoff = 0.05)

kegg_G1 <- enrichKEGG(entrez_G1, organism = "hsa",
                      universe = universe_entrez,
                      pvalueCutoff = 0.05)

kegg_G2 <- enrichKEGG(entrez_G2, organism = "hsa",
                      universe = universe_entrez,
                      pvalueCutoff = 0.05)


## ================================================================================
## 8. GSEA (ReactomePA::gsePathway)
## ================================================================================
rank_G0 <- rank_genes(res_G0_df)
rank_G1 <- rank_genes(res_G1_df)
rank_G2 <- rank_genes(res_G2_df)

gsea_G0 <- gsePathway(rank_G0, organism = "human", pvalueCutoff = 0.25)
gsea_G1 <- gsePathway(rank_G1, organism = "human", pvalueCutoff = 0.25)
gsea_G2 <- gsePathway(rank_G2, organism = "human", pvalueCutoff = 0.25)








########### =====================--------------------final
## ================================================================================
## ================================================================================
## 0. Libraries (assumed installed)
## ================================================================================
library(DESeq2)
library(tidyverse)
library(msigdbr)
library(GSVA)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(pheatmap)
library(enrichplot)
library(AnnotationDbi)

## Objects assumed to already exist from your DESeq2 pipeline:
## dds, res_G0_EV, res_G1_EV, res_G2_EV, res_G0_EV_annot, res_G1_EV_annot, res_G2_EV_annot


## ================================================================================
## 1. Reactome gene sets in ENSEMBL (for GSVA + GSEA)
## ================================================================================
reactome_df <- msigdbr(
  species       = "Homo sapiens",
  collection    = "C2",
  subcollection = "CP:REACTOME"
)

reactome_gs_ens <- split(reactome_df$ensembl_gene, reactome_df$gs_name)


## ================================================================================
## 2. Helper functions
## ================================================================================

### ENSEMBL → ENTREZ using org.Hs.eg.db (best effort)
convert_ids <- function(ensembl_ids) {
  ensembl_ids <- sub("\\..*", "", ensembl_ids)
  conv <- suppressWarnings(
    bitr(
      ensembl_ids,
      fromType = "ENSEMBL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db
    )
  )
  unique(na.omit(conv$ENTREZID))
}

### Ranked vector for GSEA (ReactomePA::gsePathway) using ENSEMBL directly
## df must contain:
##   - gene_id (ENSEMBL, possibly with version)
##   - stat   (DESeq2 Wald statistic)
rank_genes <- function(df) {
  ens_clean <- sub("\\..*", "", df$gene_id)
  ranks <- df$stat
  names(ranks) <- ens_clean
  sort(ranks, decreasing = TRUE)
}

### GSVA heatmap
plot_gsva_heatmap <- function(gsva_results, title, file) {
  png(file, width = 1000, height = 800, res = 120)
  pathway_vars <- apply(gsva_results, 1, var)
  top_n <- min(30, nrow(gsva_results))
  top_pathways <- gsva_results[order(pathway_vars, decreasing = TRUE)[1:top_n], ]
  pheatmap(
    top_pathways,
    main  = title,
    scale = "row",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    fontsize_row = 7
  )
  dev.off()
}


## ================================================================================
## 3. Prepare VST matrix for GSVA (ENSEMBL rownames)
## ================================================================================
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)

## strip ENSEMBL versions
rownames(mat) <- sub("\\..*", "", rownames(mat))


## ================================================================================
## 4. Run ssGSEA (GSVA) with ENSEMBL gene sets
## ================================================================================
param    <- ssgseaParam(mat, reactome_gs_ens)
gsva_all <- GSVA::gsva(param)


## ================================================================================
## 5. Subset GSVA scores by genotype
## ================================================================================
gsva_G0 <- gsva_all[, dds$genotype == "G0", drop = FALSE]
gsva_G1 <- gsva_all[, dds$genotype == "G1", drop = FALSE]
gsva_G2 <- gsva_all[, dds$genotype == "G2", drop = FALSE]
gsva_EV <- gsva_all[, dds$genotype == "EV", drop = FALSE]


## ================================================================================
## 6. Plot GSVA heatmaps
## ================================================================================
plot_gsva_heatmap(gsva_EV, "Reactome ssGSEA – EV", "gsva_EV.png")
plot_gsva_heatmap(gsva_G0, "Reactome ssGSEA – G0", "gsva_G0.png")
plot_gsva_heatmap(gsva_G1, "Reactome ssGSEA – G1", "gsva_G1.png")
plot_gsva_heatmap(gsva_G2, "Reactome ssGSEA – G2", "gsva_G2.png")


## ================================================================================
## 7. Universe for ORA (Reactome, GO, KEGG) in ENTREZ
## ================================================================================
universe_ens <- rownames(mat)  ## already version-stripped ENSEMBL
valid_ens    <- universe_ens[universe_ens %in% keys(org.Hs.eg.db, keytype = "ENSEMBL")]
universe_entrez <- convert_ids(valid_ens)


## ================================================================================
## 8. DEG lists in ENSEMBL → ENTREZ
## ================================================================================
## from your annotated results: res_G0_EV_annot, res_G1_EV_annot, res_G2_EV_annot
sig_G0_ens <- res_G0_EV_annot %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 0.56) %>%
  pull(ENSEMBL) %>%
  na.omit()

sig_G1_ens <- res_G1_EV_annot %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 0.56) %>%
  pull(ENSEMBL) %>%
  na.omit()

sig_G2_ens <- res_G2_EV_annot %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 0.56) %>%
  pull(ENSEMBL) %>%
  na.omit()

entrez_G0 <- convert_ids(sig_G0_ens)
entrez_G1 <- convert_ids(sig_G1_ens)
entrez_G2 <- convert_ids(sig_G2_ens)


## ================================================================================
## 9. Over-representation analysis (Reactome, GO BP, KEGG)
## ================================================================================

### ReactomePA ORA
reactome_G0 <- enrichPathway(
  gene         = entrez_G0,
  universe     = universe_entrez,
  organism     = "human",
  pvalueCutoff = 0.05,
  readable     = TRUE
)

reactome_G1 <- enrichPathway(
  gene         = entrez_G1,
  universe     = universe_entrez,
  organism     = "human",
  pvalueCutoff = 0.05,
  readable     = TRUE
)

reactome_G2 <- enrichPathway(
  gene         = entrez_G2,
  universe     = universe_entrez,
  organism     = "human",
  pvalueCutoff = 0.05,
  readable     = TRUE
)

### GO BP ORA
ego_G0 <- enrichGO(
  gene         = entrez_G0,
  universe     = universe_entrez,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  pvalueCutoff = 0.05,
  readable     = TRUE
)

ego_G1 <- enrichGO(
  gene         = entrez_G1,
  universe     = universe_entrez,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  pvalueCutoff = 0.05,
  readable     = TRUE
)

ego_G2 <- enrichGO(
  gene         = entrez_G2,
  universe     = universe_entrez,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  pvalueCutoff = 0.05,
  readable     = TRUE
)

### KEGG ORA
kegg_G0 <- enrichKEGG(
  gene         = entrez_G0,
  universe     = universe_entrez,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

kegg_G1 <- enrichKEGG(
  gene         = entrez_G1,
  universe     = universe_entrez,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

kegg_G2 <- enrichKEGG(
  gene         = entrez_G2,
  universe     = universe_entrez,
  organism     = "hsa",
  pvalueCutoff = 0.05
)


## ================================================================================
## 10. GSEA (ReactomePA::gsePathway) using ENSEMBL IDs directly
## ================================================================================
## Build DEG tables with ENSEMBL gene_id + stat
res_G0_df <- as.data.frame(res_G0_EV) %>% rownames_to_column("gene_id")
res_G1_df <- as.data.frame(res_G1_EV) %>% rownames_to_column("gene_id")
res_G2_df <- as.data.frame(res_G2_EV) %>% rownames_to_column("gene_id")

rank_G0 <- rank_genes(res_G0_df)
rank_G1 <- rank_genes(res_G1_df)
rank_G2 <- rank_genes(res_G2_df)

gsea_G0 <- gsePathway(rank_G0, organism = "human", pvalueCutoff = 0.25)
gsea_G1 <- gsePathway(rank_G1, organism = "human", pvalueCutoff = 0.25)
gsea_G2 <- gsePathway(rank_G2, organism = "human", pvalueCutoff = 0.25)


## ================================================================================
## 11. Visualization (ORA + GSEA)
## ================================================================================

## 11.1 Dotplots
dotplot(reactome_G0, showCategory = 20) + ggtitle("Reactome ORA – G0 vs EV")
dotplot(reactome_G1, showCategory = 20) + ggtitle("Reactome ORA – G1 vs EV")
dotplot(reactome_G2, showCategory = 20) + ggtitle("Reactome ORA – G2 vs EV")

dotplot(ego_G0, showCategory = 20) + ggtitle("GO BP ORA – G0 vs EV")
dotplot(ego_G1, showCategory = 20) + ggtitle("GO BP ORA – G1 vs EV")
dotplot(ego_G2, showCategory = 20) + ggtitle("GO BP ORA – G2 vs EV")

dotplot(kegg_G0, showCategory = 20) + ggtitle("KEGG ORA – G0 vs EV")
dotplot(kegg_G1, showCategory = 20) + ggtitle("KEGG ORA – G1 vs EV")
dotplot(kegg_G2, showCategory = 20) + ggtitle("KEGG ORA – G2 vs EV")

## 11.2 Enrichment maps (example for G0)
emapplot(pairwise_termsim(reactome_G0))
emapplot(pairwise_termsim(ego_G0))
emapplot(pairwise_termsim(kegg_G0))

## 11.3 Cnetplots (example for G0)
cnetplot(reactome_G0, showCategory = 10)
cnetplot(ego_G0,       showCategory = 10)
cnetplot(kegg_G0,      showCategory = 10)

## 11.4 GSEA ridgeplots
ridgeplot(gsea_G0) + ggtitle("Reactome GSEA – G0 vs EV")
ridgeplot(gsea_G1) + ggtitle("Reactome GSEA – G1 vs EV")
ridgeplot(gsea_G2) + ggtitle("Reactome GSEA – G2 vs EV")

## 11.5 GSEA running score example
gseaplot2(gsea_G0, geneSetID = 1)

