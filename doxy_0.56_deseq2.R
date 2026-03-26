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
out_dir     <- "C:/Users/User/Documents/Damola/deseq2_doxy_0.56"    ## output directory
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

# Factors and reference levels
doxy_samples$genotype  <- factor(doxy_samples$genotype,  levels = c("GL", "G0", "G1", "G2"))
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
plot_heatmap <- function(dds, res_annot, title, file, n_genes = 20) {
  # 1. Select top genes
  top_sig <- res_annot %>%
    filter(!is.na(SYMBOL)) %>%
    arrange(padj) %>%
    head(n_genes)
  
  # Ensure genes exist in the VST matrix
  vsd <- vst(dds, blind = FALSE)
  vsd_mat <- assay(vsd)
  top_sig <- top_sig[top_sig$gene_id %in% rownames(vsd_mat), ]
  
  # 2. Extract expression matrix
  plot_mat <- vsd_mat[top_sig$gene_id, ]
  rownames(plot_mat) <- top_sig$SYMBOL
  
  # 3. Annotation
  # 3. Annotation
  anno_col <- as.data.frame(colData(dds)[, "genotype", drop = FALSE])
  
  # Drop unused levels in the annotation column
  anno_col$genotype <- droplevels(anno_col$genotype)
  
  # Full palette
  full_colors <- c(
    "GL" = "#D4AF37",
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
plot_pca(dds,"PCA: doxy-only samples","pca_doxy_only_0.56.png")  ### plot PCA


res_G0_vs_GL <- results(dds, contrast = c("genotype", "G0", "GL"))
res_G1_vs_GL <- results(dds, contrast = c("genotype", "G1", "GL"))
res_G2_vs_GL <- results(dds, contrast = c("genotype", "G2", "GL"))

# G0 vs GL
res_G0_GL <- lfcShrink(dds,coef = "genotype_G0_vs_GL",res = res_G0_vs_GL,type = "apeglm")   ## shrink the log fold change
res_G0_GL_annot <- annotate_results(res_G0_GL)    ## annotate the result with the Gene symbol
write.csv(res_G0_GL_annot, "doxy_G0_vs_GL_results_0.56.csv", row.names = FALSE)  ### save to a csv file
plot_volcano(res_G0_GL_annot,"Volcano: doxy-treated, G0 vs GL","volcano_doxy_G0_vs_GL_0.56.png")   ### plot volcano

# G1 vs GL
res_G1_GL <- lfcShrink(dds,coef = "genotype_G1_vs_GL",res = res_G1_vs_GL,type = "apeglm")   ## shrink the log fold change
res_G1_GL_annot <- annotate_results(res_G1_GL)    ## annotate the result with the Gene symbol
write.csv(res_G1_GL_annot, "doxy_G1_vs_GL_results_0.56.csv", row.names = FALSE)  ### save to a csv file
plot_volcano(res_G1_GL_annot,"Volcano: doxy-treated, G1 vs GL","volcano_doxy_G1_vs_GL_0.56.png")   ### plot volcano

# G2 vs GL
res_G2_GL <- lfcShrink(dds,coef = "genotype_G2_vs_GL",res = res_G2_vs_GL,type = "apeglm")   ## shrink the log fold change
res_G2_GL_annot <- annotate_results(res_G2_GL)    ## annotate the result with the Gene symbol
write.csv(res_G2_GL_annot, "doxy_G2_vs_GL_results_0.56.csv", row.names = FALSE)  ### save to a csv file
plot_volcano(res_G2_GL_annot,"Volcano: doxy-treated, G2 vs GL","volcano_doxy_G2_vs_GL_0.56.png")   ### plot volcano

###################################################
### Venn diagram
### Significant genes within doxy treated 
sig_G0 <- res_G0_GL_annot %>% filter(padj < 0.05, abs(log2FoldChange) >= 0.56) %>% pull(SYMBOL) %>% na.omit()
sig_G1 <- res_G1_GL_annot %>% filter(padj < 0.05, abs(log2FoldChange) >= 0.56) %>% pull(SYMBOL) %>% na.omit()
sig_G2 <- res_G2_GL_annot %>% filter(padj < 0.05, abs(log2FoldChange) >= 0.56) %>% pull(SYMBOL) %>% na.omit()


### variations of venn diagram
### option 1
venn.plot <- venn.diagram( x = list(G0 = unique(sig_G0),G1 = unique(sig_G1),G2 = unique(sig_G2)),
  filename = "venn_doxy_0.56.png",
  fill = c("pink", "grey", "purple"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  main = "Overlap of DE genes (doxy-only)"
)

# Option 2. Create the plot
gene_list <- list(G0 = unique(sig_G0),G1 = unique(sig_G1),G2 = unique(sig_G2))
p_venn <- ggVennDiagram(gene_list, label_alpha = 0) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + # Professional Blue Gradient
  theme(legend.position = "none") +
  labs(title = "Overlap of DEGs (Doxy-only)",
       subtitle = "Significant genes")
#Display and Save
print(p_venn)
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
  ggtitle("Overlap of DE genes (doxy-only)")
print(p_ggvenn)


# ==============================================================================
# APPLYING HEATMAPS TO CREATE PLOTS
# ==============================================================================
# G2 vs GL
dds_G2_GL <- dds[, dds$genotype %in% c("G2", "GL")] ### create a subset
dds_G2_GL$genotype <- droplevels(dds_G2_GL$genotype)  ## droplevels removes factors that no longer appear in the subset
dds_G2_GL$replicate <- droplevels(dds_G2_GL$replicate)
### plot
plot_heatmap(dds = dds_G2_GL, 
             res_annot  = res_G2_GL_annot, 
             title      = "Heatmap: G2 vs GL @ 0.56 lfc", 
             file       = "heatmap_G2vsGL.png")

## G1 vs GL
dds_G1_GL <- dds[, dds$genotype %in% c("G1", "GL")] ### create a subset
dds_G1_GL$genotype <- droplevels(dds_G1_GL$genotype)  ## droplevels removes factors that no longer appear in the subset
dds_G1_GL$replicate <- droplevels(dds_G1_GL$replicate)
plot_heatmap(dds = dds_G1_GL, 
             res_annot  = res_G1_GL_annot, 
             title      = "Heatmap: G1 vs GL @ 0.56 lfc", 
             file       = "heatmap_G1vsGL.png")

## G0 vs GL
dds_G0_GL <- dds[, dds$genotype %in% c("G0", "GL")] ### create a subset
dds_G0_GL$genotype <- droplevels(dds_G0_GL$genotype)  ## droplevels removes factors that no longer appear in the subset
dds_G0_GL$replicate <- droplevels(dds_G0_GL$replicate)
plot_heatmap(dds = dds_G0_GL, 
             res_annot  = res_G0_GL_annot, 
             title      = "Heatmap: G0 vs GL @ 0.56 lfc", 
             file       = "heatmap_G0vsGL.png")



### ================================================================================
##  REACTOME PATHWAY ANALYSIS
## =================================================================================

# 1. Load Reactome pathways (gene symbols)
reactome_df <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcategory = "CP:REACTOME"
)

reactome_gs <- split(reactome_df$gene_symbol, reactome_df$gs_name)


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Convert ENSEMBL → SYMBOL
convert_to_symbol <- function(mat) {
  ens <- sub("\\..*", "", rownames(mat))
  
  mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = ens,
    columns = "SYMBOL",
    keytype = "ENSEMBL"
  )
  
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
  
  pheatmap(
    top_pathways,
    main = title,
    scale = "row",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    fontsize_row = 7
  )
  
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
gsva_all <- gsva(param)   # pathways × samples matrix


# ==============================================================================
# SUBSET GSVA SCORES BY GENOTYPE
# ==============================================================================

gsva_G0 <- gsva_all[, dds$genotype == "G0", drop = FALSE]
gsva_G1 <- gsva_all[, dds$genotype == "G1", drop = FALSE]
gsva_G2 <- gsva_all[, dds$genotype == "G2", drop = FALSE]
gsva_GL <- gsva_all[, dds$genotype == "GL", drop = FALSE]


# ==============================================================================
# PLOT GENOTYPE-SPECIFIC PATHWAY ACTIVITY
# ==============================================================================
plot_gsva_heatmap(gsva_GL, "Top Reactome Pathways: GL Activity", "gsva_heatmap_GL.png")
plot_gsva_heatmap(gsva_G2, "Top Reactome Pathways: G2 Activity", "gsva_heatmap_G2.png")
plot_gsva_heatmap(gsva_G1, "Top Reactome Pathways: G1 Activity", "gsva_heatmap_G1.png")
plot_gsva_heatmap(gsva_G0, "Top Reactome Pathways: G0 Activity", "gsva_heatmap_G0.png")
