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
out_dir     <- "C:/Users/User/Documents/Damola/deseq2_out_0.56"    ## output directory
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

stopifnot(all(colnames(count_mat) == meta$sample))

# Factors and reference levels
meta$treatment <- factor(meta$treatment, levels = c("untreated", "doxy"))
meta$genotype  <- factor(meta$genotype,  levels = c("GL", "G0", "G1", "G2"))

# Main DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta,
  design    = ~ genotype + treatment + genotype:treatment
)
dds <- DESeq(dds)
### select only genes with more than 10 counts
dds <- dds[rowSums(counts(dds)) >= 10, ]
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
    mutate(ENSEMBL = gsub("\\..*", "", gene_id)) %>%
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
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
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
plot_pca <- function(dds_subset, title, file) {
vsd <- vst(dds_subset, blind = FALSE)
pca <- prcomp(t(assay(vsd)))

pc_df <- as.data.frame(pca$x) %>%
  rownames_to_column("sample")

meta_df <- as.data.frame(colData(dds_subset))
meta_df$sample <- rownames(meta_df)

pc_df <- left_join(pc_df, meta_df, by = "sample")
# Convert replicate to factor so it can be used in a legend
pc_df$replicate <- as.factor(pc_df$replicate)

percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

p <- ggplot(pc_df, aes(
  x = PC1, y = PC2,
  color = treatment,      
  shape =genotype
)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(aes(label = replicate), size = 3, max.overlaps = Inf) +
  
  # Dummy layer to create replicate legend
  geom_point(aes(fill = replicate), shape = 21, alpha = 0) +
  
  labs(
    x = paste0("PC1 (", percentVar[1], "%)"),
    y = paste0("PC2 (", percentVar[2], "%)"),
    color = "Treatment",
    shape = "Genotype",
    fill = "Replicates",
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
# This function uses Z-score scaling to show 50 Top Genes
plot_heatmap <- function(dds_subset, res_annot, title, file, n_genes = 50) {
  
  # 1. Select top genes by significance (padj)
  # Filter out NAs to ensure we only have valid symbols for the plot
  top_sig <- res_annot %>%
    filter(!is.na(SYMBOL)) %>%
    arrange(padj) %>%
    head(n_genes)
  
  # 2. Prepare Expression Data (VST)
  vsd <- vst(dds_subset, blind = FALSE)
  vsd_mat <- assay(vsd)
  
  # 3. Match Matrix IDs with Symbols from res_annot
  # We subset the matrix to only the top genes identified
  plot_mat <- vsd_mat[top_sig$gene_id, ]
  rownames(plot_mat) <- top_sig$SYMBOL
  
  # 4. Define Annotation Colors (Matches your PCA colors)
  anno_col <- as.data.frame(colData(dds_subset)[, c("genotype", "treatment"), drop=FALSE])
  anno_colors <- list(
    genotype = c("GL" = "#D4AF37", "G0" = "#00CED1", "G1" = "red", "G2" = "black"),
    treatment = c("untreated" = "grey80", "doxy" = "purple")
  )
  
  # 5. Save Heatmap
  # pheatmap uses a file device (png/pdf) directly
  png(file, width = 800, height = 1000, res = 120)
  pheatmap(plot_mat, 
           main = title,
           annotation_col = anno_col,
           annotation_colors = anno_colors,
           scale = "row",                # Essential for Blue-to-Red relative view
           clustering_method = "ward.D2", 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           show_colnames = FALSE,
           fontsize_row = 8)
  dev.off()
}

## ============================================================================================
### Differential analysis
## Subset to doxy only
## genotype‑specific differential expression under doxycycline only, using GL as the reference
dds_doxy <- dds[, dds$treatment == "doxy"]
dds_doxy$treatment <- droplevels(dds_doxy$treatment)

# set GL as reference
dds_doxy$genotype <- relevel(dds_doxy$genotype, ref = "GL")

design(dds_doxy) <- ~ genotype  #### treatment is no longer  a variable since they are all treated now
dds_doxy <- DESeq(dds_doxy)
resultsNames(dds_doxy)
saveRDS(dds_doxy, "dds_doxy_only_0.56.rds")

# G0 vs GL
res_G0_GL_doxy <- results(dds_doxy, contrast = c("genotype", "G0", "GL"))   ### compare G0 to GL
res_G0_GL_doxy <- lfcShrink(dds_doxy,coef = "genotype_G0_vs_GL",res = res_G0_GL_doxy,type = "apeglm")   ## shrink the log fold change
res_G0_GL_doxy_annot <- annotate_results(res_G0_GL_doxy)    ## annotate the result with the Gene symbol
write.csv(res_G0_GL_doxy_annot, "doxy_G0_vs_GL_results_0.56.csv", row.names = FALSE)  ### save to a csv file
plot_volcano(res_G0_GL_doxy_annot,"Volcano: doxy, G0 vs GL","volcano_doxy_G0_vs_GL_0.56.png")   ### plot volcano
plot_pca(dds_doxy,"PCA: doxy-only samples","pca_doxy_only_0.56.png")  ### plot PCA

# G1 vs GL
res_G1_GL_doxy <- results(dds_doxy, contrast = c("genotype", "G1", "GL"))
res_G1_GL_doxy <- lfcShrink(dds_doxy, coef = "genotype_G1_vs_GL", res = res_G1_GL_doxy, type = "apeglm")
res_G1_GL_doxy_annot <- annotate_results(res_G1_GL_doxy)
write.csv(res_G1_GL_doxy_annot, "doxy_G1_vs_GL_results_0.56.csv", row.names = FALSE)
plot_volcano(res_G1_GL_doxy_annot,"Volcano: doxy, G1 vs GL","volcano_doxy_G1_vs_GL_0.56.png")

# G2 vs GL
res_G2_GL_doxy <- results(dds_doxy, contrast = c("genotype", "G2", "GL"))
res_G2_GL_doxy <- lfcShrink(dds_doxy, coef = "genotype_G2_vs_GL", res = res_G2_GL_doxy, type = "apeglm")
res_G2_GL_doxy_annot <- annotate_results(res_G2_GL_doxy)
write.csv(res_G2_GL_doxy_annot, "doxy_G2_vs_GL_results_0.56.csv", row.names = FALSE)
plot_volcano(res_G2_GL_doxy_annot,"Volcano: doxy, G2 vs GL","volcano_doxy_G2_vs_GL_0.56.png")

## ============================================================
### Comparing each treated to untreated within each genotype
## ==============================================================
# ---- G0: doxy vs untreated ----
dds_G0 <- dds[, dds$genotype == "G0"]
dds_G0$genotype <- droplevels(dds_G0$genotype)   ### isolates the treatment effect within the G0 genotype
### refine the design
design(dds_G0) <- ~ treatment   ### genotype is not pat of the model anymore since the subset all belong to G0
dds_G0 <- DESeq(dds_G0)
saveRDS(dds_G0, "dds_G0_0.56.rds")
##Extract the differential expression results (i.e. doxy - untreated)
res_G0_doxy_vs_untr <- results(dds_G0, contrast = c("treatment", "doxy", "untreated"))
res_G0_doxy_vs_untr <- lfcShrink(dds_G0,coef = "treatment_doxy_vs_untreated", res = res_G0_doxy_vs_untr,type = "apeglm")
### annotate result with gene symbols
res_G0_doxy_vs_untr_annot <- annotate_results(res_G0_doxy_vs_untr)
write.csv(res_G0_doxy_vs_untr_annot, "G0_doxy_vs_untreated_results_0.56.csv", row.names = FALSE)
### plots
plot_volcano(res_G0_doxy_vs_untr_annot, "Volcano: G0, doxy vs untreated", "volcano_G0_doxy_vs_untreated_0.56.png")
plot_pca(dds_G0,"PCA: G0 samples","pca_G0_0.56.png")

#====================================================
## ---- G1: doxy vs untreated ----
dds_G1 <- dds[, dds$genotype == "G1"]
dds_G1$genotype <- droplevels(dds_G1$genotype)
design(dds_G1) <- ~ treatment
dds_G1 <- DESeq(dds_G1)
saveRDS(dds_G1, "dds_G1_0.56.rds")

res_G1_doxy_vs_untr <- results(dds_G1, contrast = c("treatment", "doxy", "untreated"))
res_G1_doxy_vs_untr <- lfcShrink(dds_G1, coef = "treatment_doxy_vs_untreated", res = res_G1_doxy_vs_untr,type = "apeglm")

res_G1_doxy_vs_untr_annot <- annotate_results(res_G1_doxy_vs_untr)
write.csv(res_G1_doxy_vs_untr_annot, "G1_doxy_vs_untreated_results_0.56.csv", row.names = FALSE)

plot_volcano(res_G1_doxy_vs_untr_annot, "Volcano: G1, doxy vs untreated","volcano_G1_doxy_vs_untreated_0.56.png")

plot_pca(dds_G1,"PCA: G1 samples", "pca_G1_0.56.png")


###  ==================================================
# ---- G2: doxy vs untreated ----
dds_G2 <- dds[, dds$genotype == "G2"]
dds_G2$genotype <- droplevels(dds_G2$genotype)

design(dds_G2) <- ~ treatment
dds_G2 <- DESeq(dds_G2)
saveRDS(dds_G2, "dds_G2_0.56.rds")

res_G2_doxy_vs_untr <- results(dds_G2, contrast = c("treatment", "doxy", "untreated"))
res_G2_doxy_vs_untr <- lfcShrink(dds_G2, coef = "treatment_doxy_vs_untreated", res = res_G2_doxy_vs_untr, type = "apeglm")

res_G2_doxy_vs_untr_annot <- annotate_results(res_G2_doxy_vs_untr)
write.csv(res_G2_doxy_vs_untr_annot, "G2_doxy_vs_untreated_results_0.56.csv", row.names = FALSE)

plot_volcano(res_G2_doxy_vs_untr_annot, "Volcano: G2, doxy vs untreated", "volcano_G2_doxy_vs_untreated_0.56.png")

plot_pca(dds_G2,"PCA: G2 samples", "pca_G2_0.56.png")


###################################################
### Venn diagram
### Significant genes within doxy treated 
sig_G0 <- res_G0_GL_doxy_annot %>% filter(padj < 0.05) %>% pull(SYMBOL) %>% na.omit()
sig_G1 <- res_G1_GL_doxy_annot %>% filter(padj < 0.05) %>% pull(SYMBOL) %>% na.omit()
sig_G2 <- res_G2_GL_doxy_annot %>% filter(padj < 0.05) %>% pull(SYMBOL) %>% na.omit()


### variations of venn diagram
### option 1
venn.plot <- venn.diagram( x = list(G0 = sig_G0,G1 = sig_G1,G2 = sig_G2),
  filename = "venn_doxy_0.56.png",
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  main = "Overlap of DE genes (doxy-only)"
)
gene_list <- list(G0 = sig_G0, G1 = sig_G1, G2 = sig_G2)

# Option 2. Create the plot
p_venn <- ggVennDiagram(gene_list, label_alpha = 0) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + # Professional Blue Gradient
  theme(legend.position = "none") +
  labs(title = "Overlap of DEGs (Doxy-only)",
       subtitle = "Significant genes (padj < 0.05) vs GL")
#Display and Save
print(p_venn)
ggsave("venn_modern_doxy.png", p_venn, width = 7, height = 6, dpi = 300)

## option 3
# Use your specific genotype colors from the PCA
genotype_colors <- c("G0" = "#00CED1", "G1" = "red", "G2" = "black")

p_ggvenn <- ggvenn(
  gene_list, 
  fill_color = genotype_colors,
  stroke_size = 0.5, 
  set_name_size = 5,
  text_size = 4
) +
  ggtitle("Overlap of DE genes (doxy-only)")

print(p_ggvenn)



### significant genes within the genotype (treated vs untreated)
sig_G0_treat <- res_G0_doxy_vs_untr_annot %>% filter(padj < 0.05) %>% pull(SYMBOL) %>% na.omit()
sig_G1_treat <- res_G1_doxy_vs_untr_annot %>% filter(padj < 0.05) %>% pull(SYMBOL) %>% na.omit()
sig_G2_treat <- res_G2_doxy_vs_untr_annot %>% filter(padj < 0.05) %>% pull(SYMBOL) %>% na.omit()

venn.diagram(
  x = list(
    G0 = sig_G0_treat,
    G1 = sig_G1_treat,
    G2 = sig_G2_treat
  ),
  filename = "venn_treated_vs_untreated_0.56.png",
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  main = "Overlap of DE genes (treated vs untreated)"
)
gene_list2 = list(
  G0 = sig_G0_treat,
  G1 = sig_G1_treat,
  G2 = sig_G2_treat
)
# Option 2. Create the plot
p_venn <- ggVennDiagram(gene_list2, label_alpha = 0) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + # Professional Blue Gradient
  theme(legend.position = "none") +
  labs(title = "Overlap of DEGs (Doxy vs untreated)",
       subtitle = "Significant genes (padj < 0.05) vs GL")
#Display and Save
print(p_venn)
ggsave("venn_treated_vs_untreated_0.56_modern.png", p_venn, width = 7, height = 6, dpi = 300)

## option 3
# Use your specific genotype colors from the PCA
genotype_colors <- c("G0" = "#00CED1", "G1" = "red", "G2" = "black")

p_ggvenn <- ggvenn(
  gene_list2, 
  fill_color = genotype_colors,
  stroke_size = 0.5, 
  set_name_size = 5,
  text_size = 4
) +
  ggtitle("Overlap of DE genes (doxy vs untreated)")

print(p_ggvenn)
ggsave("venn_treated_vs_untreated_0.56_variation3.png", p_venn, width = 7, height = 6, dpi = 300)



### ===========================================================================================================
### Reactome
### Reactome on doxy response within each genotype
### significant genes
sig_G0_entrez <- bitr(sig_G0_treat, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID  ### convert gene symbol to entrez ID
sig_G1_entrez <- bitr(sig_G1_treat, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
sig_G2_entrez <- bitr(sig_G2_treat, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID

### G0 reactome
react_G0 <- enrichPathway(gene = sig_G0_entrez,organism = "human",pvalueCutoff = 0.05,readable = TRUE)

### G1 
react_G1 <- enrichPathway(gene = sig_G1_entrez,organism = "human",pvalueCutoff = 0.05,readable = TRUE)

### G2
react_G2 <- enrichPathway(gene = sig_G2_entrez,organism = "human",pvalueCutoff = 0.05,readable = TRUE)

### dotplot
### ==> toDO : improve graphics
p_G0 <- dotplot(react_G0, showCategory = 20) + ggtitle("Reactome: G0 doxy vs untreated")
ggsave("reactome_G0_dotplot_doxy_vs_untreated.png", p_G0, width = 7, height = 6, dpi = 300)
p_G1 <- dotplot(react_G1, showCategory = 20) + ggtitle("Reactome: G1 doxy vs untreated")
ggsave("reactome_G1_dotplot_doxy_vs_untreated.png", p_G1, width = 7, height = 6, dpi = 300)
p_G2 <- dotplot(react_G2, showCategory = 20) + ggtitle("Reactome: G2 doxy vs untreated")
ggsave("reactome_G2_dotplot_doxy_vs_untreated.png", p_G2, width = 7, height = 6, dpi = 300)

## barplot
p_bar_G2 <- barplot(react_G2, showCategory = 20)
ggsave("reactome_G2_barplot_doxy_vs_untreated.png", p_bar_G2, width = 7, height = 6, dpi = 300)

## pathway gene network
p_cnet_G2 <- cnetplot(react_G2, showCategory = 10)
ggsave("reactome_G2_cnetplot.png", p_cnet_G2, width = 8, height = 7, dpi = 300)

## enrichment map
react_G2_sim <- enrichplot::pairwise_termsim(react_G2)
p_emap_G2 <- emapplot(react_G2_sim)
ggsave("reactome_G2_emapplot.png", p_emap_G2, width = 8, height = 7, dpi = 300)

### save result
write.csv(as.data.frame(react_G0), "reactome_G0.csv", row.names = FALSE)
write.csv(as.data.frame(react_G1), "reactome_G1.csv", row.names = FALSE)
write.csv(as.data.frame(react_G2), "reactome_G2.csv", row.names = FALSE)



#### ===================================
### Reactome
###get Reactome gene sets
reactome_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")  ## download the canonical pathways from Reactome
reactome_gs <- split(reactome_df$gene_symbol, reactome_df$gs_name)   ###convert it into a GSVA‑ready list
# ## Convert ENSEMBL → SYMBOL
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
 
   keep <- mapping$ENSEMBL %in% ens
   mapping <- mapping[keep, ]
 
   mat2 <- mat[match(mapping$ENSEMBL, ens), , drop = FALSE]
   rownames(mat2) <- mapping$SYMBOL
   mat2
 }
 
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
 
# ## apply
norm_G0_sym <- collapse_duplicates(convert_to_symbol(counts(dds_G0, normalized = TRUE)))
norm_G1_sym <- collapse_duplicates(convert_to_symbol(counts(dds_G1, normalized = TRUE)))
norm_G2_sym <- collapse_duplicates(convert_to_symbol(counts(dds_G2, normalized = TRUE)))

# 1. Ensure your input is a matrix (GSVA prefers matrix over data frames)
mat_G0 <- as.matrix(norm_G0_sym)
mat_G1 <- as.matrix(norm_G1_sym)
mat_G2 <- as.matrix(norm_G2_sym)

# 2. Define the parameter using the NEW syntax
# We use ssgseaParam and explicitly pass the gene sets list
# Note: kcdf = "Gaussian" is appropriate for log-transformed or VST data. 
# For raw normalized counts, some prefer "Poisson", but "Gaussian" is standard for ssGSEA.

param_G0 <- ssgseaParam(mat_G0, reactome_gs)
param_G1 <- ssgseaParam(mat_G1, reactome_gs)
param_G2 <- ssgseaParam(mat_G2, reactome_gs)

# 3. Run GSVA
gsva_G0 <- gsva(param_G0)
gsva_G1 <- gsva(param_G1)
gsva_G2 <- gsva(param_G2)

message("GSVA/ssGSEA successfully completed for all genotypes.")

# ==============================================================================
# APPLYING HEATMAPS 
# ==============================================================================
# Doxy-Only Progression (GL, G0, G1, G2) ---
# This shows how the 50 most significant genes change across all genotypes
plot_heatmap(dds_subset = dds_doxy, 
             res_annot  = res_G2_GL_doxy_annot, 
             title      = "Top 50 DEGs: Genotype Progression (Doxy)", 
             file       = "heatmap_doxy_progression_0.56.png")

# Treatment Effect within G2 ---
# This highlights the "Doxy-response" signature specifically for the G2 genotype
plot_heatmap(dds_subset = dds_G2, 
             res_annot  = res_G2_doxy_vs_untr_annot, 
             title      = "Top 50 DEGs: G2 Doxy vs Untreated", 
             file       = "heatmap_G2_treatment_effect_0.56.png")

# Treatment Effect within G1 ---
# This highlights the "Doxy-response" signature specifically for the G1 genotype
plot_heatmap(dds_subset = dds_G1, 
             res_annot  = res_G1_doxy_vs_untr_annot, 
             title      = "Top 50 DEGs: G1 Doxy vs Untreated", 
             file       = "heatmap_G1_treatment_effect_0.56.png")

# Treatment Effect within G0 ---
# This highlights the "Doxy-response" signature specifically for the G0 genotype
plot_heatmap(dds_subset = dds_G0, 
             res_annot  = res_G0_doxy_vs_untr_annot, 
             title      = "Top 50 DEGs: G0 Doxy vs Untreated", 
             file       = "heatmap_G0_treatment_effect_0.56.png")

# ==============================================================================
# REACTOME / GSVA INTEGRATION 
# ==============================================================================
# To visualize the GSVA results as a heatmap (Pathway Activity):

plot_gsva_heatmap <- function(gsva_results, title, file) {
  png(file, width = 1000, height = 800, res = 120)
  # Pick the top 30 pathways by variance to keep the plot readable
  pathway_vars <- apply(gsva_results, 1, var)
  top_pathways <- gsva_results[order(pathway_vars, decreasing = TRUE)[1:30], ]
  
  pheatmap(top_pathways, 
           main = title,
           scale = "row",
           color = colorRampPalette(c("blue", "white", "red"))(100),
           fontsize_row = 7)
  dev.off()
}

# Apply to GSVA outputs
plot_gsva_heatmap(gsva_G2, "Top Reactome Pathways: G2 Activity", "gsva_heatmap_G2.png")
plot_gsva_heatmap(gsva_G1, "Top Reactome Pathways: G1 Activity", "gsva_heatmap_G1.png")
plot_gsva_heatmap(gsva_G0, "Top Reactome Pathways: G0 Activity", "gsva_heatmap_G0.png")

