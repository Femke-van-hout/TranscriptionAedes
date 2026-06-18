TranscriptionAedes_midgut-expression

# Compare expression of previously described midgut markers in the Aag2 vs midgut nascent RNA libraries

####################################################################################################
# Compare enterocyte marker expression in Aag2 versus midgut PRO-seq libraries
#
# Aim:
#   Compare nascent RNA signal for previously described midgut/enterocyte marker genes
#   between Aag2 and midgut PRO-seq libraries.
#
# Input:
#   Enterocytes.tsv
#   genecounts_Aag2_midgut
#
# Output:
#   Heatmap of top 100 enterocyte marker genes ordered by Mosquito Cell Atlas logFC.
####################################################################################################


# Load libraries
library(data.table)
library(dplyr)
library(tibble)
library(pheatmap)

# Genes of interest
# From Mosquito Single Cell Atlas: https://cells.ucsc.edu/?ds=mosquito and then gut
# Then extract the top-100 enterocyte genes by logFC 
ent <- fread("Enterocytes.tsv")
ent$logFC <- ent$`logFC|float`
ent$gene <- ent$id
top_genes <- ent %>%
  arrange(desc(logFC)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  slice_head(n = 100) %>%
  pull(gene)


# Filter the Aag2 and midgut PRO-seq 3pnt coverage dataframe for genes of interest
# Coverage pooled per replicate (only PBS/unstimulated samples)

# Filter for top_genes and convert to numeric matrix
mat_pooled <- genecounts_Aag2_midgut %>%
  filter(gene %in% top_genes) %>%
  distinct(gene, .keep_all = TRUE) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# small pseudocount to avoid zeros before z-scoring
mat_pooled <- mat_pooled + 1e-6

# z-score calculations
# for the matrix with per-replicate pools
zmat_pool <- t(scale(t(mat_pooled)))
zmat_pool[is.na(zmat_pool)] <- 0

# annotate columns: Aag2 and midgut data each 2 replicates
cond_pool <- data.frame(
  group = factor(c("Aag2","Aag2","Midgut","Midgut"), levels = c("Aag2","Midgut"))
)
row.names(cond_pool) <- colnames(zmat_pool)

# Ensure that heatmap is plotted with gene order according to descending logFC from the Mosquito Cell Atlas
gene_order <- top_genes  
zmat_pool_ordered <- zmat_pool[gene_order, , drop = FALSE]


# Plot heatmap
pheatmap(
  zmat_pool_ordered,
  cluster_rows = FALSE,  # turned off row clustering to keep to gene order
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = cond_pool,
  color = colorRampPalette(c("#6059B6","#E5E7E9", "#CB4335"))(101),
  main = "Top 100 Enterocyte genes - ordered by logFC"
)


