# Exact TSN as based on PRO-seq 5 pnt coverage

# 3 steps
	# 1: defining TSS search regions (R script)
		# 01_TSS_regions.R
	# 2: calling 5pnt (bash)
		# 02_TSN_calls.sh
	# 3: filter and update TSS/TSN annotation (R script)
		# 03_TSN_filter_update.R



####################################################################################################
# 01_TSS_regions.R
#
# Language: R
#
# Aim:
#   Generate search regions to define TSNs from the VectorBase GFF annotation.
#
# Output:
#   tss_region_plus.bed
#   tss_region_minus.bed
#
# Notes:
#   Plus strand:  annotated TSS -100 to +400
#   Minus strand: annotated TSS -400 to +100, relative to strand direction
####################################################################################################

library(GenomicRanges)
library(rtracklayer)
library(tidyverse)


# Directories and input files

local_dir <- "../TSNcalling"
setwd(local_dir)

gff_file <- file.path(local_dir, "VectorBase-60_AaegyptiLVP_AGWG.gff")

plus_outfile  <- file.path(local_dir, "tss_region_plus.bed")
minus_outfile <- file.path(local_dir, "tss_region_minus.bed")


# Load and filter GFF file

gff_annot <- rtracklayer::import(gff_file) %>%
  as.data.frame() %>%
  tibble::as_tibble()

keep_biotypes <- c(
  "lncRNA", "ncRNA", "pre_miRNA", "protein_coding", "pseudogene",
  "RNase_MRP_RNA", "RNase_P_RNA", "rRNA", "snoRA", "snRNA",
  "SRP_RNA", "tRNA"
)

Aedes_refGene <- gff_annot %>%
  filter(gene_ebi_biotype %in% keep_biotypes) %>%
  select(seqnames, start, end, strand, gene, description, ID)

colnames(Aedes_refGene) <- c(
  "chr", "txStart", "txEnd", "strand",
  "geneName", "geneDescription", "txID"
)

# Remove genes too close to chromosome ends
Aedes_refGene <- Aedes_refGene %>%
  filter(txStart >= 500, txEnd >= 500)

# Standardize chromosome names: AaegL5_1 -> 1
Aedes_refGene$chr <- gsub("AaegL5_", "", as.character(Aedes_refGene$chr))


# Create TSS regions separately for plus and minus strand genes
Aedes_refGene_pl <- Aedes_refGene %>%
  filter(strand == "+") %>%
  mutate(
    TSS = txStart,
    TSSregion_start = TSS - 100,
    TSSregion_end = TSS + 400
  )

Aedes_refGene_mn <- Aedes_refGene %>%
  filter(strand == "-") %>%
  mutate(
    TSS = txEnd,
    TSSregion_start = TSS - 400,
    TSSregion_end = TSS + 100
  )

# Export as BED files

bed_pl <- data.frame(
  chr = Aedes_refGene_pl$chr,
  start = Aedes_refGene_pl$TSSregion_start,
  end = Aedes_refGene_pl$TSSregion_end,
  name = Aedes_refGene_pl$txID
)

bed_mn <- data.frame(
  chr = Aedes_refGene_mn$chr,
  start = Aedes_refGene_mn$TSSregion_start,
  end = Aedes_refGene_mn$TSSregion_end,
  name = Aedes_refGene_mn$txID
)

write.table(
  bed_pl,
  plus_outfile,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  bed_mn,
  minus_outfile,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

message("Plus-strand TSS regions written: ", nrow(bed_pl))
message("Minus-strand TSS regions written: ", nrow(bed_mn))
