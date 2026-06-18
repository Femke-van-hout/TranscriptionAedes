# TranscriptionAedes-core-promoter-motifs

# Steps
  # 01_prepare_promoter_regions.R
  # 02_search_motifs.sh
  # 03_plot_motifs.R

####################################################################################################
# 03_plot_motifs.R
#
# Language: R
#
# Aim:
#   Combine plus/minus motif search outputs and plot motif positions relative to old TSS and PRO-seq updated TSNs.
#
# Input:
#   Motif output files from 02_search_motifs.sh
#
# Output:
#   Histogram and density plots for INR, TATA, and DPE motif position, for either old/current TSS or PRO-seq updated TSN intervals
####################################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

# Directories
local_dir <- "../TSNcalling"
motif_dir <- file.path(local_dir, "core_promoter_motifs")
setwd(motif_dir)

# Function to read and combine motif results
read_motif_results <- function(old_plus_file,
                               old_minus_file,
                               new_plus_file,
                               new_minus_file,
                               motif_name,
                               position_correction = 0) {

  oldplus <- read.delim(old_plus_file, header = TRUE)
  oldminus <- read.delim(old_minus_file, header = TRUE)
  newplus <- read.delim(new_plus_file, header = TRUE)
  newminus <- read.delim(new_minus_file, header = TRUE)

  old <- bind_rows(oldplus, oldminus) %>%
    select(genename_isoform, motif_position, motif_sequence) %>%
    rename(
      pos_old = motif_position,
      motif_old = motif_sequence
    )

  new <- bind_rows(newplus, newminus) %>%
    select(genename_isoform, motif_position, motif_sequence) %>%
    rename(
      pos_new = motif_position,
      motif_new = motif_sequence
    )

  newold <- full_join(new, old, by = "genename_isoform") %>%
    mutate(
      motif = motif_name,
      pos_new = pos_new + position_correction,
      pos_old = pos_old + position_correction
    )

  return(newold)
}

# Function to plot motif positions
plot_motif_positions <- function(df, motif_name) {

  df_long <- df %>%
    pivot_longer(
      cols = c("pos_new", "pos_old"),
      names_to = "pos_type",
      values_to = "position"
    ) %>%
    filter(!is.na(position))

  p_hist <- ggplot(df_long, aes(x = position)) +
    geom_histogram(
      aes(color = pos_type, fill = pos_type),
      alpha = 0.4,
      position = "identity",
      binwidth = 1
    ) +
    theme_classic() +
    labs(
      title = paste0(motif_name, " motif positions"),
      x = "Motif position relative to TSS/TSN",
      y = "Count"
    )

  p_density <- ggplot(df_long, aes(x = position, color = pos_type)) +
    geom_density(linewidth = 1) +
    theme_classic() +
    labs(
      title = paste0(motif_name, " motif position density"),
      x = "Motif position relative to TSS/TSN",
      y = "Density"
    )

  print(p_hist)
  print(p_density)

  ggsave(
    filename = paste0("histogram_", motif_name, "_oldTSS_vs_newTSN.pdf"),
    plot = p_hist,
    width = 6,
    height = 4
  )

  ggsave(
    filename = paste0("density_", motif_name, "_oldTSS_vs_newTSN.pdf"),
    plot = p_density,
    width = 6,
    height = 4
  )
}

####################################################
# INR
####################################################

inr_results <- read_motif_results(
  old_plus_file = "motifs_INR_oldTSS_plus.txt",
  old_minus_file = "motifs_INR_oldTSS_minus.txt",
  new_plus_file = "motifs_INR_newTSS_plus.txt",
  new_minus_file = "motifs_INR_newTSS_minus.txt",
  motif_name = "INR",
  position_correction = 2
)

message("INR motifs in old TSS regions: ", sum(!is.na(inr_results$pos_old)))
message("INR motifs in new TSN regions: ", sum(!is.na(inr_results$pos_new)))

plot_motif_positions(inr_results, "INR")

####################################################
# TATA
####################################################

tata_results <- read_motif_results(
  old_plus_file = "motifs_TATA_oldTSS_plus.txt",
  old_minus_file = "motifs_TATA_oldTSS_minus.txt",
  new_plus_file = "motifs_TATA_newTSS_plus.txt",
  new_minus_file = "motifs_TATA_newTSS_minus.txt",
  motif_name = "TATA",
  position_correction = 0
)

message("TATA motifs in old TSS regions: ", sum(!is.na(tata_results$pos_old)))
message("TATA motifs in new TSN regions: ", sum(!is.na(tata_results$pos_new)))

plot_motif_positions(tata_results, "TATA")

####################################################
# DPE
####################################################

dpe_results <- read_motif_results(
  old_plus_file = "motifs_DPE_oldTSS_plus.txt",
  old_minus_file = "motifs_DPE_oldTSS_minus.txt",
  new_plus_file = "motifs_DPE_newTSS_plus.txt",
  new_minus_file = "motifs_DPE_newTSS_minus.txt",
  motif_name = "DPE",
  position_correction = 1
)

message("DPE motifs in old TSS regions: ", sum(!is.na(dpe_results$pos_old)))
message("DPE motifs in new TSN regions: ", sum(!is.na(dpe_results$pos_new)))

plot_motif_positions(dpe_results, "DPE")
