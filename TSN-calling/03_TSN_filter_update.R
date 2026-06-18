# Exact TSN as based on PRO-seq 5 pnt coverage

# 3 steps
	# 1: defining TSS search regions (R script)
		# 01_TSS_regions.R
	# 2: calling 5pnt (bash)
		# 02_TSN_calls.sh
	# 3: filter and update TSS/TSN annotation (R script)
		# 03_TSN_filter_update.R


####################################################################################################
# 03_TSN_filter_update.R
#
# Language: R
#
# Aim:
#   Filter newly called TSN positions based on:
#     1. RPM coverage above the median in both replicates
#     2. Identical newly called TSN position in rep1 and rep2
#
# Output:
#   TSSlist_updated.txt
#
####################################################################################################

# Load libraries
library(dplyr)
library(tibble)
library(ggplot2)

# Input files
local_dir <- "../TSNcalling"
setwd(local_dir)

tss_plus_file  <- file.path(local_dir, "tss_region_plus.bed")
tss_minus_file <- file.path(local_dir, "tss_region_minus.bed")

rep1_file <- file.path(local_dir, "TSS_max_cov_and_pos_perregion_merged_rep1-all.txt")
rep2_file <- file.path(local_dir, "TSS_max_cov_and_pos_perregion_merged_rep2-all.txt")

outfile_updated_tss <- file.path(local_dir, "TSSlist_updated.txt")

# The TSS search regions
tss_pl <- read.table(tss_plus_file, header = FALSE)
tss_mn <- read.table(tss_minus_file, header = FALSE)

colnames(tss_pl) <- c("chr", "start", "end", "genename_isoform")
colnames(tss_mn) <- c("chr", "start", "end", "genename_isoform")

# The called positions from script 02_TSN_calls.sh
calledTSS_rep1 <- read.delim(rep1_file)
calledTSS_rep2 <- read.delim(rep2_file)
same_transcript_order <- all(calledTSS_rep1[, 6] == calledTSS_rep2[, 6])


# Build big dataframe containing for each current TSS also the search region, the newly called TSNs and their coverage

df_TSS <- data.frame(
  genename_isoform = calledTSS_rep1$genename.isoform,
  chr = calledTSS_rep1$chr,
  TSSregion_start = calledTSS_rep1$startsite.TSS,
  TSSregion_end = calledTSS_rep1$endsite.TSS,
  newTSS_rep1 = calledTSS_rep1$maxcoverage_position,
  newTSS_rep2 = calledTSS_rep2$maxcoverage_position,
  newTSScov_rep1 = calledTSS_rep1$maxcoverage_coverage,
  newTSScov_rep2 = calledTSS_rep2$maxcoverage_coverage
)

# Calculate distances
plus_df_TSS <- df_TSS %>%
  filter(genename_isoform %in% tss_pl$genename_isoform) %>%
  add_column(oldTSS = .$TSSregion_start + 100, .after = 2) %>%
  mutate(
    distTSS_rep1 = newTSS_rep1 - oldTSS,
    distTSS_rep2 = newTSS_rep2 - oldTSS,
    strand = "+"
  )

minus_df_TSS <- df_TSS %>%
  filter(genename_isoform %in% tss_mn$genename_isoform) %>%
  add_column(oldTSS = .$TSSregion_end - 100, .after = 2) %>%
  mutate(
    distTSS_rep1 = oldTSS - newTSS_rep1,
    distTSS_rep2 = oldTSS - newTSS_rep2,
    strand = "-"
  )

df_TSS <- bind_rows(plus_df_TSS, minus_df_TSS)


# For each called TSN also want to report normalized coverage
# Normalization factors based on wc -l of 5pnt bed files:
	# rep1: 32,708,222
	# rep2: 46,331,677

nf_rep1 <- 1 / 32708222 * 1000000
nf_rep2 <- 1 / 46331677 * 1000000

df_TSS <- df_TSS %>%
  mutate(
    newTSScov_rep1_RPM = newTSScov_rep1 * nf_rep1,
    newTSScov_rep2_RPM = newTSScov_rep2 * nf_rep2
  )

# Filter for coverage above median
median_rep1 <- median(df_TSS$newTSScov_rep1_RPM, na.rm = TRUE)
median_rep2 <- median(df_TSS$newTSScov_rep2_RPM, na.rm = TRUE)

# Only assign exact TSN for those sites where (1) coverage above median and (2) same position called in rep1 and rep2
df_TSS <- df_TSS %>%
  mutate(
    filter_rep1_RPM_median = newTSScov_rep1_RPM > median_rep1,
    filter_rep2_RPM_median = newTSScov_rep2_RPM > median_rep2,
    filter_same_RPM_median =
      filter_rep1_RPM_median == TRUE &
      filter_rep2_RPM_median == TRUE &
      newTSS_rep1 == newTSS_rep2
  )

message("TSSs passing median RPM filter and same-position filter: ", sum(df_TSS$filter_same_RPM_median, na.rm = TRUE))


#############################################################
# Dataframes for the total and only the changes TSS sites
#############################################################

# Changed TSSs only
# Only keep those sites where newly called TSS = TSN was assigned
changedTSS_median <- df_TSS %>%
  filter(filter_same_RPM_median == TRUE) %>%
  mutate(
    newTSS = newTSS_rep1,
    distTSS = distTSS_rep1
  )

# Total TSS list:
# If a transcript passes the filter, report the newly called TSS = TSN
# Otherwise keep the old annotated TSS
totalTSS_median <- df_TSS %>%
  mutate(
    newTSS = if_else(filter_same_RPM_median == TRUE, newTSS_rep1, as.integer(oldTSS)),
    distTSS = if_else(filter_same_RPM_median == TRUE, distTSS_rep1, 0)
  )

# Export
write.table(
  totalTSS_median,
  file = outfile_updated_tss,
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

# Plot distance between old and updated TSS
ggplot(changedTSS_median, aes(x = "", y = distTSS)) +
  stat_boxplot(geom = "errorbar", linewidth = 1.2) +
  geom_boxplot(
    linewidth = 1.2,
    fatten = 1.5,
    outlier.color = "grey",
    outlier.size = 0.4
  ) +
  theme_classic() +
  labs(
    x = "",
    y = "Distance between old and updated TSS",
    title = "Changed TSSs after median RPM filter"
  )

message("Median distance, changed TSSs: ", median(changedTSS_median$distTSS, na.rm = TRUE))
print(quantile(changedTSS_median$distTSS, na.rm = TRUE))

