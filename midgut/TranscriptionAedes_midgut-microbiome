TranscriptionAedes_midgut-microbiome

# Processing the taxonomic classification, ran by Galaxy pipeline (Kraken2 + Bracken)
# Then processing and visualising in R

####################################################################################################
# Midgut microbiome composition from Bracken/Kraken output
#
# Language: R
#
# Aim:
#   Visualize bacterial composition in midgut PRO-seq libraries.
#
# Input:
#   Bracken total table files from Galaxy Kraken2/Bracken workflow.
#
# Final figure:
#   Show the per-replicate averages
#     Rep1 = mean(rep1_pool1, rep1_pool2)
#     Rep2 = mean(rep2_pool1, rep2_pool2)
#
####################################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# Choose parameters/settings 
taxonomy_level <- "G"      # G = genus
taxonomy_label <- "Genus"

RPM_threshold <- 500           # taxa below this RPM threshold are grouped as Other
scale_to_100 <- TRUE       # TRUE = relative abundance (bar plot scaling to 100%); FALSE = raw RPM

########################################################## 
# Input files
##########################################################

# Specify bracken file names
bracken_files <- c(
  rep1_pool1 = "./rep1_pool1_TotalTable.tabular",
  rep1_pool2 = "./rep1_pool2_TotalTable.tabular",
  rep2_pool1 = "./rep1_pool1_TotalTable.tabular",
  rep2_pool2 = "./rep2_pool2_TotalTable.tabular"
)

# For normalisation
seq_depths <- c(
  rep1_pool1 = 10.57e6,
  rep1_pool2 = 8.52e6,
  rep2_pool1 = 12.05e6,
  rep2_pool2 = 12.95e6
)


# Load in Bracken data, extract all bacteria, and normalise the counts
data_list <- list()

for (sample_name in names(bracken_files)) {

  bracken_data <- read.table(
    bracken_files[[sample_name]],
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )

  colnames(bracken_data) <- c(
    "percentage", "counts", "col3", "tax_level", "taxid", "name"
  )

  domain_indices <- which(bracken_data$tax_level == "D")

  first_domain_index <- domain_indices[1]

  if (length(domain_indices) > 1) {
    last_index <- domain_indices[2] - 1
  } else {
    last_index <- nrow(bracken_data)
  }

  bacteria_data <- bracken_data[first_domain_index:last_index, ]

  filtered_data <- bacteria_data %>%
    filter(tax_level == taxonomy_level) %>%
    select(name, counts) %>%
    mutate(
      name = trimws(name),
      RPM = counts / seq_depths[[sample_name]] * 1e6,
      sample = sample_name
    )

  data_list[[sample_name]] <- filtered_data
}

bracken_total <- bind_rows(data_list)

message("Samples included:")
print(unique(bracken_total$sample))

########################################################## 
# Make sample x genus table
##########################################################

all_samples <- names(bracken_files)

df_bracken_RPM <- bracken_total %>%
  group_by(sample, name) %>%
  summarise(RPM = sum(RPM, na.rm = TRUE), .groups = "drop")

df_bracken_relative <- df_bracken_RPM %>%
  complete(sample = all_samples, name, fill = list(RPM = 0))

# now have the 4 indivudal samples but want to pool/average per replicate

sample_to_rep <- c(
  rep1_pool1 = "Rep1",
  rep1_pool2 = "Rep1",
  rep2_pool1 = "Rep2",
  rep2_pool2 = "Rep2"
)

df_bracken_RPM_reps <- df_bracken_RPM %>%
  mutate(replicate = recode(sample, !!!sample_to_rep)) %>%
  group_by(replicate, name) %>%
  summarise(RPM = mean(RPM), .groups = "drop")

# Thresholding: above which RPM values to plot Genus name or instead group to Other
# Set above at 
	# RPM_threshold <- 500 

df_thr <- df_bracken_RPM_reps %>%
  mutate(name = if_else(RPM < RPM_threshold, "Other", name)) %>%
  group_by(replicate, name) %>%
  summarise(RPM = sum(RPM), .groups = "drop")

######################################################### 
# Converting to relative abundances/scaling to 100%
##########################################################

if (scale_to_100) {
  df_plot <- df_thr %>%
    group_by(replicate) %>%
    mutate(yval = RPM / sum(RPM, na.rm = TRUE)) %>%
    ungroup()

  y_lab <- "Relative abundance (%)"
  y_scale <- scale_y_continuous(labels = percent_format(accuracy = 1))

} else {
  df_plot <- df_thr %>%
    mutate(yval = RPM)

  y_lab <- "RPM"
  y_scale <- scale_y_continuous(labels = comma)
}

######################################################### 
# Plotting
##########################################################

order_vec <- c(
  "Elizabethkingia", "Riemerella",
  "Bacillus", "Staphylococcus",
  "Citrobacter", "Escherichia", "Klebsiella", "Metakosakonia",
  "Phytobacter", "Salmonella", "Stenotrophomonas", "Vibrio",
  "Delftia", "Other"
)

df_plot$name <- factor(df_plot$name, levels = order_vec, ordered = TRUE)
df_plot$replicate <- factor(df_plot$replicate, levels = c("Rep1", "Rep2"))

genus_cols <- c(
  Elizabethkingia = "#fa9fb5",
  Riemerella = "#dd3497",
  Bacillus = "#fc4e2a",
  Staphylococcus = "#bd0026",
  Citrobacter = "#edf8b1",
  Escherichia = "#c7e9b4",
  Klebsiella = "#7fcdbb",
  Metakosakonia = "#41b6c4",
  Phytobacter = "#1d91c0",
  Salmonella = "#225ea8",
  Stenotrophomonas = "blue3",
  Vibrio = "black",
  Delftia = "#41ab5d",
  Other = "grey80"
)

p <- ggplot(df_plot, aes(x = replicate, y = yval, fill = name)) +
  geom_bar(stat = "identity", position = "stack", width = 0.85) +
  y_scale +
  scale_fill_manual(values = genus_cols, drop = FALSE) +
  labs(
    title = "Midgut microbial composition",
    subtitle = paste0(
      "Mean per replicate; taxa with RPM < ",
      threshold,
      " grouped as Other"
    ),
    x = NULL,
    y = y_lab,
    fill = taxonomy_label
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )

print(p)
