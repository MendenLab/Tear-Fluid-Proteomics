library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggrepel)
library(ggpubr)

# protein data
data_proteins = X20210125_142432_Slot1_45_Report
data_genes_proteins = data_proteins[,c(2,3,6)]
data_proteins_original = data_proteins
data_proteins = data_proteins[,c(7:216)]
rownames(data_proteins) = data_genes_proteins$PG.ProteinGroups
colnames = colnames(data_proteins)
colnames_cleaned = gsub("^\\[[0-9]+\\]\\s*", "", colnames)       # Remove [xxx] at start
colnames_cleaned = gsub("\\.PG\\.Quantity$", "", colnames_cleaned) 
colnames(data_proteins) = colnames_cleaned

# metadata
metadata =  Slot1_45_IdentificationsOverview_all_samples_Ana
metadata = metadata[!is.na(metadata$RunNr),]
metadata = metadata[,c(2,4)]

# analysis
data_proteins = data_proteins[,metadata$FileName]

# change "Filtered" to NA
non_missing_counts = rowSums(data_proteins != "Filtered")

completeness_df <- data.frame(
  protein_id = rownames(data_proteins),
  num_quantified_replicates = non_missing_counts)

# Create histogram using ggplot2
completeness_plot = ggplot(completeness_df, aes(x = num_quantified_replicates)) +
  geom_histogram(binwidth = 1, fill = "#3182bd", color = "white", boundary = 0) +
  scale_x_continuous(breaks = seq(0, ncol(data_proteins), by = 50)) +
  labs(
    title = "Protein Data Completeness",
    x = "Number of samples with quantified replicates",
    y = "Number of proteins"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.title.position = "plot",
    panel.grid.minor = element_blank()
  )
pdf("plot_completeness.pdf",width = 8)
completeness_plot
dev.off()

# meatadata data cleaning
extract_dia_marker <- function(x) {
  # Remove index like [41] and quotes
  x <- gsub("\\[\\d+\\]\\s*\"", "", x)
  x <- gsub("\"", "", x)
  
  # Try to match pattern: DIA_X_YY or DIA_YY
  m <- regexec("DIA(?:_[A-Z]+)?_(\\d+)", x)
  reg_result <- regmatches(x, m)
  
  # Extract number and return as DIA_XX format
  if (length(reg_result[[1]]) >= 2) {
    return(paste0("DIA_", reg_result[[1]][2]))
  } else {
    return(NA)
  }
}
metadata$sample = sapply(metadata$FileName,extract_dia_marker)
metadata = metadata %>% rename(replicate_id = FileName)

# coefficient of variation
# Step 1: Make long-format expression data
data_proteins_long <- data_proteins %>%
  as.data.frame() %>%
  mutate(protein = data_genes_proteins$PG.ProteinGroups) %>%
  pivot_longer(
    cols = -protein,
    names_to = "replicate_id",
    values_to = "intensity"
  )

# Step 2: Join with metadata
data_proteins_annotated <- data_proteins_long %>%
  left_join(metadata, by = "replicate_id")

# Step 3: Calculate CV per protein per sample
cv_per_sample <- data_proteins_annotated %>%
  group_by(protein, sample) %>%
  mutate(intensity = ifelse(intensity == "Filtered",NA,as.numeric(intensity))) %>%
  filter(sum(!is.na(intensity)) >= 2) %>%  # Only keep if ≥ 2 replicates
  filter(!is.na(intensity)) %>%
  summarise(
    mean_intensity = mean(intensity, na.rm = TRUE),
    sd_intensity = sd(intensity, na.rm = TRUE),
    cv = (sd_intensity / mean_intensity) * 100,
    .groups = "drop"
  )

# mean per protein
mean_cv_per_protein <- cv_per_sample %>%
  group_by(protein) %>%
  summarise(mean_cv = mean(cv, na.rm = TRUE)) %>%
  arrange(desc(mean_cv))  # Sorted by variability

# integrate protein names
mean_cv_per_protein <- mean_cv_per_protein %>%
  left_join(data_proteins_original %>% 
              select(PG.ProteinGroups,PG.Genes,PG.ProteinNames) %>%
              rename(protein = PG.ProteinGroups))

proteins_interest = c("CRYM","PFKL", "CAPZA2", "ALDH16A1", "SERPINC1", "HP")

mean_cv_per_protein <- mean_cv_per_protein %>%
  mutate(highlight = ifelse(PG.Genes %in% proteins_interest, "Protein signature", "Other"))

# order proteins
mean_cv_per_protein <- mean_cv_per_protein %>%
  arrange(desc(mean_cv)) %>%
  mutate(Protein = factor(protein, levels = protein))

# Create the scatter plot
plot_CV = ggplot(mean_cv_per_protein, aes(x = Protein, y = mean_cv, color = highlight)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Protein signature" = "#d89027", "Other" = "#3182bd")) +
  geom_text_repel(
    data = mean_cv_per_protein %>% filter(highlight == "Protein signature"),
    aes(label = PG.Genes),
    #color = "black",
    size = 3.5,
    box.padding = 0.8,
    point.padding = 0.3,
    segment.color = 'grey50',
    max.overlaps = 20,
    fontface = "bold"
  ) +
  # geom_text(data = mean_cv_per_protein %>% filter(highlight == "Protein signature"),
  #           aes(label = PG.Genes), vjust = -1, size = 3.5, color = "black") +
  labs(
    title = "Coefficient of variation (CV) per protein",
    x = "Protein",
    y = "CV (%)",
    color = "Protein Type"
  ) +
  theme_minimal(base_size = 12.5) +
  theme(
    axis.text.x = element_blank(),  # Optional: hide if too many proteins
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.ticks.x = element_line(color = "gray70"),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) 


pdf("plot_CV.pdf",width = 10)
plot_CV
dev.off()

plot_all = ggarrange(completeness_plot, plot_CV, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)

pdf("plot_completeness_CV.pdf",width = 10, height = 8)
plot_all
dev.off()