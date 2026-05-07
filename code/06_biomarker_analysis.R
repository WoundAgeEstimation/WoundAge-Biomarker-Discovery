#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro 
# @Editor: Annica Gosch
# @Description: Parse the metadata and perform exploratory analysis
# @software version: R=4.4.2/4.5.3

# Library
sh <- suppressPackageStartupMessages
sh(library(tidyverse))
library(data.table)
library(ComplexHeatmap)
library(circlize) 


#Set res 
out_data <- "results_new/06_biomarkers"
if (!dir.exists(out_data)){
  dir.create(out_data)
}

# Load the data----
# Load data (RPM). RPM is the reads normalized the total reads in the sample 
rpm1 <- fread("data/rpm/Auto_IvaGomes_THGE_29-05-2024_B_Chip1_torrent-server_6371.rpm.bcmatrix.txt") %>% column_to_rownames("Gene")
rpm2 <- fread("data/rpm/Auto_IvaGomes_THGE_25-06-2024_A_Chip2_torrent-server_21.rpm.bcmatrix.txt")%>% column_to_rownames("Gene")
rpm3 <- fread("data/rpm/Auto_IvaGomes_THGE_25-06-2024_B_Chip3_torrent-server_23.rpm.bcmatrix.xls.txt")%>% column_to_rownames("Gene")
rpm4 <- fread("data/rpm/Auto_IvaGomes_THGE_01-07-2024_D_Chip4_torrent-server_46.rpm.bcmatrix.txt")%>% column_to_rownames("Gene")

#Modify rpm4 names
colnames(rpm4)[2:5] <- paste0(colnames(rpm4)[2:5], "2")

a <- list(rpm1, rpm2, rpm3, rpm4)
rpms <- lapply(a, function(x) x %>% select(-c(Target, COSMIC_CGC_FLAG, NCBI_NAME, HGNC_SYMBOL_ACC, MIM_MORBID_DESCRIPTION, ENTREZ_GENE_ID, U133PLUS2_PSID)))
rpms <- do.call("cbind", rpms)

# Remove two duplicated samples 

rpms <- rpms %>% select(-c( IonXpress_037, IonXpress_038))

# Load matching samples to replace colnames 
match <- list(fread("data/expression_data/chip1_matching_samples.txt"), 
              fread("data/expression_data/chip2_matching_samples.txt"), 
              fread("data/expression_data/chip3_matching_samples.txt"), 
              fread("data/expression_data/chip4_matching_samples_mod.txt")
)

match <- read.csv("data/match.csv") %>% select(-Sample.Name)
match_vector <- match$Sample_name
names(match_vector) <- match$Barcode.ID

# Replace the names in the expression table
colnames(rpms) <- match_vector[colnames(rpms)]

# Load metadata
metadata_sample.full <- read.csv("data/metadata_wound_time_assigned.csv", row.names = 1, sep = ",") 
# load DEGs
degs_data <- readRDS(file.path("results_new", "02_degs", "DGA_filtered.rds"))

DEGs_control_late <- degs_data$`control vs prolonged`
DEGs_control_early <- degs_data$`control vs acute`
DEGs_control_intermediary <- degs_data$`control vs intermediate`

# Compute Ratio of expression ----
rpms.w <- rpms[,c(1:ncol(rpms)) %% 2 == 0]
rpms.c <- rpms[,c(1:ncol(rpms)) %% 2 == 1]

ratio.log <- log2((rpms.w + 0.1) / (rpms.c +0.1))
ratio <- (rpms.w + 0.1) / (rpms.c +0.1)

# Find biomarkers ----

#Differential expressed exclusively in that group
#Ratio >= 2 in all the sample (of that group). 
#Basal expression in the wound high or Basal expression in Wound high (positive and negative markers) 
#In all, I tried to get ~20 genes per group. 

degs_early <- DEGs_control_early %>% filter(FDR.P.val < 0.05)
degs_int <- DEGs_control_intermediary %>% filter(FDR.P.val < 0.05)
degs_late <- DEGs_control_late %>% filter(FDR.P.val < 0.05)


### Acute ----
#### Get exclusive genes
degs_early_exclusive <- setdiff(degs_early$ID, union(degs_int$ID, degs_late$ID))
degs_early_exclusive <- degs_early %>% filter(ID %in% degs_early_exclusive)

### Get ratio for those genes 
metadata_sample.full.early <- metadata_sample.full %>% filter(wound_time == "Acute")
early_samples <- metadata_sample.full.early %>% pull(Sample_name)
early_samples <- early_samples[grepl("^W", early_samples)]

ratio_early_genes <- ratio[degs_early_exclusive$ID,early_samples]

###Median ratios 
median_early_degs_ratio <-  apply(log(ratio_early_genes), 1, median)
plot(hist(median_early_degs_ratio, 100)) 

# Filter for genes with a high effect
degs_early_exclusive.high.effect <- degs_early_exclusive %>% filter(abs(FC) > 2)
rpms.w.early.genes.higheffect <- rpms.w[degs_early_exclusive.high.effect$ID, early_samples]
# Filter some lowexpressed genes
apply(rpms.w.early.genes.higheffect,1,median)[order(apply(rpms.w.early.genes.higheffect,1,median))]

median_rpm <- apply(rpms.w.early.genes.higheffect,1,median)
keep <- names(median_rpm[median_rpm > 10])

degs_early_exclusive.high.expression <- degs_early_exclusive.high.effect %>% filter(ID %in% keep)

# Ratio for these genes
ratio.early.candidate.markers <- ratio_early_genes[degs_early_exclusive.high.expression$ID,]
# order by median ratio
ratio.early.candidate.markers <- ratio.early.candidate.markers[order(apply(ratio.early.candidate.markers, 1, median), decreasing = T), ]

# Number of samples with ratio >2 
good_samples <- rowSums(ratio.early.candidate.markers > 2)

data_sorted <- sort(good_samples, decreasing = TRUE)

svg(file.path(out_data, "01_early_biomarkers_sample_number.svg"), w = 8, h = 6)
barplot(
  data_sorted,
  col = "skyblue",
  border = "black",
  las = 2, # rotate labels
  xlab = "",
  ylab = "Number of Samples with ratio > 2 (out of 11)"
)
dev.off()

#plot(good_samples, 
#        degs_early_exclusive.high.expression$FC)

# Biomarkers to plot
ratio.early.markers <- ratio.early.candidate.markers[good_samples >= 8,] # Select genes with ratio in more that 8 samples

### Intermediary ----
#### Get exclusive genes
degs_int_exclusive <- setdiff(degs_int$ID, union(degs_early$ID, degs_late$ID))
degs_int_exclusive <- degs_int %>% filter(ID %in% degs_int_exclusive)

### Get ratio for those genes 
#### I am mostly interested in positive markers, so no ploting in Log is fine here.
metadata_sample.full.int <- metadata_sample.full %>% filter(wound_time == "Intermediary")
int_samples <- metadata_sample.full.int %>% pull(Sample_name)
int_samples <- int_samples[grepl("^W", int_samples)]

ratio_int_genes <- ratio[degs_int_exclusive$ID,int_samples]

###Median ratios 
median_int_degs_ratio <-  apply(log(ratio_int_genes), 1, median)
plot(hist(median_int_degs_ratio, 100)) # Lot of genes with a high effect, without taking into account confounding factors. 

# Filter for genes with a high effect
degs_int_exclusive.high.effect <- degs_int_exclusive %>% filter(abs(FC) > 2)
# This time I got both positive and negative markers.
degs_int_exclusive.high.effect.w <- degs_int_exclusive.high.effect %>% filter(FC < 0)
rpms.w.int.genes.higheffect <- rpms.w[degs_int_exclusive.high.effect.w$ID, int_samples]

degs_int_exclusive.high.effect.c <- degs_int_exclusive.high.effect %>% filter(FC > 0)
rpms.c.int.genes.higheffect <- rpms.c[degs_int_exclusive.high.effect.c$ID, ]


# Filter some lowexpressed genes
#apply(rpms.w.int.genes.higheffect,1,median)[order(apply(rpms.w.int.genes.higheffect,1,median))]
#apply(rpms.c.int.genes.higheffect,1,median)[order(apply(rpms.w.int.genes.higheffect,1,median))]

median_rpm <- apply(rpms.w.int.genes.higheffect,1,median)
keep <- names(median_rpm[median_rpm > 10])
degs_int_exclusive.high.expression <- degs_int_exclusive.high.effect %>% filter(ID %in% keep)

median_rpm <- apply(rpms.c.int.genes.higheffect,1,median)
keep <- names(median_rpm[median_rpm > 10])
degs_int_exclusive.high.expression <- rbind(
  degs_int_exclusive.high.expression,
  degs_int_exclusive.high.effect %>% filter(ID %in% keep)
)

degs_int_exclusive.high.expression.candidate.markers <- degs_int_exclusive.high.expression %>% 
  arrange(-abs(FC)) %>%
  slice_head(n = 20)

# Ratio for these genes
ratio.int.candidate.markers <- ratio_int_genes[degs_int_exclusive.high.expression.candidate.markers$ID,]
# order by median ratio
ratio.int.candidate.markers <- ratio.int.candidate.markers[order(apply(ratio.int.candidate.markers, 1, median), decreasing = T), ]

# Number of samples with ratio >2 
good_samples.pos <- rowSums(ratio.int.candidate.markers > 2)
good_sample.neg <- rowSums(ratio.int.candidate.markers < 0.5)


svg(file.path(out_data, "01_intermediate_biomarkers_sample_number.svg"), w = 8, h = 6)
barplot(
  c(good_samples.pos, good_sample.neg),
  col = "skyblue",
  border = "black",
  las = 2, # rotate labels
  xlab = "",
  ylab = "Number of Samples with ratio > 2 (out of 6)"
)
dev.off()

good_samples.pos <- good_samples.pos[good_samples.pos >= 5]
good_sample.neg <- good_sample.neg[good_sample.neg >= 5]

good_samples <- c(good_sample.neg, good_samples.pos)


# Biomarkers to plot
ratio.int.markers <- ratio.int.candidate.markers[names(good_samples),] # Select genes with ratio in more that 5 samples

### Prolonged ----
#### Get exclusive genes
degs_late_exclusive <- setdiff(degs_late$ID, union(degs_int$ID, degs_early$ID))
degs_late_exclusive <- degs_late %>% filter(ID %in% degs_late_exclusive)

### Get ratio for those genes 
metadata_sample.full.late <- metadata_sample.full %>% filter(wound_time == "Prolonged")
late_samples <- metadata_sample.full.late %>% pull(Sample_name)
late_samples <- late_samples[grepl("^W", late_samples)]

ratio_late_genes <- ratio[degs_late_exclusive$ID,late_samples]

###Median ratios 
median_late_degs_ratio <-  apply(log(ratio_late_genes), 1, median)
plot(hist(median_late_degs_ratio, 100)) # Lot of genes with a high effect, without taking into account confounding factors. 

# Filter for genes with a high effect
degs_late_exclusive.high.effect <- degs_late_exclusive %>% filter(abs(FC) > 2)
rpms.w.late.genes.higheffect <- rpms.w[degs_late_exclusive.high.effect$ID, late_samples]

#Most of these are positive markers, only 1 is a negative marker.
#Since the expression of that gene is high seems fine, I do not need to filter for it in separate
# Filter some lowexpressed genes
apply(rpms.w.late.genes.higheffect,1,median)[order(apply(rpms.w.late.genes.higheffect,1,median))]

median_rpm <- apply(rpms.w.late.genes.higheffect,1,median)
keep <- names(median_rpm[median_rpm > 10])

degs_late_exclusive.high.expression <- degs_late_exclusive.high.effect %>% filter(ID %in% keep)

# Ratio for these genes
ratio.late.candidate.markers <- ratio_late_genes[degs_late_exclusive.high.expression$ID,]
# order by median ratio
ratio.late.candidate.markers <- ratio.late.candidate.markers[order(apply(ratio.late.candidate.markers, 1, median), decreasing = T), ]


# Number of samples with ratio >2 
good_samples <- rowSums(ratio.late.candidate.markers > 2)
#AMOTL2 is the exception, because its a negative marker
data_sorted <- sort(good_samples, decreasing = TRUE)
data_sorted[14] <- 3 # Visual inspection of the gene reveals 3 samples where the expression is higher in control with a ratio of 0.5 or lower. 

svg(file.path(out_data, "01_late_biomarkers_sample_number.svg"), w = 8, h = 6)
barplot(
  data_sorted,
  col = "skyblue",
  border = "black",
  las = 2, # rotate labels
  xlab = "",
  ylab = "Number of Samples with ratio > 2"
)
dev.off()

# Biomarkers to plot
ratio.late.markers <- ratio.late.candidate.markers[good_samples >= 5,] # Select genes with ratio in more that 6 samples

# Individual plots per genes ----
all_markers <- c(row.names(ratio.early.markers),
                 row.names(ratio.int.candidate.markers),
                 row.names(ratio.late.candidate.markers)
)

ratio.markers <- ratio[all_markers,] %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "ratio"
  )

# Merge with metadata
ratio.candidate.markers <- ratio.markers %>%
  merge(metadata_sample.full, by.x = "sample", by.y = "Sample_name")

wound_age_order <- c(
  "seconds",
  "seconds/ minutes",
  "minutes",
  "min",
  "15 min",
  "minutes - < 1 hour",
  "hours (< 10h)",
  "hours",
  "24 hours",
  "2 days",
  "55 hours (2,3 days)", 
  "70 hours (2,9 days)",
  "7 days",
  "12 days",
  "12/13 days", 
  "22 days", 
  "28 days"
)

ratio.candidate.markers.long <- ratio.candidate.markers %>% 
  mutate(`wound age` = factor(`wound.age`, levels = wound_age_order)) %>% 
  arrange(`wound age`) %>% 
  mutate(sample = factor(sample, levels = unique(sample)))

ratio.candidate.markers.long$ratio <- log2(ratio.candidate.markers.long$ratio)

plots <- lapply(unique(ratio.candidate.markers.long$gene), function(x){
  ratio.candidate.markers.long %>% 
    filter(gene == x) %>%
    ggplot(aes(x = sample, y = ratio)) +
    geom_point(size = 5, color = "steelblue") +
    geom_text(
      aes(label = paste(localisation, `wound age`, age, Gender, sep = ", ")),
      vjust = -0.5, size = 2
    ) +
    labs(
      title = paste("Expression for", x),
      x = "",
      y = "Log Expression Ratio"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "blue") +
    geom_hline(yintercept = -1, linetype = "dashed", color = "blue") +
    #expand_limits(y = max(ratio) * 1.1) + 
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(angle = 45, hjust = 1, size = 14),
      plot.margin = margin(t = 10, r = 40, b = 10, l = 20) 
    ) + 
    coord_cartesian(clip = "off")
})

for (i in seq_along(plots)) {
  ggsave(file.path(out_data, paste0("02_plot_", unique(ratio.candidate.markers.long$gene)[i], ".pdf")), plots[[i]], width = 12, height = 6)
}


# Final heatmap ----
ratio.candidate.markers.long <- ratio.candidate.markers.long %>% 
  arrange(`wound age`) 

sample_order <- unique(ratio.candidate.markers.long$sample)
genes<-c(rownames(ratio.early.markers),rownames(ratio.int.markers),rownames(ratio.late.markers))
ratio.candiates <- ratio[genes, ] %>% select(all_of(sample_order))

early_samples <- ratio.candidate.markers.long %>% 
  filter(wound_time == "Acute") %>% 
  pull(sample)

intermediate_samples <- ratio.candidate.markers.long %>% 
  filter(wound_time == "Intermediary") %>% 
  pull(sample)

late_samples <- ratio.candidate.markers.long %>% 
  filter(wound_time == "Prolonged") %>% 
  pull(sample) 


mat <- ratio.candiates
mat[mat < 1] <- -1 / mat[mat < 1]

mat[mat > 20] <- 20
mat[mat < -20] <- -20


stage_anno <- HeatmapAnnotation(
  Woundage = factor(
    c(rep("acute", length(unique(early_samples))),
      rep("intermediate", length(unique(intermediate_samples))),
      rep("prolonged", length(unique(late_samples)))),
    levels = c("acute", "intermediate", "prolonged")
  ),
  col = list(Woundage = c(
    "acute" = "#DDDDDD",
    "intermediate" = "#999999",
    "prolonged" = "#555555"
  ))
)

col_fun <- colorRamp2(
  c(min(mat), 0, max(mat)),
  c("royalblue", "white", "tomato1")
)



early_genes        <- row.names(ratio.early.markers) 
intermediary_genes <- row.names(ratio.int.markers)  
late_genes         <- row.names(ratio.late.markers) 


# Assign group to each gene in the same order as in mat
gene_stage <- ifelse(rownames(mat) %in% early_genes,        "acute",
                     ifelse(rownames(mat) %in% intermediary_genes, "intermediate",
                            ifelse(rownames(mat) %in% late_genes,         "prolonged", NA)))


row_anno <- rowAnnotation(
  Marker = factor(gene_stage, levels = c("acute", "intermediate", "prolonged")),
  col = list(Marker = c(
    "acute"        = "#9E9E9E", 
    "intermediate" = "#6E6E6E", 
    "prolonged"         = "#3B3B3B"
  ))
)

svg(file.path(out_data, "05_marker_heatmap.svg"), w = 8, h =10)
Heatmap(
  mat,
  name = "Ratio",
  top_annotation = stage_anno,
  left_annotation = row_anno,  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_order = sample_order,
  col = col_fun
)
dev.off()


# Save data 
## Data required for heatmap----

# --- Save all data required to rebuild the heatmap ---
heatmap_data <- list(
  mat = mat,
  sample_order = sample_order,
  early_samples = early_samples,
  intermediate_samples = intermediate_samples,
  late_samples = late_samples,
  early_genes = early_genes,
  intermediary_genes = intermediary_genes,
  late_genes = late_genes,
  gene_stage = gene_stage,
  col_fun = col_fun
)

saveRDS(heatmap_data, file = file.path(out_data, "heatmap_input_data.rds"))
write.csv(as.data.frame(mat), file = file.path(out_data, "heatmap_matrix.csv"), row.names = TRUE)

# Save data required for the gene ratio 
saveRDS(ratio.candidate.markers.long, file = file.path(out_data, "gene_ratio.rds"))
