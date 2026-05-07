#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro 
# @Editor: Annica Gosch
# @Description: Parse the metadata and perform exploratory analysis
# @software version: R=4.4.2/4.5.3

# Load libraries ----
spsm <- suppressPackageStartupMessages
spsm(library(tidyverse))
spsm(library(data.table))
spsm(library(reshape2))
spsm(library(ggcorrplot))
spsm(library(ComplexHeatmap))
spsm(library(circlize))
spsm(library(RColorBrewer))
spsm(library(caret))
spsm(library(lme4))

# Create output dir ----
out_data <- "results_new/01_data_exploration"
if (!dir.exists(out_data)){
  dir.create(out_data, recursive  =T )
}

# Load data ----
# (RPM). RPM is the reads normalized the total reads in the sample
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
metadata_sample.full <- read.csv("data/metadata.csv", row.names = 1, sep = ";") 
metadata_sample.full <- metadata_sample.full[-c(21,22),]


# Basic transformation and selection ----

# location mapping
location_mapping <- function(raw_locations) {
  normalized <- tolower(trimws(raw_locations))
  dplyr::case_when(
    grepl("chest|thoraco|upper chest", normalized) ~ "Chest",
    grepl("abdomen|body center|flank|groin|hip|pelvis", normalized) ~ "Abdomen/Pelvis",
    grepl("arm|hand", normalized) ~ "Upper Limb",
    grepl("leg|thigh|tight|femoral", normalized) ~ "Lower Limb",
    grepl("neck", normalized) ~ "Neck",
    grepl("\\?|or", normalized) ~ "Ambiguous",
    TRUE ~ "Unknown"
  )
}

metadata_sample.full$body_location <- location_mapping(metadata_sample.full$localisation)
metadata_sample.full$Chip <- as.character(metadata_sample.full$Chip)

metadata_sample.full.selected <- metadata_sample.full %>% 
  mutate(group = ifelse(substr(metadata_sample.full$Sample, 1, 1) == "R", "Control", "Wound")) %>%
  select(c(Subject, Sample_name, Age, Gender, RIN, Chip, PMI, body_location, case.type,group))

categorical_vars <- metadata_sample.full.selected %>%
  select(where(~ is.character(.x) || is.factor(.x))) %>%
  names()
categorical_vars <- categorical_vars[!categorical_vars %in% c("Subject", "Sample_name")]

# Replace empty strings by NA
metadata_sample.full.selected <- metadata_sample.full.selected %>% as.data.frame() %>%
  mutate(across(where(~ is.character(.) | is.factor(.)), ~ na_if(.x, "")))

encoded_df <- metadata_sample.full.selected

# One-hot encode categorical variables with NAs preserved
for (var in categorical_vars) {
  dummies <- dummyVars(reformulate(var), data = metadata_sample.full.selected)
  one_hot_data <- predict(dummies, newdata = metadata_sample.full.selected)
  colnames(one_hot_data) <- gsub("_", " ", colnames(one_hot_data))
  
  # Mark rows that were all-zero => original NA; replace with NA rows
  NA_rows <- rowSums(one_hot_data) == 0
  if (any(NA_rows)) one_hot_data[NA_rows, 1:ncol(one_hot_data)] <- NA
  
  if (ncol(one_hot_data) == 2){
    one_hot_data <- one_hot_data[,1, drop = FALSE]
  }
  
  encoded_df <- encoded_df %>%
    select(-all_of(var)) %>%
    bind_cols(as.data.frame(one_hot_data))
}

# Unsupervised analysis ----
## PCA (Wound samples only)
# Log transform full rpm matrix 
rpms_log <- log(rpms + 1)

# Subset to wound samples only
encoded_df.wound <- encoded_df[encoded_df$groupControl == "0", ]
encoded_df.wound <- encoded_df.wound[encoded_df.wound$Sample %in% colnames(rpms_log), ]

# Prepare rpms for wound samples (log version)
rpms.log.wound <- rpms_log %>% select(all_of(encoded_df.wound$Sample))

## Keep most variable genes across wound samples
sd_genes_wound <- apply(rpms.log.wound, 1, sd)
top5000_var_genes_wound <- order(sd_genes_wound, decreasing = TRUE)[1:5000]
rpms_wound_top <- rpms.log.wound[top5000_var_genes_wound, ] %>% t()
rpms_wound_top <- scale(rpms_wound_top, center = TRUE, scale = TRUE)

# Run PCA
pca_data <- prcomp(rpms_wound_top)
pc_scores <- pca_data$x[, 1:4]
pc_var <- pca_data$sdev^2
pc_var_explained <- round(100 * pc_var / sum(pc_var), 1)  # % variance
pc_names <- paste0("PC", 1:4, " (", pc_var_explained[1:4], "%)")
colnames(pc_scores) <- pc_names
rownames(pc_scores) <- rownames(rpms_wound_top)

metadata_aligned <-  encoded_df.wound[match(rownames(pc_scores), encoded_df.wound$Sample_name), ]

# Align metadata to PCA scores
pc_meta <- cbind(pc_scores, metadata_aligned)

# remove columns not intended to be used in PCA
metadata_aligned <- metadata_aligned %>% 
  select(-c(Subject, Sample_name, `body locationChest`, `body locationLower Limb`, `body locationNeck`, `body locationUpper Limb`, `body locationAbdomen/Pelvis`,`groupControl`))

# Rename metadata columns (hardcoded as before)
colnames(metadata_aligned) <- c("Age", "RIN", "PMI", "Sex (Female)", "Chip 01", "Chip 02", "Chip 03", "Chip 04", "Case type (Criminal)")
metadata_aligned <- metadata_aligned[,c(1,4,2,3,9,5:8)]

# compute correlations and p-values
cor_matrix <- matrix(NA, nrow = ncol(pc_scores), ncol = ncol(metadata_aligned))
rownames(cor_matrix) <- colnames(pc_scores)
colnames(cor_matrix) <- colnames(metadata_aligned)

cor_pvalues_matrix <- matrix(NA, nrow = ncol(pc_scores), ncol = ncol(metadata_aligned))
rownames(cor_pvalues_matrix) <- colnames(pc_scores)
colnames(cor_pvalues_matrix) <- colnames(metadata_aligned)

for (i in 1:ncol(pc_scores)) {
  for (j in 1:ncol(metadata_aligned)) {
    x <- pc_scores[,i]
    y <- metadata_aligned[,j]
    cor_matrix[i, j] <- cor(x, y, use = "complete.obs", method = "spearman")
    cor_pvalues_matrix[i,j] <- cor.test(x, y, use = "complete.obs", method = "spearman")$p.value
  }
}

# adjust p-values across all tests (FDR)
cor_pvalues_adj_matrix <- matrix(
  p.adjust(as.vector(cor_pvalues_matrix), method = "fdr"),
  nrow = nrow(cor_pvalues_matrix),
  ncol = ncol(cor_pvalues_matrix)
)
rownames(cor_pvalues_adj_matrix) <- rownames(cor_pvalues_matrix)
colnames(cor_pvalues_adj_matrix) <- colnames(cor_pvalues_matrix)

# create data.frame for plotting
cor_df <- melt(cor_matrix, varnames = c("PC", "Metadata"), value.name = "cor")
pvalue_df <- melt(cor_pvalues_adj_matrix, varnames = c("PC", "Metadata"), value.name = "FDR")
cor_df$FDR <- pvalue_df$FDR

cor_df$Metadata <- factor(cor_df$Metadata, levels = c("Age", "Sex (Female)","RIN", "PMI",  "Case type (Criminal)", "Chip 01", "Chip 02", "Chip 03", "Chip 04"))

svg(file.path(out_data, "01_PC_vs_metadata_Wound_final.svg"), w = 10, h = 8)
ggplot(cor_df, aes(x = Metadata, y = PC, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", cor)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "Correlation / R", limits = c(-1, 1)) +
  theme_minimal() +
  coord_fixed() +
  scale_y_discrete(limits = rev) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Association between First 4 PCs and Metadata (Wound Samples)") + 
  xlab("") + ylab("")
dev.off()

## Repeat PCA (Reference samples only)
# Log transform full rpm matrix 
rpms_log <- log(rpms + 1)

# Subset to reference samples only
encoded_df.Ref <- encoded_df[encoded_df$groupControl == "1", ]
encoded_df.Ref <- encoded_df.Ref[encoded_df.Ref$Sample %in% colnames(rpms_log), ]

# Prepare rpms for reference samples (log version)
rpms.log.Ref <- rpms_log %>% select(all_of(encoded_df.Ref$Sample))

## Keep most variable genes across wound samples
sd_genes_Ref <- apply(rpms.log.Ref, 1, sd)
top5000_var_genes_Ref <- order(sd_genes_Ref, decreasing = TRUE)[1:5000]
rpms_Ref_top <- rpms.log.Ref[top5000_var_genes_Ref, ] %>% t()
rpms_Ref_top <- scale(rpms_Ref_top, center = TRUE, scale = TRUE)

# Run PCA
pca_data <- prcomp(rpms_Ref_top)
pc_scores <- pca_data$x[, 1:10]
pc_var <- pca_data$sdev^2
pc_var_explained <- round(100 * pc_var / sum(pc_var), 1)  # % variance
pc_names <- paste0("PC", 1:10, " (", pc_var_explained[1:10], "%)")
colnames(pc_scores) <- pc_names
rownames(pc_scores) <- rownames(rpms_Ref_top)

metadata_aligned <-  encoded_df.Ref[match(rownames(pc_scores), encoded_df.Ref$Sample_name), ]

# Align metadata to PCA scores
pc_meta <- cbind(pc_scores, metadata_aligned)

# remove columns not intended to be used in PCA
metadata_aligned <- metadata_aligned %>% 
  select(-c(Subject, Sample_name, `body locationChest`, `body locationLower Limb`, `body locationNeck`, `body locationUpper Limb`, `body locationAbdomen/Pelvis`,`groupControl`))

# Rename metadata columns (hardcoded as before)
colnames(metadata_aligned) <- c("Age", "RIN", "PMI", "Sex (Female)", "Chip 01", "Chip 02", "Chip 03", "Chip 04", "Case type (Criminal)")

# compute correlations and p-values
cor_matrix <- matrix(NA, nrow = ncol(pc_scores), ncol = ncol(metadata_aligned))
rownames(cor_matrix) <- colnames(pc_scores)
colnames(cor_matrix) <- colnames(metadata_aligned)

cor_pvalues_matrix <- matrix(NA, nrow = ncol(pc_scores), ncol = ncol(metadata_aligned))
rownames(cor_pvalues_matrix) <- colnames(pc_scores)
colnames(cor_pvalues_matrix) <- colnames(metadata_aligned)

for (i in 1:ncol(pc_scores)) {
  for (j in 1:ncol(metadata_aligned)) {
    x <- pc_scores[,i]
    y <- metadata_aligned[,j]
    cor_matrix[i, j] <- cor(x, y, use = "complete.obs", method = "spearman")
    cor_pvalues_matrix[i,j] <- cor.test(x, y, use = "complete.obs", method = "spearman")$p.value
  }
}

# adjust p-values across all tests (FDR)
cor_pvalues_adj_matrix <- matrix(
  p.adjust(as.vector(cor_pvalues_matrix), method = "fdr"),
  nrow = nrow(cor_pvalues_matrix),
  ncol = ncol(cor_pvalues_matrix)
)
rownames(cor_pvalues_adj_matrix) <- rownames(cor_pvalues_matrix)
colnames(cor_pvalues_adj_matrix) <- colnames(cor_pvalues_matrix)

# create data.frame for plotting
cor_df <- melt(cor_matrix, varnames = c("PC", "Metadata"), value.name = "cor")
pvalue_df <- melt(cor_pvalues_adj_matrix, varnames = c("PC", "Metadata"), value.name = "FDR")
cor_df$FDR <- pvalue_df$FDR

cor_df$Metadata <- factor(cor_df$Metadata, levels = c("Age", "RIN", "PMI", "Sex (Female)", "Chip 01", "Chip 02", "Chip 03", "Chip 04", "Case type (Criminal)", "Institute A", "Institute B", "Institute C"))

svg(file.path(out_data, "01_PC_vs_metadata_Reference_Final.svg"), w = 10, h = 8)
ggplot(cor_df, aes(x = Metadata, y = PC, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", cor)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "Correlation / R", limits = c(-1, 1)) +
  theme_minimal() +
  coord_fixed() +
  scale_y_discrete(limits = rev) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Association between First 10 PCs and Metadata (Reference Samples)") + 
  xlab("") + ylab("")
dev.off()


## Hierarchical clustering (Wound samples only)
# Log transform full rpm matrix
rpms_log <- log(rpms + 1)

# Subset to wound samples only
encoded_df.wound <- encoded_df[encoded_df$groupControl == "0", ]
encoded_df.wound <- encoded_df.wound[encoded_df.wound$Sample %in% colnames(rpms_log), ]

# Prepare rpms for wound samples (log version)
rpms.log.wound <- rpms_log %>% select(all_of(encoded_df.wound$Sample))

# Create a wound age factor, dividing the current samples
metadata_sample.wound_age <- metadata_sample.full %>% select(Subject,  wound.age) %>% distinct()
metadata_sample.wound_age <- metadata_sample.wound_age %>% rename(wound_age = `wound.age`)

metadata_sample.wound_age$wound_age <- c("2-3 days", "2-3 days",
                                         "7 days", "12-13 days", 
                                         "minutes","minutes", 
                                         "2-3 days", "hours", 
                                         "1 day", ">20 days", 
                                         "12-13 days", "2-3 days", 
                                         "minutes", "minutes", 
                                         "hours", "minutes", 
                                         ">20 days","minutes", 
                                         "minutes", "minutes", 
                                         "minutes","2-3 days", 
                                         "7 days")

metadata_sample.wound_age$wound_age <- factor(
  metadata_sample.wound_age$wound_age,
  levels = c("minutes", "hours", "1 day", "2-3 days", "7 days", "12-13 days", ">20 days"),
  ordered = TRUE
)

# Merge wound_age into sample-level metadata for wound samples
metadata_sample.full.selected.wound <- metadata_sample.full.selected %>%
  filter(group == "Wound") %>%
  left_join(metadata_sample.wound_age, by = "Subject")


## Hierarchical clustering  ----
## Remove the effect of other variables with a linear regression

# Compute residuals per gene using lmer (preserve random effect Chip)
residuals_mat <- apply(rpms.log.wound, 1, function(x){
  all_data <- metadata_sample.full.selected.wound
  all_data$expression <- x
  resid(lmer(expression ~ Age + Gender + RIN + PMI + case.type + (1|Chip), data = all_data))
})
residuals_mat <- t(residuals_mat)
saveRDS(residuals_mat, file = file.path(out_data, "residuals_model.rds"))


# select top 5000 variable genes from residuals
sd_genes_res <- apply(residuals_mat, 1, sd)
top5000_var_genes_res <- order(sd_genes_res, decreasing = TRUE)[1:5000]
residuals.top <- residuals_mat[top5000_var_genes_res, ] %>% t()
residuals.top.scaled <- scale(residuals.top, center = TRUE, scale = TRUE)
row.names(residuals.top) <- colnames(rpms.log.wound)

dist_mat_res <- dist(residuals.top.scaled)
hc_res <- hclust(dist_mat_res, method = "ward.D2")

# Colors / annotation for heatmap 
wound_time_cols <- brewer.pal(n = 7, "Greys")[c(4,5,6,1,2,3,7)]

names(wound_time_cols) <- unique(metadata_sample.full.selected.wound$wound_age)

ha_res <- HeatmapAnnotation(
  `Wound age` = metadata_sample.full.selected.wound$wound_age,
  col = list(
    `Wound age` = wound_time_cols
  ),
  annotation_name_side = "left"
)

row.names(residuals.top.scaled) <- substr(row.names(residuals.top), 2, 4)

ht_res <- Heatmap(
  t(residuals.top.scaled),
  cluster_columns = as.dendrogram(hc_res),
  cluster_rows = FALSE,
  show_row_names = FALSE,
  show_column_names = T,
  top_annotation = ha_res,
  column_title = "HC on residuals (top 5000 variable genes in wound samples)",
  heatmap_legend_param = list(title = "Residual expression"),
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
)

svg(file.path(out_data, "02_HC_residuals_top5000_genes_final.svg"), w = 12)
draw(ht_res)
dev.off()

# Save residuals-based HC data
saveRDS(list(
  residuals = residuals_mat,
  residuals_top = residuals.top,
  residuals_top_scaled = residuals.top.scaled,
  hc_res = hc_res
), file = file.path(out_data, "hc_residuals_inputs_outputs.rds"))
