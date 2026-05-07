#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro 
# @Editor: Annica Gosch
# @Description: Parse the metadata and perform exploratory analysis
# @software version: R=4.4.2/4.5.3

# Load libraries ----
shhh <- suppressPackageStartupMessages
shhh(library(tidyverse))
shhh(library(ggrepel))
shhh(library(UpSetR))

#Note: On the Tables the down genes are upregulated on the Wound vs control. 
#On the plot its the opposite since its easier to interpret

# Set directory
out.dir <- "results_new/02_degs"

# SET FDR threshold
fdr_treshold <- 0.05
degs_early <- read.csv("results_new/02_degs/control_vs_acute.txt", sep = "\t", dec = ",", skip = 4) 
degs_intermediary <- read.csv("results_new/02_degs/control_vs_intermediate.txt", sep = "\t", dec = ",", skip = 4)
degs_late <- read.csv("results_new/02_degs/control_vs_prolonged.txt", sep = "\t", dec = ",", skip = 4) 
# Between timepoint comparison
degs_early_intermediary <- read.csv("results_new/02_degs/acute_vs_intermediate.txt", sep = "\t", dec = ",", skip = 4) 
degs_intermediary_late <- read.csv("results_new/02_degs/intermediate_vs_prolonged.txt", sep = "\t", dec = ",", skip = 4) 

DGA <- list(
  `control vs acute` = degs_early,
  `control vs intermediate` = degs_intermediary,
  `control vs prolonged` = degs_late,
  `acute vs intermediate` = degs_early_intermediary,
  `intermediate vs prolonged` = degs_intermediary_late
)

# Remove lowly expressed genes
DGA <- lapply(
  DGA, function(x) {
    remove <- x[,2] <= 1 | x[,3] <= 1
    x <- x[!remove, ]
    print(sum(remove))
    # Adjust p-values (FDR)
    x$FDR.P.val <- p.adjust(x$P.val, method = "fdr")
    
    return(x)
  }
)

# Reverse the sign of logFC
DGA <- lapply(DGA, function(x) x %>% mutate(Fold.Change = -Fold.Change))

# get DEGs 
DEGS <- lapply(DGA, function(x) x %>% filter(FDR.P.val < fdr_treshold))

number_of_degs <- lapply(names(DEGS), function(x) {
  degs <- DEGS[[x]]
  down <- nrow(degs %>% filter(Fold.Change < 0 ))
  up <- nrow(degs %>% filter(Fold.Change > 0 ))
  
  data.frame("Comparison" = x, "Downregulated" = up, "Upregulated" = down)
} 
)

number_of_degs <- do.call("rbind", number_of_degs)

degs_summary_long <- number_of_degs %>%
  pivot_longer(cols = c(Upregulated, Downregulated), names_to = "Regulation", values_to = "Count")

degs_summary_long$Regulation <- factor(degs_summary_long$Regulation, levels = c("Upregulated", "Downregulated"))
degs_summary_long$Comparison <- factor(degs_summary_long$Comparison, levels = c(
  "control vs acute",
  "control vs intermediate",
  "control vs prolonged",
  "acute vs intermediate",
  "intermediate vs prolonged"
))

# Plot number of DEGs ----
a <- ggplot(degs_summary_long, aes(x = Comparison, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(x = "", y = "Number of DEGs", title = "", fill = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
  ) +
  scale_fill_manual(values = c("Upregulated" = "#CC6677", "Downregulated" = "#88CCEE"))  

svg(paste0(out.dir, "/01_degs.svg"), w = 8, h = 5)
print(a)
dev.off()

# Overlap between comparison----

#Not taken into account directions here, because early vs intermediary 
# upregulated might be down in early and it would be meaningful. 
# Although, we have to take care in control_vs intermediary, as that comparison has a lot of downregulated genes. 

# Create a list of the gene sets
gene_sets <- list(
  "control vs acute" = DEGS$`control vs acute`$ID,
  "control vs intermediate" = DEGS$`control vs intermediate`$ID,
  "control vs prolonged" = DEGS$`control vs prolonged`$ID,
  "acute vs intermediary" = DEGS$`acute vs intermediate`$ID,
  "intermediary vs prolonged" = DEGS$`intermediate vs prolonged`$ID
)

# Generate the upset plot
a <- upset(fromList(gene_sets), 
           sets = c("control vs acute", "control vs intermediate",  "control vs prolonged" , 
                    "acute vs intermediary" ,"intermediary vs prolonged"),
           order.by = "freq",  
           main.bar.color = "#56B4E9",  
           sets.bar.color = "#009E73",
           text.scale = 1.5)

svg(paste0(out.dir, "/02_upSetR.svg"), w = 11, h = 4)
print(a)
dev.off()


# Volcano plots ----

# Load and filter all files
DGA <- lapply(DGA, function(x) {
  x <- x %>% mutate(deg_status = ifelse(FDR.P.val < fdr_treshold, "DEG", "Not DEG"))
  colnames(x)[4] <- "FC"
  return (x)
})

plot_volcano <- function(data, title, n_genes_to_be_labeled = 10) {
  top_degs <- data %>%
    filter(FDR.P.val < 0.05) %>%
    arrange(desc(abs(FC))) %>%
    slice_head(n = 10) %>%
    pull(ID)  # Extract gene IDs
  
  data <- data %>%
    mutate(Description = ifelse(ID %in% top_degs, ID, NA)) %>% 
    mutate(logFC = ifelse(FC > 0, log2(FC), -log2(-FC)))
  
  # Step 3: Create the plot
  ggplot(data, aes(x = logFC, y = -log10(P.val), colour = deg_status, label = Description)) +
    geom_point(alpha = 0.6) +
    geom_text_repel(aes(fontface=2))  +
    scale_alpha_manual(values = c(0.7, 1)) + 
    theme_minimal() + 
    scale_color_manual(values = c("DEG" = "black", "Not DEG" = "grey78")) +
    theme_minimal() +
    labs(title = title, x = "Log Fold Change", y = "-log10(P.val)") +
    theme(legend.position = "none")
}


for (name in names(DGA)) {
  plot_name <- file.path(out.dir, paste0("03_volcano_", gsub(" |-", "_", name), ".svg"))
  a <- plot_volcano(DGA[[name]], name)
  svg(plot_name) # Saves in disk
  print(a)
  dev.off() 
}


# P-value distribution across comparison ----

for (comp in names(DGA)){
  degs <- DGA[[comp]]
  
  a <- ggplot(degs, aes(x = P.val)) + 
    geom_histogram(fill = "lightblue", colour = "grey30", binwidth = 0.01, show.legend = FALSE) +
    geom_vline(xintercept = c(0.05), color = "red", linetype = "dashed", alpha = 0.5, linewidth = 0.5) +
    xlab(bquote("pvalue")) +
    ylab("Counts") +
    theme_bw() +
    theme(strip.text = element_text(size = 12)) + 
    ggtitle(comp) 
  
  svg(paste0(out.dir, "/04_pvalue_distribution", comp, ".svg"))
  print(a)
  dev.off()
  
}


# Save the filtered genes, both in .rds and .csv format----

# Save each comparison as CSV
for (comp in names(DGA)) {
  file_name <- gsub(" |-", "_", comp)
  file_path <- file.path(out.dir, paste0("filtered_", file_name, ".csv"))
  
  write.csv(DGA[[comp]], file = file_path, row.names = FALSE)
}

# Save the whole object (all comparisons) as RDS
saveRDS(DGA, file = file.path(out.dir, "DGA_filtered.rds"))

# This is the data used for the barplot of number of DEGs
write.csv(degs_summary_long,
          file = file.path(out.dir, "DEG_summary_counts.csv"),
          row.names = FALSE)
