#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @Description: Plot Enrichment results from ORSUM and other ontologies
# @software version: R=4.4.2

# Library
sh <- suppressPackageStartupMessages

sh(library(tidyverse))
sh(library(DOSE))

out.dir.results <- "results/05_enrichement_plots/"

if (!dir.exists(out.dir.results)){
  dir.create(out.dir.results)
}

load_and_process_orsum <- function(filepath, comparison_name) {
  df <- read.csv(filepath, sep = "\t")
  colnames(df)[6] <- "Comparison Rank"
  df %>% mutate(comparison = comparison_name)
  #df %>% slice_head(n = 15) 
}

# Define a helper function to simplify loading with direction label
load_orsum <- function(path, contrast, direction) {
  load_and_process_orsum(path, contrast) %>% mutate(direction = direction)
}


orsum_terms <- read.csv("results/03_enrichement/Summary_All_Contrasts/filteredResult-Summary.tsv", sep = "\t") %>% 
  filter(Representing.term.rank < 11) %>% 
  rename(repr_rank = `Representing.term.rank`)

orsum_terms <- orsum_terms %>%
  pivot_longer(
    cols = ends_with(".term.rank"),  # Select all rank columns
    names_to = "contrast",
    values_to = "rank"
  )

orsum_terms <- orsum_terms %>% filter(rank != "None")


# Add direction
orsum_terms$direction <- ifelse(grepl("^Up", orsum_terms$contrast), "Upregulated", "Downregulated")

#Replace name in the contrast
contrast_vector <- c("control vs intemediate", "intemediate vs prolonged",  "control vs acute", 
                     "acute vs intemediate", "control vs intemediate", "intemediate vs prolonged", 
                     "control vs prolonged")

names(contrast_vector) <- names(table(orsum_terms$contrast))


orsum_terms$contrast <- unname(contrast_vector[orsum_terms$contrast])

# Order
orsum_terms$contrast <- factor(orsum_terms$contrast, levels = c("control vs acute", "control vs intemediate", "acute vs intemediate",
                                                                "control vs prolonged",
                                                                "intemediate vs prolonged"))
orsum_terms$direction <- factor(orsum_terms$direction, levels = c("Upregulated", "Downregulated"))


orsum_terms <- orsum_terms %>% 
  arrange(contrast, direction)

orsum_terms$`Representing.term.name` <- factor(orsum_terms$`Representing.term.name`,
                                               levels = unique(orsum_terms$`Representing.term.name`)[length( unique(orsum_terms$`Representing.term.name`)):1])

#Filter out comparisons for plotting
orsum_terms<-orsum_terms%>%
  filter (contrast %in% c("control vs acute", "control vs intermediate", "control vs prolonged"))


## Plot (GO) Orsum data ----
pdf(paste0(out.dir.results, "GeneOntology_BP_Orsum.pdf"), w = 13, h = 10)
plot.go.bp.dir <- ggplot(orsum_terms, aes(x = contrast, y = `Representing.term.name`, colour = direction)) +
  geom_point(size = 4) +
  scale_fill_manual(values = c("Upregulated" = "#CC6677", "Downregulated" = "#88CCEE")) +
  theme_dose(10) +
  ylab("") +
  theme(axis.text.x = element_text(size = 18, angle = 90, hjust = 1, vjust = 0.1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        plot.margin = margin(l = 100), 
        legend.text = element_text(size = 18),
        legend.title = element_blank()) + 
  xlab("")
plot(plot.go.bp.dir)
dev.off()


# Save data for GO terms

saveRDS(
  orsum_terms,
  file = file.path(out.dir.results, "GO_BP_Orsum_terms.rds")
)



## Load KEGG/DO/DNG data----
load_enrichement_data <- function(filepath, comparison_name){
  df <- read.csv(filepath, sep = ",")
  df %>% mutate(comparison = comparison_name) %>% slice_head(n = 10)
}

## KEGG ----
kegg_terms_acute.up <- load_enrichement_data(
  "results/03_enrichement/control vs acute/kegg.up.results.csv",
  "control vs acute"
) %>% mutate(direction = "Up")


kegg_terms_intemediate.down <- load_enrichement_data(
  "results/03_enrichement/control vs intermediate/kegg.down.results.csv",
  "control vs intermediate"
) %>% mutate(direction = "Down")

kegg_terms_prolonged.up <- load_enrichement_data(
  "results/03_enrichement/control vs prolonged/kegg.up.results.csv",
  "control vs prolonged"
) %>% mutate(direction = "Up")

kegg_terms_intemediate_prolonged.up <- load_enrichement_data(
  "results/03_enrichement/intermediate vs prolonged/kegg.up.results.csv",
  "intermediate vs prolonged"
) %>% mutate(direction = "Up")

kegg_terms_intemediate_prolonged.down <- load_enrichement_data(
  "results/03_enrichement/intermediate vs prolonged/kegg.down.results.csv",
  "intermediate vs prolonged"
) %>% mutate(direction = "Down")

KEGG <- do.call("rbind", 
                list(kegg_terms_acute.up,
                     kegg_terms_intemediate.down,
                     kegg_terms_prolonged.up,
                     kegg_terms_intemediate_prolonged.up,
                     kegg_terms_intemediate_prolonged.down
                )
)


KEGG <- KEGG %>% 
  mutate(comparison = factor(comparison, levels = c("control vs acute", "control vs intermediate",
                                                     "acute vs intermediate", "control vs prolonged", "intermediate vs prolonged")),
         direction = factor(direction, levels = c("Up", "Down"))) %>% 
  arrange(comparison,direction) %>% 
  mutate(Description = factor(Description, levels = unique(Description)[length(unique(Description)):1]))

plot.KEGG.dir <- ggplot(KEGG, aes(x = comparison, y = Description, colour = direction)) +
  geom_point(size = 4) +
  scale_fill_manual(values = c("Upregulated" = "#CC6677", "Downregulated" = "#88CCEE"))  + 
  theme_dose(10) +
  ylab("") +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        plot.margin = margin(l = 100), 
        legend.text = element_text(size = 16),
        legend.title = element_blank()) + 
  xlab("")

pdf(paste0(out.dir.results, "KEGG.direction.pdf"), w = 15, h = 12)
plot(plot.KEGG.dir)
dev.off()
