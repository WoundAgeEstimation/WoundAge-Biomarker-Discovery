#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @Description: Perform enrichment analysis on the data. Run using Up + Down in separate and together. Use several ontologies
# @software version: R=4.4.2

# Library
sh <- suppressPackageStartupMessages

sh(library(clusterProfiler))
sh(library(org.Hs.eg.db))
sh(library(tidyverse))
sh(library(DOSE))


#  Make data dir 
results.dir <- "results/"
out.dir.results <- paste0(results.dir, "03_enrichement/")

if (!dir.exists(out.dir.results)) {dir.create(out.dir.results)}

# Load results ----

degs_data <- readRDS(file.path(results.dir, "02_degs", "DGA_filtered.rds"))

# Helper save function

save_results <- function(go_obj, prefix, out_dir, sum_results = T) {
  if (is.null(go_obj)){
    return ()
  }
  
  results <- go_obj@result
  sig_ids <- results %>% filter(p.adjust < 0.05) %>% pull(ID)
  
  if (length(sig_ids) > 0) {
    saveRDS(go_obj, file.path(out_dir, paste0(prefix, ".results.rds")))
    write.csv(results, file.path(out_dir, paste0(prefix, ".results.csv")), row.names = FALSE)
    
    if (sum_results == T){
      write.table(sig_ids, file.path(out_dir, paste0("orsum.", prefix, ".results.csv")),
                  quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
    
  }
}

# Enrichment function 

perform_enrichement <- function(degs, name){
  degs.up <- degs %>% 
    filter(FDR.P.val < 0.05 & FC > 1) %>% 
    pull(ID)
  
  degs.down <- degs %>% 
    filter(FDR.P.val < 0.05 & FC < 1) %>% 
    pull(ID)
  
  degs.all <- c(degs.up, degs.down)
  
  out.dir <- file.path(out.dir.results, name)
  out.dir <- paste0(out.dir, "/")
  if (!dir.exists(out.dir)) {dir.create(out.dir)}
  
  
  degs_id <- bitr(degs$ID, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  

  # Enrichment in UP
  if (length(degs.up) != 0){
    go.bp.up <- enrichGO(gene = degs.up, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "SYMBOL", universe = degs$ID, minGSSize = 5)
    go.mf.up <- enrichGO(gene = degs.up, OrgDb = org.Hs.eg.db, ont = "MF", keyType = "SYMBOL", universe = degs$ID, minGSSize = 5)
    go.cc.up <- enrichGO(gene = degs.up, OrgDb = org.Hs.eg.db, ont = "CC", keyType = "SYMBOL", universe = degs$ID, minGSSize = 5)
    
    save_results(go.bp.up, "go.bp.up", out.dir)
    save_results(go.mf.up, "go.mf.up", out.dir)
    save_results(go.cc.up, "go.cc.up", out.dir)
    
    
  }
  # Enrichment in UP (entrez)
  degs.up.entrz <- bitr(degs.up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (length(degs.up.entrz) != 0){
    kegg.up <- enrichKEGG(gene = degs.up.entrz$ENTREZID, universe = degs_id$ENTREZID, minGSSize = 5)
    do.up <- enrichDO(gene = degs.up.entrz$ENTREZID, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = degs_id$ENTREZID)
    dgn.up <- enrichDGN(gene = degs.up.entrz$ENTREZID, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = degs_id$ENTREZID)
    
    save_results(kegg.up, "kegg.up", out.dir, sum_results = F)
    save_results(do.up, "do.up", out.dir, sum_results = F)
    save_results(dgn.up, "dgn.up", out.dir, sum_results = F)
    
  }
  
  # Enrichment in Down
  if (length(degs.down) != 0){
    go.bp.down <- enrichGO(gene = degs.down, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "SYMBOL", universe = degs$ID, minGSSize = 5)
    go.mf.down <- enrichGO(gene = degs.down, OrgDb = org.Hs.eg.db, ont = "MF", keyType = "SYMBOL", universe = degs$ID, minGSSize = 5)
    go.cc.down <- enrichGO(gene = degs.down, OrgDb = org.Hs.eg.db, ont = "CC", keyType = "SYMBOL", universe = degs$ID, minGSSize = 5)
    
    save_results(go.bp.down, "go.bp.down", out.dir)
    save_results(go.mf.down, "go.mf.down", out.dir)
    save_results(go.cc.down, "go.cc.down", out.dir)
    
  }
  # Enrichment in Down (entrez)
  degs.down.entrz <- bitr(degs.down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (length(degs.down.entrz) != 0){
    kegg.down <- enrichKEGG(gene = degs.down.entrz$ENTREZID, universe = degs_id$ENTREZID, minGSSize = 5)
    do.down <- enrichDO(gene = degs.down.entrz$ENTREZID, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = degs_id$ENTREZID)
    dgn.down <- enrichDGN(gene = degs.down.entrz$ENTREZID, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = degs_id$ENTREZID)
    
    save_results(kegg.down, "kegg.down", out.dir, sum_results = F)
    save_results(do.down, "do.down", out.dir, sum_results = F)
    save_results(dgn.down, "dgn.down", out.dir, sum_results = F)
    
  }
  
  
  # Enrichment in all
  if (length(degs.all) != 0){
    go.bp.all <- enrichGO(gene = degs.all, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "SYMBOL", universe = degs$ID, minGSSize = 5)
    go.mf.all <- enrichGO(gene = degs.all, OrgDb = org.Hs.eg.db, ont = "MF", keyType = "SYMBOL", universe = degs$ID, minGSSize = 5)
    go.cc.all <- enrichGO(gene = degs.all, OrgDb = org.Hs.eg.db, ont = "CC", keyType = "SYMBOL", universe = degs$ID, minGSSize = 5)
    
    save_results(go.bp.all, "go.bp.all", out.dir)
    save_results(go.mf.all, "go.mf.all", out.dir)
    save_results(go.cc.all, "go.cc.all", out.dir)
    
  }
  # Enrichment in all
  degs.all.entrz <- bitr(degs.all, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (length(degs.all.entrz) != 0){
    kegg.all <- enrichKEGG(gene = degs.all.entrz$ENTREZID, universe = degs_id$ENTREZID, minGSSize = 5)
    do.all <- enrichDO(gene = degs.all.entrz$ENTREZID, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = degs_id$ENTREZID)
    dgn.all <- enrichDGN(gene = degs.all.entrz$ENTREZID, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = degs_id$ENTREZID)
    
    save_results(kegg.all, "kegg.all", out.dir, sum_results = F)
    save_results(do.all, "do.all", out.dir, sum_results = F)
    save_results(dgn.all, "dgn.all", out.dir, sum_results = F)
    
  }
}



for (i in names(degs_data)){
  name <- i 
  data_degs <- degs_data[[name]]
  perform_enrichement(data_degs, name)
}

