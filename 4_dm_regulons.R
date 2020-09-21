# This script aims to identify highly methylated regulons among
# the ones identified as master regulators.

options(stringsAsFactors = F)
library(tidyverse)
library(gdata)
library(ggplot2)
library(RTN)

#---- Loading and preprocessing data ----

subgroups <- c("Group3", "Group4", "SHH")
names(subgroups) <- c("g3", "g4", "shh")

for(i in 1:length(subgroups)){
  print(subgroups[i])
  print(names(subgroups)[i])


  load(paste0("./rdata_files/", subgroups[i], "/3_probewised.RData"))
  
  load(paste0("../expression_analysis/rdata_files/network/", 
              names(subgroups)[i],"_rtn.RData"))
  
  rm(rtni)
  
  regulons <- tna.get(rtna, what = "regulons")
  
  tfs <- rtna@regulatoryElements
  
  #dm_genes <- gene_coords$hgnc_symbol
  
  
  #---- Table creation for each regulon ----
  
  # Indexing converts all unexistent values to NA, instead of repeating the vector
  #regulons <- sapply(regulons, '[', seq(max(lengths(regulons))))
  #regulons <- as.data.frame(regulons)
  
  human_genes_count <- length(rtna@phenotype)
  
  cont_table <- data.frame(matrix(0, nrow = 2, ncol = 2),
                           row.names = c("in_regulon", "out_regulon"))
  
  colnames(cont_table) = c("met", "not_met")
  
  tables_creation <- function(x){
    
    unlist(x)
    
    cont_table["in_regulon", "met"] <- sum(dm_genes %in% x)
    
    
    cont_table["in_regulon", "not_met"] <- sum(!is.na(x)) - sum(dm_genes %in% x)
    
    cont_table["out_regulon", "met"] <- sum(!(dm_genes %in% x))
    
    cont_table["out_regulon", "not_met"] <- 
      human_genes_count - sum(cont_table)
    
    return(cont_table)
    
  }
  
  #---- Significance test ----
  
  tables <- lapply(regulons, tables_creation)
  
  p_values <- lapply(tables, fisher.test)
  
  p_values <- sapply(p_values, function(x) return(x$p.value))
  
  adjusted_p <- p.adjust(p_values, method = "BH")
  
  hm_regulons <- as.data.frame(adjusted_p)

  hm_regulons <- hm_regulons %>% 
    arrange(adjusted_p) %>% 
    filter(adjusted_p < 0.05)
  
  
  #---- Writing results ----
  
  mrs <- tna.get(rtna, what = "mra")
  
  mr_and_hm <- intersect(rownames(hm_regulons), mrs$Regulon)
  
  write.csv(mr_and_hm, paste0("./rdata_files/", subgroups[i], "/mr_and_hm.csv"))
  write.csv(hm_regulons, paste0("./rdata_files/", subgroups[i], "/hm_regulons.csv"))
  
  save(tables, p_values, adjusted_p, hm_regulons, #mr_and_hm,
       file = paste0("./rdata_files/", subgroups[i], "/4_dm_regulons.RData"))
  
  #gdata::keep(subgroups, sure = T)

}
