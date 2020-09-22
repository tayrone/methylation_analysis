library(tidyverse)

#---- Loads all data requied to identify hipo and hyper methylated regulons ----

load("/data4/tayrone25/medulloblastoma/methylation/rdata_files/
     Group3/4_dm_regulons.RData")

hm_regulons <- rownames(hm_regulons)

load("/data4/tayrone25/medulloblastoma/expression_analysis/
     rdata_files/network/g3_rtn.RData")

regulons_list <- rtna@listOfRegulons


load("/data4/tayrone25/medulloblastoma/methylation/rdata_files/
     Group3/3_probewised.RData")


#---- By summing up all results from fit model, define hipo and hyper 
# methylated genes ----

results <- decideTests(fit, p.value = 0.01)
results <- data.frame(results@.Data, name = rownames(results@.Data),
                      row.names = NULL)

results <- select(results, name, up_down = control.Group3)

results <- results %>% 
  inner_join(diff_methylated, by = "name") %>% 
  select(name, up_down, ucsc_refgene_name) %>% 
  drop_na()

results <- results %>% 
  group_by(gene = ucsc_refgene_name) %>% 
  transmute(sum = sum(up_down)) %>% 
  distinct()


#---- Now, define portion of hipo and hyper methylation in each regulon ----

gdata::keep(hm_regulons, regulons_list, dm_genes, sure = T)

sapply(regulons_list, function(x) sum(dm_genes %in% x))

       