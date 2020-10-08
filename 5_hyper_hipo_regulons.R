library(tidyverse)
library(limma)

#---- Loads all data requied to identify hipo and hyper methylated regulons ----

load("./rdata_files/Group3/4_dm_regulons.RData")

hm_regulons <- rownames(hm_regulons)

load("../expression_analysis/rdata_files/signature/g3_signature.RData")

tfs <- names(rtna@listOfRegulons)

load("./rdata_files/Group3/3_probewised.RData")


methylation_results <- decideTests(fit, p.value = 0.01)
methylation_results <- data.frame(methylation_results@.Data, 
                                  name = rownames(methylation_results@.Data),
                                  row.names = NULL)

methylation_results <- select(methylation_results, name, up_down = control.Group3)

methylation_results <- methylation_results %>% 
  inner_join(diff_methylated, by = "name") %>% 
  select(name, up_down, gene = ucsc_refgene_name) %>% 
  drop_na()


length(unique(methylation_results$gene))
length(intersect(unique(methylation_results$gene), probe_regulation$gene))

# Since there are more genes in the methylation_results data frame, in 
# comparison to probe_regulation, it is impossible to select the methylation
# probes by comparing them to the expression regulation results. The next step 
# might be to select all probes by their location on the gene, before doing
# the investigation, for both expression and methylation analyses.






# methylation_results <- methylation_results[methylation_results$gene %in% tfs, ]
# 
# #---- By summing up all results from fit model, define hipo and hyper 
# # methylated genes ----
# 
# results <- decideTests(fit, p.value = 0.01)
# results <- data.frame(results@.Data, name = rownames(results@.Data),
#                       row.names = NULL)
# 
# results <- select(results, name, up_down = control.Group3)
# 
# results <- results %>% 
#   inner_join(diff_methylated, by = "name") %>% 
#   select(name, up_down, ucsc_refgene_name) %>% 
#   drop_na()
# 
# results <- results %>% 
#   group_by(gene = ucsc_refgene_name) %>% 
#   transmute(sum = sum(up_down)) %>% 
#   distinct()
# 
# 
# ggplot(results) +
#   geom_tile()
# 
# #---- Now, define portion of hipo and hyper methylation in each regulon ----
# 
# function(regulon){
#   return(results[results$gene %in% regulon, ])
#   
# }
# 
# genes_in_regulon <- 
#   lapply(regulons_list, function(regulon) results[results$gene %in% regulon, ])
# 
# 
# gdata::keep(hm_regulons, regulons_list, dm_genes, sure = T)
# 
# sapply(regulons_list, function(x) sum(dm_genes %in% x))

       