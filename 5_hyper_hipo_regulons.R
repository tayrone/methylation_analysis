load("/data4/tayrone25/medulloblastoma/methylation/rdata_files/Group3/4_dm_regulons.RData")

hm_regulons <- rownames(hm_regulons)

load("/data4/tayrone25/medulloblastoma/expression_analysis/rdata_files/network/g3_rtn.RData")

regulons_list <- rtna@listOfRegulons


load("/data4/tayrone25/medulloblastoma/methylation/rdata_files/Group3/3_probewised.RData")

dm_genes

gdata::keep(hm_regulons, regulons_list, dm_genes, sure = T)

 <- sapply(regulons_list, function(x) sum(dm_genes %in% x))

       