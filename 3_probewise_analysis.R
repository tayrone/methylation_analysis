library(limma)
library(minfi)
library(gdata)
library(stringr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

subgroups <- c("SHH", "Group3", "Group4")

for(subgroup in subgroups){

  load(paste0("./rdata_files/", subgroup, "/2_filtered.RData"))
  
  targets <- as.data.frame(filtered_quantile_normalized@colData)
  
  annotation_epic <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  
  
  #---- Statistical analysis for further evaluation ----
  
  contrast <- factor(targets$Subgroup)
  
  design <- model.matrix(~0+contrast, data = targets)
  
  colnames(design) <- levels(contrast)
  
  m_values <- minfi::getM(filtered_quantile_normalized)
  
  fit <- lmFit(m_values, design)
  
  contrast_matrix <- 
    makeContrasts(paste0(colnames(design)[1], "-", colnames(design)[2]), 
                  levels = design)
  
  fit <- contrasts.fit(fit, contrast_matrix)
  fit <- eBayes(fit)
  
  summary(decideTests(fit, p.value = 0.01)) #let's peek at the amount of DM probes
  
  rm(contrast)
  
  
  #---- Assembles all information (hence number = "Inf"), 
  # ordered by p-value, for further inspection ---
  
  row_positions <- match(rownames(m_values), annotation_epic$Name)
  
  subbed_annotation_epic <- as.data.frame(annotation_epic[row_positions, ])
  
  all_cpgs <- limma::topTable(fit, coef = 1, number = Inf, 
                              sort.by = "p", genelist = subbed_annotation_epic)
  
  empty_columns <- which(sapply(all_cpgs, function(x) all(is.na(x) | x == '')))
  
  all_cpgs <- all_cpgs[, -empty_columns]
  
  colnames(all_cpgs) <- tolower(colnames(all_cpgs))
  
  all_cpgs$ucsc_refgene_name <- 
    sapply(strsplit(all_cpgs$ucsc_refgene_name, split = ";", fixed = T), `[`, 1)
  
  
  #---- First, selection of significant probes 
  # is done by logfc and ajusted p values.
  # Then, count the number of DM probes related to each gene, 
  # to check the level of methylation of each one ---- 
  
  diff_methylated <- 
    all_cpgs[abs(all_cpgs$logfc) >= 1.5 & all_cpgs$adj.p.val <= 0.01, ]
  
  genes_table <-  
    sort(table(as.character(diff_methylated$ucsc_refgene_name), useNA = "no"), 
         decreasing = T)
  
  dm_genes <- names(genes_table)
  
  
  save(diff_methylated, diff_methylated, fit, dm_genes, genes_table,
       file = paste0("./rdata_files/", subgroup, "/3_probewised.RData"))
  
  gdata::keep(subgroups, subgroup, sure = T)
  
}
