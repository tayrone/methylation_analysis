library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(gdata)
library(maxprobes)

load("rdata_files/1_loaded.RData")

#subgroups <- c("SHH", "Group3", "Group4")
#case_targets <- case_targets[case_targets$Subgroup = subgroup, ]


rg_set <- combineArrays(raw_case_data, raw_control_data, 
                        outType = "IlluminaHumanMethylationEPIC")

gdata::keep(rg_set, sure = T)

#---- Detection p-values are calculated for both case and control data ----


p_values <- detectionP(rg_set)


p_colmeans <- data.frame(col_means = colMeans(p_values), 
                         samples = colnames(p_values))

range(p_colmeans$col_means)

quantile_normalized <- preprocessQuantile(rg_set)


#---- Removes probes that are flawed for one or more samples ----

p_values <- p_values[match(featureNames(quantile_normalized), 
                                             rownames(p_values)), ] 

keep <- rowSums(p_values < 0.01) == ncol(quantile_normalized) 
table(keep)

filtered_quantile_normalized <- quantile_normalized[keep, ]

rm(p_values, keep)


#---- Sex influences DNA methylation, so we remove probes on sex chromosomes ----

annotation_epic <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

sex_chromo_probes <- 
  annotation_epic$Name[annotation_epic$chr %in% c("chrX","chrY")]

keep <- !(featureNames(filtered_quantile_normalized) %in% sex_chromo_probes)

table(keep)

filtered_quantile_normalized <- filtered_quantile_normalized[keep, ]

rm(sex_chromo_probes, keep)


#---- Removal of probes where common SNPs may affect CpGs ----

filtered_quantile_normalized <- 
  minfi::dropLociWithSnps(filtered_quantile_normalized)


#---- Cross reactive probes can mistakenly infer some sex-associated methylation ----

cross_reactive_probes <- unlist(xreactive_probes(array_type = "EPIC"))

keep <- !(featureNames(filtered_quantile_normalized) %in% 
            cross_reactive_probes)

table(keep)

filtered_quantile_normalized <- filtered_quantile_normalized[keep, ]

rm(cross_reactive_probes, keep)

save.image("./rdata_files/2_filtering.RData")

