# Loads data of contrast groups under investigation. 
library(GEOquery)
library(tidyverse)
library(minfi)
library(gdata)


#---- Creates control targets object, which contains clinical data for 
# each control sample ----

control_files <- list.files("control_data", "GSM*", full.names = T)
control_files <- str_extract(control_files, "^.*[0-9]")

trimmed_filenames <- str_extract(control_files, "GSM[0-9]*")

trimmed_filenames <- unique(trimmed_filenames)

control_targets <- sapply(trimmed_filenames, getGEO)

control_targets <- 
  sapply(control_targets, function(x) x@header$characteristics_ch1)

control_targets <- t(control_targets)
#colnames(control_targets) <- c("tissue", "age", "sex", "subgroup")

# Removes unnecessary info from data frame cells
for(i in row.names(control_targets)){
  control_targets[i, ] <- sapply(strsplit(control_targets[i, ], 
                                          split = ": ", fixed = T), `[`, 2)
}


control_targets <- transmute(as.data.frame(control_targets),
                             Gender = as.character(V3), #define as not a factor
                             Age = tolower(V2),
                             Basename = control_files[seq(1, 30, 2)],
                             Subgroup = "control")

control_targets$Gender[control_targets$Gender == "UNKNOWN"] <- NA_character_


raw_control_data <- read.metharray.exp(targets = control_targets, force = T)


#---- Cancer-related targets object is built, in order to complement
# control_targets ----

case_data_directory <- "./cancer_data"

case_targets <- read.metharray.sheet(case_data_directory, 
                                     pattern = "targets.csv")

rownames(case_targets) <- case_targets$Acession

#case_targets$Contrast <- "g34"

case_targets <- select(case_targets, Gender, Age, Basename, Subgroup = Type)

case_targets <- case_targets[!(case_targets$Subgroup == "WNT") &
                              (case_targets$Basename != "character(0)"), ]

raw_case_data <- read.metharray.exp(targets = case_targets)

gdata::keep(case_targets, control_targets, 
            raw_control_data, raw_case_data, sure = T)


save(case_targets, control_targets, raw_case_data, raw_control_data,
     file = paste0("./rdata_files/1_loaded.RData"))
