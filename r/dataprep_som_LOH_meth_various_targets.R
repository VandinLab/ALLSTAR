# Install packages --------------------------------------------------------
#install.packages("renv")

renv::init()
renv::install(c("base", "readxl", "MASS", "mgcv","dplyr", "tidyr",
                "cluster", "foreign", "nnet", "rpart", "Hmisc", "readr", "utils"))
renv::snapshot()

# Load packages -----------------------------------------------------------
library(utils)
library(base)
library(readxl)
library(dplyr)
library(tidyr)
library(Hmisc)
library(readr)


# Set threshold to filter most frequent mutations -------------------------
top_threshold <- 300


# Germline mutations ------------------------------------------------------
source("r/germline_load.R")
germline_mut <- germline_load()


# Somatic mutations -------------------------------------------------------
source("r/somatic_load.R")
somatic_mut <- somatic_load()

# create data frame with top frequent mutations
somatic_mut_sorted <- somatic_mut %>% select(-TCGA_ID) %>% colSums() %>% sort(decreasing = T, index.return = T)
somatic_mut_top <- somatic_mut %>% select(-TCGA_ID) %>% select(somatic_mut_sorted$ix[1:top_threshold])
somatic_mut_top <- cbind.data.frame(somatic_mut$TCGA_ID, somatic_mut_top)
colnames(somatic_mut_top)[1] <- "TCGA_ID"


# LOH information ------------------------------------------------------
source("r/LOH_load.R")
LOH <- LOH_load()

# filter for top LOH genes
LOH_sorted <- LOH %>% select(-TCGA_ID) %>% colSums(na.rm = T) %>% sort(decreasing = T, index.return = T)
LOH_top <- LOH %>% select(-TCGA_ID) %>% select(LOH_sorted$ix[1:top_threshold])
LOH_top <- cbind.data.frame(LOH$TCGA_ID, LOH_top)
colnames(LOH_top)[1] <- "TCGA_ID"

# remove rows with one NA
LOH_top <- LOH_top[complete.cases(LOH_top), ]


# Methylation data ------------------------------------------------------
source("r/methylation_load.R")
methylation <- methylation_load()


# Aggregated mutation data ------------------------------------------------
source("r/aggregated_mut_load.R")
aggregated_mut <- aggregated_mut_load()

# create data frame with top frequent mutations
aggregated_mut_sorted <- aggregated_mut %>% select(-TCGA_ID) %>% colSums() %>% sort(decreasing = T, index.return = T)
aggregated_mut_top <- aggregated_mut %>% select(-TCGA_ID) %>% select(aggregated_mut_sorted$ix[1:top_threshold])
aggregated_mut_top <- cbind.data.frame(aggregated_mut$TCGA_ID, aggregated_mut_top)
colnames(aggregated_mut_top)[1] <- "TCGA_ID"

# Clinical ----------------------------------------------------------------
source("r/clinical_load.R")
clinical_factor <- clinical_load()


# Molecular subtype -------------------------------------------------------
source("r/molecular_subtype_load.R")
molecular_subtype <- molecular_subtype_load()


# Triple negative ---------------------------------------------------------
source("r/triple_negative_load.R")
triple_negative <- triple_negative_load()


# LVI ---------------------------------------------------------------------
LVI <- clinical_factor %>%
  select(TCGA_ID, LVI)


# Histological type -------------------------------------------------------
histological_type <- clinical_factor %>%
  select(TCGA_ID, histological_type)


# Covariates --------------------------------------------------------------
selected_clinical <- clinical_factor %>%
  select(TCGA_ID, gender, race, history_other_malignancy, age_at_diagnosis, menopause_status) %>%
  left_join(., germline_mut, by = "TCGA_ID")
selected_clinical[is.na.data.frame(selected_clinical)] <- 0


# Aggregate: molecular subtype --------------------------------------------
### Just somatic mutations
data_somatic <- inner_join(selected_clinical, somatic_mut_top, by = "TCGA_ID")
data_somatic_subtype <- inner_join(data_somatic, molecular_subtype, by = "TCGA_ID")

`%!in%` = Negate(`%in%`)
data_somatic_subtype <- data_somatic_subtype %>%
  filter(molecular_subtype %!in% c("CLOW", "Not available") & !is.na(molecular_subtype))

data_somatic_subtype$molecular_subtype <- as.numeric(factor(data_somatic_subtype$molecular_subtype))

write.csv(data_somatic_subtype %>% 
            select(-TCGA_ID) , "inputs/som_subtype.csv", 
          quote = F, sep = ",", row.names = F)


### Somatic, LOH and mutations
data_somatic <- inner_join(selected_clinical, somatic_mut_top, by = "TCGA_ID")
data_somatic_LOH <- inner_join(data_somatic, LOH_top, by = "TCGA_ID")
data_somatic_LOH_meth <- inner_join(data_somatic_LOH, methylation, by = "TCGA_ID")
data_somatic_LOH_meth_subtype <- inner_join(data_somatic_LOH_meth, molecular_subtype, by = "TCGA_ID")

`%!in%` = Negate(`%in%`)
data_somatic_LOH_meth_subtype <- data_somatic_LOH_meth_subtype %>%
  filter(molecular_subtype %!in% c("CLOW", "Not available") & !is.na(molecular_subtype))

data_somatic_LOH_meth_subtype$molecular_subtype <- as.numeric(factor(data_somatic_LOH_meth_subtype$molecular_subtype))

write.csv(data_somatic_LOH_meth_subtype %>% 
            select(-TCGA_ID) , "inputs/som_LOH_meth_subtype.csv", 
          quote = F, sep = ",", row.names = F)


### Aggregated
data_aggregated_mut <- inner_join(selected_clinical, aggregated_mut_top, by = "TCGA_ID")
data_aggregated_mut_subtype <- inner_join(data_aggregated_mut, molecular_subtype, by = "TCGA_ID")

`%!in%` = Negate(`%in%`)
data_aggregated_mut_subtype <- data_aggregated_mut_subtype %>%
  filter(molecular_subtype %!in% c("CLOW", "Not available") & !is.na(molecular_subtype))

data_aggregated_mut_subtype$molecular_subtype <- as.numeric(factor(data_aggregated_mut_subtype$molecular_subtype))

write.csv(data_aggregated_mut_subtype %>% 
            select(-TCGA_ID) , "inputs/aggregated_mut_subtype.csv", 
          quote = F, sep = ",", row.names = F)


# Aggregate: LVI ----------------------------------------------------------
### Just somatic mutations
data_somatic <- inner_join(selected_clinical, somatic_mut_top, by = "TCGA_ID")
data_somatic_LVI <- inner_join(data_somatic, LVI, by = "TCGA_ID")

write.csv(data_somatic_LVI %>% 
            select(-TCGA_ID) , "inputs/som_LVI.csv", 
          quote = F, sep = ",", row.names = F)


### Somatic, LOH and mutations
data_somatic <- inner_join(selected_clinical, somatic_mut_top, by = "TCGA_ID")
data_somatic_LOH <- inner_join(data_somatic, LOH_top, by = "TCGA_ID")
data_somatic_LOH_meth <- inner_join(data_somatic_LOH, methylation, by = "TCGA_ID")
data_somatic_LOH_meth_LVI <- inner_join(data_somatic_LOH_meth, LVI, by = "TCGA_ID")

write.csv(data_somatic_LOH_meth_LVI %>% 
            select(-TCGA_ID) , "inputs/som_LOH_meth_LVI.csv", 
          quote = F, sep = ",", row.names = F)

### Aggregated
data_aggregated_mut <- inner_join(selected_clinical, aggregated_mut_top, by = "TCGA_ID")
data_aggregated_mut_LVI <- inner_join(data_aggregated_mut, LVI, by = "TCGA_ID")

write.csv(data_aggregated_mut_LVI %>% 
            select(-TCGA_ID) , "inputs/aggregated_mut_LVI.csv", 
          quote = F, sep = ",", row.names = F)



# Aggregate: triple negative ----------------------------------------------

### Just somatic mutations
data_somatic <- inner_join(selected_clinical, somatic_mut_top, by = "TCGA_ID")
data_somatic_3N <- inner_join(data_somatic, triple_negative, by = "TCGA_ID")

write.csv(data_somatic_3N %>% 
            select(-TCGA_ID) , "inputs/som_3N.csv", 
          quote = F, sep = ",", row.names = F)


### Somatic, LOH and mutations
data_somatic <- inner_join(selected_clinical, somatic_mut_top, by = "TCGA_ID")
data_somatic_LOH <- inner_join(data_somatic, LOH_top, by = "TCGA_ID")
data_somatic_LOH_meth <- inner_join(data_somatic_LOH, methylation, by = "TCGA_ID")
data_somatic_LOH_meth_3N <- inner_join(data_somatic_LOH_meth, triple_negative, by = "TCGA_ID")

write.csv(data_somatic_LOH_meth_3N %>% 
            select(-TCGA_ID) , "inputs/som_LOH_meth_3N.csv", 
          quote = F, sep = ",", row.names = F)

### Aggregated
data_aggregated_mut <- inner_join(selected_clinical, aggregated_mut_top, by = "TCGA_ID")
data_aggregated_mut_3N <- inner_join(data_aggregated_mut, triple_negative, by = "TCGA_ID")

write.csv(data_aggregated_mut_3N %>% 
            select(-TCGA_ID) , "inputs/aggregated_mut_3N.csv", 
          quote = F, sep = ",", row.names = F)


# Aggregate: histologic type -------------------------------------------
### Just somatic mutations
data_somatic <- inner_join(selected_clinical, somatic_mut_top, by = "TCGA_ID")
data_somatic_histo <- inner_join(data_somatic, histological_type, by = "TCGA_ID")

write.csv(data_somatic_histo %>% 
            select(-TCGA_ID) , "inputs/som_histo.csv", 
          quote = F, sep = ",", row.names = F)


### Somatic, LOH and mutations
data_somatic <- inner_join(selected_clinical, somatic_mut_top, by = "TCGA_ID")
data_somatic_LOH <- inner_join(data_somatic, LOH_top, by = "TCGA_ID")
data_somatic_LOH_meth <- inner_join(data_somatic_LOH, methylation, by = "TCGA_ID")
data_somatic_LOH_meth_histo <- inner_join(data_somatic_LOH_meth, histological_type, by = "TCGA_ID")

write.csv(data_somatic_LOH_meth_histo %>% 
            select(-TCGA_ID) , "inputs/som_LOH_meth_histo.csv", 
          quote = F, sep = ",", row.names = F)

### Aggregated
data_aggregated_mut <- inner_join(selected_clinical, aggregated_mut_top, by = "TCGA_ID")
data_aggregated_mut_histo <- inner_join(data_aggregated_mut, histological_type, by = "TCGA_ID")

write.csv(data_aggregated_mut_histo %>% 
            select(-TCGA_ID) , "inputs/aggregated_mut_histo.csv", 
          quote = F, sep = ",", row.names = F)
