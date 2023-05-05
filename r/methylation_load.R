methylation_load <- function() {
  methylation <- read_table("data/Methylation_Analysis_Data.txt")
  
  # filter important info only
  methylation <- methylation %>%
    select(Patient_ID, Gene, Hypermethylated)
  methylation <- methylation[!duplicated(methylation), ]
  
  # change colnames
  colnames(methylation) <- c("TCGA_ID", "gene", "hypermethylation")
  # add germline identifier
  methylation$gene <- paste(methylation$gene, "_meth", sep="")
  # spread info on gene level
  methylation$hypermethylation <- ifelse(methylation$hypermethylation == TRUE, 1, 0)
  methylation <- methylation %>% spread(gene, hypermethylation, fill = 0)
  
  return(methylation)
}