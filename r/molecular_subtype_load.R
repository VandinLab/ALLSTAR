molecular_subtype_load <- function() {
  
  molecular_annotation <- read_excel("data/TCGA_BRCA_histologic_features.xlsx", na = "NA")
  molecular_subtype <- molecular_annotation %>% select(`Sample CLID`, `PAM50 and Claudin-low (CLOW) Molecular Subtype`)
  colnames(molecular_subtype) <- c("TCGA_ID", "molecular_subtype")
  
  return(molecular_subtype)
  
}