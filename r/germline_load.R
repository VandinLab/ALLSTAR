germline_load <- function() {
  
  BRCA_germline_somatic_classification <- read_excel("data/BRCA_germline_somatic_classification.xlsx", 
                                                     skip = 1)
  # filter germline mutations only
  germline_mut <- BRCA_germline_somatic_classification %>% filter(`Alteration Type` == "germline mutation")
  # remove Penn data
  germline_mut <- germline_mut %>% filter(Site != "Penn")
  # filter important info only
  germline_mut <- germline_mut %>% select(`Tumor ID`, `Gene Altered`)
  # change colnames
  colnames(germline_mut) <- c("TCGA_ID", "gene")
  # add germline identifier
  germline_mut$gene <- paste(germline_mut$gene, "_germline", sep="")
  # spread info on gene level
  germline_mut$mut_presence <- 1
  germline_mut <- germline_mut %>% spread(gene, mut_presence, fill = 0)
  
  # select only data rich mutations
  germline_mut <- germline_mut %>%
    select(TCGA_ID, BRCA1_germline, BRCA2_germline)
  
  return(germline_mut)
  
}